#' Poisson process factorization for topographic-dependent mutational signatures analysis
#'
#' @param gr_Mutations GenomicRanges object contatining covariates and samples
#' @param gr_SignalTrack Matrix containing signals for the whole genome
#' @param method Which method to use for the estimation. Available options are 'mle' and 'map'
#' @param bin_weight Weights to give to each region
#' @param K Upper bound to the number of signatures expected in the data
#' @param prior_params The hyperparameters of the prior.
#' @param init Initialization parameters. See \code{SigPoisProcess.init}
#' @param controls Control parameters for each method. See \code{SigPoisProcess.contol}
#'
#' @returns An object of class SigPoisProcess_mult
#' @export
#' @useDynLib SigPoisProcess, .registration = TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom magrittr %>%
SigPoisProcess_mult <- function(gr_Mutations,
                                SignalTrack,
                                bin_weight,
                                K = 10,
                                method = "map",
                                prior_params = SigPoisProcess_mult.PriorParams(),
                                controls = SigPoisProcess_mult.controls(),
                                init = SigPoisProcess_mult.init(),
                                update_Signatures = TRUE,
                                update_Theta = TRUE,
                                update_Betas = TRUE){

  # Extract the data and the covariates
  X <- as.data.frame(GenomicRanges::mcols(gr_Mutations)) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()

  # Check that colnames of SignalTrack and colnames of X match
  if(any(colnames(SignalTrack) != colnames(X))){
    stop("Must have colnames(SignalTrack) == colnames(X))")
  }

  # Extract the matrix of aggregated mutations
  MutMatrix <- getTotalMutations(gr_Mutations)

  # Extract mutation data
  Mutations <- as.data.frame(GenomicRanges::mcols(gr_Mutations)) %>%
    dplyr::select(sample, channel)
  Mutations$sample <- as.factor(Mutations$sample) #, levels = unique(Mutations$sample))
  Mutations$channel <- as.factor(Mutations$channel)

  # Extract channel ID and mutation ID
  channel_id <- as.numeric(Mutations$channel)
  sample_id <- as.numeric(Mutations$sample)

  #----------------------------------------------- Model dimensions
  I <- length(levels(Mutations$channel)) # Number of mutational channels
  J <- length(unique(sample_id)) # Number of patients
  p <- ncol(SignalTrack) # Number of covariates (excluding the intercept)
  N <- nrow(Mutations) # Total number of mutations

  ## ----------------------------------------------- Model settings
  prior_params <- do.call("SigPoisProcess_mult.PriorParams", prior_params)
  a <- prior_params$a
  alpha <- prior_params$alpha
  epsilon <- prior_params$epsilon
  if (is.null(prior_params$a0) | is.null(prior_params$b0)) {
    a0 <- a * J + 1
    b0 <- epsilon * a * J
  } else {
    a0 <- prior_params$a0
    b0 <- prior_params$b0
  }

  #---- Shrinkage method
  shrinkage = prior_params$shrinkage
  if (is.null(prior_params$c0) | is.null(prior_params$d0)) {
    c0 <- K / 2 + 1
    d0 <- epsilon * K / 2
    } else {
    c0 <- prior_params$c0
    d0 <- prior_params$d0
  }

  #----------------------------------------------- Prior for the signatures
  SigPrior <- prior_params$SigPrior
  if (!is.null(SigPrior)) {
    if (!is.matrix(SigPrior)) {
      stop("SigPrior must be a matrix")
    }
    if (length(setdiff(rownames(SigPrior), levels(Mutations$channel))) > 0) {
      stop("Must have rownames(SigPrior) equal to levels(Mutations$channel)")
    }
    Kpre <- ncol(SigPrior)
    SigPrior <- SigPrior[levels(Mutations$channel), ] # reorder
    if (K > 0) {
      SigPrior_new <- matrix(alpha, nrow = I, ncol = K)
      colnames(SigPrior_new) <- paste0("SigN", sprintf("%02d", 1:K))
      rownames(SigPrior_new) <- rownames(SigPrior)
      SigPrior <- cbind(SigPrior, SigPrior_new)
      K <- Kpre + K
    } else {
      K <- Kpre
    }
  } else {
    SigPrior <- matrix(alpha, nrow = I, ncol = K)
    colnames(SigPrior) <- paste0("SigN", sprintf("%02d", 1:K))
    rownames(SigPrior) <- levels(Mutations$channel)
  }

  #----------------------------------------------- Initial values
  # Initialize the model
  init <- do.call("SigPoisProcess_mult.init", init)
  #---- Initial value for the signatures
  if(is.null(init$R_start)){
    R_start <- sample_signatures(SigPrior)
  } else {
    R_start <- init$R_start
    if(ncol(R_start) != K){
      stop("ncol(R_start) must be equal to specified K")
    }
  }
  #---- Initial value for the baseline loadings
  if(is.null(init$Theta_start)){
    Theta_start <- matrix(rgamma(J * K, colSums(MutMatrix)), nrow = K, byrow = TRUE)/sum(bin_weight)
  } else {
    Theta_start <- init$Theta_start
  }
  #---- Initial value for the relevance weights
  if(is.null(init$Mu_start)){
    Mu_start <- rowMeans(Theta_start)#1 / rgamma(K, a * J  + 1, epsilon * a * J)
  } else {
    Mu_start <- init$Mu_start
  }

  #---- Initial value for the regression coefficients (start at 0)
  if(is.null(init$Betas_start)){
    Betas_start <- matrix(rnorm(p * K, mean = 0, sd = 1e-6), nrow = p, ncol = K)#matrix(0, nrow = p, ncol = K)
  } else {
    Betas_start <- init$Betas_start
  }

  #----------------------------------------------- Control parameters for the model
  if (!method %in% c("mle", "map")) {
    stop("method must be either 'map' or 'mle'")
  }
  controls <- do.call("SigPoisProcess_mult.controls", controls)

  #----------------------------------------------- Estimate the model
  # List to append results
  results <- list(sol = NULL, method = method)

  ######################################################
  # Optimization
  ######################################################
  cat("Start optimization \n")
  maxiter <- controls$maxiter
  n_iter_betas <- controls$n_iter_betas
  tol <- controls$tol
  # Find the MLE using the multiplicative rules
  out <- .PoissonProcess_optim(R_start = R_start,
                               Theta_start = Theta_start,
                               Betas_start = Betas_start,
                               X = X,
                               SignalTrack = SignalTrack,
                               bin_weight = bin_weight,
                               channel_id = channel_id - 1,
                               sample_id = sample_id - 1,
                               SigPrior = SigPrior,
                               a = a,
                               a0 = a0,
                               b0 = b0,
                               c0 = c0,
                               d0 = d0,
                               shrinkage = shrinkage,
                               update_R = update_Signatures,
                               update_Theta = update_Theta,
                               update_Betas = update_Betas,
                               method = method,
                               maxiter = maxiter,
                               n_iter_betas = n_iter_betas,
                               tol = tol)

  #--- Postprocess results
  # Signatures
  Sigs <- out$R
  colnames(Sigs) <- colnames(SigPrior)
  rownames(Sigs) <- rownames(SigPrior)
  # Loadings
  Thetas <- out$Theta
  rownames(Thetas) <- colnames(SigPrior)
  colnames(Thetas) <- colnames(MutMatrix)
  # Coefficients
  Betas <- out$Betas
  dimnames(Betas) <- list(colnames(X), colnames(SigPrior))

  #---- Store results in a list
  results$sol <- out
  results$Thetas <- Thetas
  results$Signatures <- Sigs
  results$Betas <- Betas

  if(method == "map"){
    results$Mu <- out$Mu
  } else {
    results$Mu <- NULL
  }

  results$init <- list(R_start = R_start,
                       Theta_start = Theta_start,
                       Betas_start = Betas_start,
                       Mu_start = Mu_start)

  results$PriorParams <- list(a = a, alpha = alpha, epsilon = epsilon,
                         a0 = a0, b0 = b0, SigPrior = SigPrior)


  class(results) <- "SigPoisProcess_mult"
  return(results)
}


#---- Initialization functions
#' @export
SigPoisProcess_mult.PriorParams <- function(a = 1.01, alpha = 1.01,
                                            epsilon = 0.001,
                                            a0 = NULL,
                                            b0 = NULL,
                                            c0 = NULL,
                                            d0 = NULL,
                                            shrinkage = "none",
                                            SigPrior = NULL){
  list(a = a, alpha = alpha, epsilon = epsilon,
       a0 = a0, b0 = b0, c0 = c0, d0 = d0,
       shrinkage = shrinkage,
       SigPrior = SigPrior)
}


#' @export
SigPoisProcess_mult.controls <- function(nsamples = 2000,
                                        burnin = 5000,
                                        maxiter = 1e5,
                                        n_iter_betas = 2,
                                        tol = 1e-7,
                                        ncores = 1){
  list(nsamples = nsamples, burnin = burnin,
       maxiter = maxiter, tol = tol,
       ncores = ncores, n_iter_betas = n_iter_betas)
}

#' @export
SigPoisProcess_mult.init <- function(R_start = NULL,
                                     Theta_start = NULL,
                                     Betas_start = NULL,
                                     Mu_start = NULL){
  list(R_start = R_start,
       Theta_start = Theta_start,
       Betas_start = Betas_start,
       Mu_start = Mu_start)
}
