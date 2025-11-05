#' Poisson process factorization for mutational signatures analysis with genomic covariates
#'
#' @param gr_Mutations GenomicRanges object where \code{mcols} contains
#' @param SignalTrack A matrix whose columns containing the genomic covariates for the whole genome.
#'                    Each row indicates the value of the covariate in a bin.
#' @param CopyTrack A matrix whole columns contain the number of copies multiplied by
#'                  the size of the bin. For example, if Patient_01 has 3 copies in the
#'                  first bin and the length of the bin is 2000 base pairs, then set
#'                  \code{'CopyTrack[1,1] <- 3/2* 2000'}. By default, we suggest dividing the copy
#'                  number by 2, so that 1 indicates normal conditions.
#' @param K Upper bound to the number of signatures expected in the data
#' @param method Which method to use for the estimation. Available options are \code{'mle', 'map', 'mcmc'}
#'               Beware of time
#' @param prior_params The hyperparameters of the prior. \code{SigPoisProcess.PriorParams}
#' @param init Initialization parameters. See \code{SigPoisProcess.init}
#' @param controls Control parameters for each method. See \code{SigPoisProcess.contol}
#' @param verbose Weather to track the output and the logposterior. Setting \code{verbose = TRUE} prints every 10 iterations.
#' @param update_Signatures Update the signature matrix in the model. Default is \code{TRUE}
#' @param update_Theta Update the Baseline exposures in the model. Default is \code{TRUE}
#' @param update_Betas Update the regression coefficients in the model. Default is \code{TRUE}
#'
#' @returns An object of class SigPoisProcess
#' @export
#' @useDynLib SigPoisProcess, .registration = TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom magrittr %>%
SigPoisProcess <- function(gr_Mutations,
                           SignalTrack,
                           CopyTrack,
                           K = 10,
                           method = "map",
                           prior_params = SigPoisProcess.PriorParams(),
                           controls = SigPoisProcess.controls(),
                           init = SigPoisProcess.init(),
                           verbose = FALSE,
                           update_Signatures = TRUE,
                           update_Theta = TRUE,
                           update_Betas = TRUE){

  # Extract copy number
  #copy_number <- gr_Mutations$copy_number
  #gr_Mutations$copy_number <- NULL
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

  # Shuffle copy tracks
  CopyTrack <- CopyTrack[, colnames(MutMatrix)]

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
  shrinkage <- "mu_sigma"
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
    Area_patients <- colSums(CopyTrack)
    Theta_start <- matrix(rgamma(J * K, colSums(MutMatrix)), nrow = K, byrow = TRUE)#/sum(bin_weight)
    Theta_start <- Theta_start / t(Area_patients)[rep(1, K), ]
  } else {
    Theta_start <- init$Theta_start
  }
  #---- Initial value for the relevance weights
  if(is.null(init$Mu_start)){
    Mu_start <- rowMeans(Theta_start)#1 / rgamma(K, a * J  + 1, epsilon * a * J)
  } else {
    Mu_start <- init$Mu_start
  }

  #---- Initial value for the variances
  if(is.null(init$Sigma2_start)){
    Sigma2_start <- rep(d0/(c0 + 1), K)#1 / rgamma(K, a * J  + 1, epsilon * a * J)
  } else {
    Sigma2_start <- init$Sigma2_start
  }

  #---- Initial value for the regression coefficients (start at 0)
  if(is.null(init$Betas_start)){
    Betas_start <- matrix(rnorm(p * K, mean = 0, sd = 1e-7), nrow = p, ncol = K)#matrix(0, nrow = p, ncol = K)
  } else {
    Betas_start <- init$Betas_start
  }

  #----------------------------------------------- Control parameters for the model
  if (!method %in% c("mle", "map", 'mcmc')) {
    stop("method must be either 'map', 'mle' or 'mcmc'")
  }
  controls <- do.call("SigPoisProcess_mult.controls", controls)

  #----------------------------------------------- Estimate the model
  # List to append results
  results <- list(sol = NULL, method = method, chain = NULL)

  ######################################################
  # Optimization
  ######################################################
  if(method == "map") {
    if(verbose){
      cat("Start Optimization \n")
    }
    maxiter <- controls$maxiter
    n_iter_betas <- controls$n_iter_betas
    tol <- controls$tol
    # Find the MLE using the multiplicative rules
    out <- .PoissonProcess_optim_CN(R_start = R_start,
                                    Theta_start = Theta_start,
                                    Betas_start = Betas_start,
                                    Mu_start = Mu_start,
                                    Sigma2_start = Sigma2_start,
                                    X = X,
                                    SignalTrack = SignalTrack,
                                    CopyTrack = CopyTrack,
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
                                    tol = tol,
                                    verbose = verbose)

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
  }


  ######################################################
  # MCMC
  ######################################################
  if(method == "mcmc") {
    if(verbose){
      cat("Start MCMC solution \n")
    }
    nsamples <- controls$nsamples
    burnin <- controls$burnin
    out <- .PoissonProcess_BayesCN(R_start = R_start,
                                   Theta_start = Theta_start,
                                   Betas_start = Betas_start,
                                   Mu_start = Mu_start,
                                   Sigma2_start = Sigma2_start,
                                   X = X,
                                   SignalTrack = SignalTrack,
                                   CopyTrack = CopyTrack,
                                   nsamples = nsamples,
                                   burn_in = burnin,
                                   channel_id = channel_id - 1,
                                   sample_id = sample_id - 1,
                                   SigPrior = SigPrior,
                                   update_R = update_Signatures,
                                   update_Theta = update_Theta,
                                   update_Betas = update_Betas,
                                   a = a,
                                   a0 = a0,
                                   b0 = b0,
                                   c0 = c0,
                                   d0 = d0,
                                   verbose = verbose)

    # Add dimension names
    dimnames(out$SIGSchain) <- list(NULL,  rownames(SigPrior), colnames(SigPrior))
    dimnames(out$THETAchain) <- list(NULL,  colnames(SigPrior), colnames(MutMatrix))
    dimnames(out$BETASchain) <- list(NULL,  colnames(X), colnames(SigPrior))
    dimnames(out$MUchain) <- list(NULL, colnames(SigPrior))
    dimnames(out$SIGMA2chain) <- list(NULL, colnames(SigPrior))

    #--- Postprocess results
    # Signatures
    Sigs <- apply(out$SIGSchain[-c(1:burnin), , ], c(2,3), mean)
    # Loadings
    Thetas <- apply(out$THETAchain[-c(1:burnin), , ], c(2,3), mean)
    # Coefficients
    Betas <- apply(out$BETASchain[-c(1:burnin), , ], c(2,3), mean)
    # Mu and sigma
    Mu <- colMeans(out$MUchain[-c(1:burnin), ])
    Sigma2 <- colMeans(out$SIGMA2chain[-c(1:burnin), ])

    #---- Store results in a list
    results$chain <- out
  }

  #---- Store results in a list
  results$Thetas <- Thetas
  results$Signatures <- Sigs
  results$Betas <- Betas

  if(method == "map"){
    results$Mu <- c(out$Mu)
    names(results$Mu) <- colnames(SigPrior)
    results$Sigma2 <- c(out$Sigma2)
    names(results$Sigma2) <- colnames(SigPrior)
  } else if (method == "mcmc") {
    results$Mu <- Mu
    results$Sigma2 <- Sigma2
  } else {
    results$Mu <- NULL
    results$Sigma2 <- NULL
  }

  results$init <- list(R_start = R_start,
                       Theta_start = Theta_start,
                       Betas_start = Betas_start,
                       Mu_start = Mu_start,
                       Sigma2_start = Sigma2_start)

  results$PriorParams <- list(a = a, alpha = alpha, epsilon = epsilon,
                              a0 = a0, b0 = b0, c0 = c0, d0 = d0,
                              SigPrior = SigPrior)

  results$controls <- controls

  class(results) <- "SigPoisProcess"
  return(results)
}


#' Prior parameters in the model
#'
#' @param a Prior parameter in the baseline coefficients (gamma prior).
#' @param alpha Prior parameter for the signature (Dirichlet prior).
#' @param epsilon Defaults to 0.001.
#' @param a0 Prior for mu ~ InvGamma(a0, b0). Defaults to NULL (compressive).
#' @param b0 Prior for mu ~ InvGamma(a0, b0). Defaults to NULL (compressive).
#' @param c0 Prior over sigma2 ~ InvGamma(c0, d0).
#' @param d0 sigma2 ~ InvGamma(c0, d0).
#' @param SigPrior Informative prior matrix for the Dirichlet.
#'
#' @export
SigPoisProcess.PriorParams <- function(a = 1.01, alpha = 1.01,
                                            epsilon = 0.001,
                                            a0 = NULL,
                                            b0 = NULL,
                                            c0 = 100,
                                            d0 = 1,
                                            SigPrior = NULL){
  list(a = a, alpha = alpha, epsilon = epsilon,
       a0 = a0, b0 = b0, c0 = c0, d0 = d0,
       SigPrior = SigPrior)
}


#' Control values for the estimation algorithm
#'
#' @param nsamples Number of mcmc iterations
#' @param burnin Number of burnin iterations
#' @param maxiter Maximum number of iteratiors for the optimization
#' @param n_iter_betas Number of fisher scoring per iteration in the optimization
#' @param tol Tolerance for convergence
#'
#' @export
SigPoisProcess.controls <- function(nsamples = 2000,
                                    burnin = 2000,
                                    maxiter = 5000,
                                    n_iter_betas = 2,
                                    tol = 1e-6){
  list(nsamples = nsamples, burnin = burnin,
       maxiter = maxiter, tol = tol,
       n_iter_betas = n_iter_betas)
}

#' Initial values for the estimation algorithm
#'
#' @param R_start Start value for the signature matrix.
#' @param Theta_start Start value for the baseline  matrix.
#' @param Betas_start Start value for the regression coefficient matrix.
#' @param Mu_start Start value for the relevance weights.
#' @param Sigma2_start Start value for the variance of the regression coefficients.
#'
#' @export
SigPoisProcess.init <- function(R_start = NULL,
                               Theta_start = NULL,
                               Betas_start = NULL,
                               Mu_start = NULL,
                               Sigma2_start = NULL){
  list(R_start = R_start,
       Theta_start = Theta_start,
       Betas_start = Betas_start,
       Mu_start = Mu_start,
       Sigma2_start = Sigma2_start)
}

#' Sample a mutational signature matrix
#'
#' @param Alpha Matrix of dimension I x K.
#' @export
sample_signatures <- function(Alpha) {
  # Alpha is a matrix of I x K
  I <- nrow(Alpha)
  Ktot <- ncol(Alpha)
  R <- matrix(rgamma(I * Ktot, Alpha), nrow = I) + 1e-10 # Small nugget to avoid degeneracies of the gamma prior
  R <- apply(R, 2, function(x) x/sum(x))
  colnames(R) <- colnames(Alpha)
  return(R)
}
