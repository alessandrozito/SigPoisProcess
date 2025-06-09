#' Poisson process factorization for topographic-dependent mutational signatures analysis
#'
#' @param gr_Mutations GenomicRanges object contatining covariates and samples
#' @param df_areas data.frame containing the total areas of each covariate for each sample
#' @param add_baseline Wheter to add a baseline equal to 1 to the whole model.
#' @param method Which method to use for the estimation. Available options are 'mle', 'map' and 'gibbs'
#' @param prior_params The hyperparameters of the prior.
#' @param controls Control parameters for each method. See \code{SigPoisProcess.contol}
#' @returns An object of class SigPoisProcess.
#' @export
#' @useDynLib SigPoisProcess
#' @import Rcpp
#' @import RcppArmadillo
SigPoisProcess <- function(gr_Mutations,
                           df_areas,
                           method = "map",
                           add_baseline = TRUE,
                           K = 10,
                           prior_params = SigPoisProcess.PriorParams(),
                           controls = SigPoisProcess.controls(),
                           init = SigPoisProcess.init()) {

  # Extract the data and the covariates
  X <- as.data.frame(mcols(gr_Mutations)) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()

  X_totals <- df_areas %>%
    column_to_rownames(var = "sample") %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()

  # Extract mutation data
  Mutations <- as.data.frame(mcols(gr_Mutations)) %>%
    dplur::select(sample, channel)
  Mutations$sample <- factor(Mutations$sample, levels = rownames(X_totals))
  Mutations$channel <- as.factor(Mutations$channel)

  # Add baseline
  if (add_baseline) {
    # Pad baseline to the data
    X <- cbind("baseline" = 1, X)
    # Load the Genome and get the chromosome length
    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    std_chrs <- paste0("chr", c(1:22, "X", "Y"))
    chrom_lengths <- seqlengths(genome)[1:24]
    X_totals <- cbind("baseline" = sum(1 * chrom_lengths), X_totals)
  }

  # Extract channel and sample Id from the mutation list
  if (!is.factor(Mutations$channel)) {
    stop("Variable channel in Mutations has to be a factor")
  }

  if (!is.factor(Mutations$sample)) {
    stop("Variable sample in Mutations has to be a factor")
  }

  # Extract channel ID and mutation ID
  channel_id <- as.numeric(Mutations$channel)
  sample_id <- as.numeric(Mutations$sample)

  #----------------------------------------------- Model dimensions
  I <- length(unique(channel_id)) # Number of mutational channels
  J <- max(sample_id) # Number of patients
  L <- ncol(X) # Number of covariates
  N <- nrow(Mutations) # Total number of mutations

  ## ----------------------------------------------- Model settings
  prior_params <- do.call("SigPoisProcess.PriorParams", prior_params)
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
  Betas_scaled <- prior_params$Betas_scaled
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

  #----------------------------------------------- Total area
  X_bar <- array(NA, dim = c(K, J, L))
  for (k in 1:K) X_bar[k, , ] <- X_total

  #----------------------------------------------- Initial values
  # Initialize the model
  init <- do.call("SigPoisProcess.init", init)
  #---- Initial value for the signatures
  if(is.null(init$R_start)){
    R_start <- sample_signatures(SigPrior)
  } else {
    R_start <- init$R_start
    if(ncol(R_start) != K){
      stop("ncol(R_start) must be equal to specified K")
    }
  }
  #---- Initial value for the relevance weights
  if(is.null(init$Mu_start)){
    Mu_start <- matrix(1 / rgamma(K * L, 2 * a * J + 1, epsilon * a * J), nrow = K, ncol = L)
  } else {
    Mu_start <- init$Mu_start
  }
  #---- Initial value for the loadings
  if(is.null(init$Betas_start)){
    Betas_start <- array(rgamma(K * J * L, a, a / Mu_start), dim = c(K, J, L))
  } else {
    Betas_start <- init$Betas_start
  }

  #----------------------------------------------- Control parameters for the model
  if (!method %in% c("cavi", "gibbs", "map", "mle")) {
    stop("method must be either 'gibbs', 'map', 'cavi' or 'mle'")
  }
  controls <- do.call("SigPoisProcess.controls", controls)

  #----------------------------------------------- Estimate the model
  # List to append results
  results <- list(cavi = NULL, gibbs = NULL, map = NULL, mle = NULL, method = method)
  ######################################################
  # MLE
  ######################################################
  if (method == "mle") {
    cat("Calculating the MLE for the model \n")
    maxiter <- controls$maxiter
    tol <- controls$tol
    # Find the MLE using the multiplicative rules
    out_mle <- compute_SigPoisProcess_MLE(R_start = R_start, Betas_start = Betas_start,
                                          X = X, X_bar = X_bar,
                                          channel_id = channel_id - 1, sample_id = sample_id - 1,
                                          maxiter = maxiter, tol = tol)
    #--- Postprocess results
    Sigs <- out_mle$R
    colnames(Sigs) <- colnames(SigPrior)
    rownames(Sigs) <- rownames(SigPrior)
    Loads <- out_mle$Betas
    dimnames(Loads) <- list(colnames(SigPrior), levels(Mutations$sample), colnames(X))

    #---- Store results in a list
    results$mle <- out_mle
    results$Signatures <- Sigs
    results$Betas <- Loads
    results$Mu <- NULL
    results$SigPrior <- SigPrior

  }
  ######################################################
  # MAP
  ######################################################
  if (method == "map") {
    cat("Calculating the maximum-a-posteriori for the model \n")
    maxiter <- controls$maxiter
    tol <- controls$tol
    merge_move <- controls$merge_move
    out_map <- compute_SigPoisProcess_MAP(R_start = R_start, Betas_start = Betas_start,
                                            Mu_start = Mu_start, SigPrior = SigPrior,
                                            a = a, a0 = a0, b0 = b0, X = X, X_bar = X_bar,
                                            channel_id = channel_id - 1, sample_id = sample_id - 1,
                                            maxiter = maxiter, tol = tol, scaled = Betas_scaled, merge_move = merge_move)

    #--- Postprocess results
    Sigs <- out_map$R
    colnames(Sigs) <- colnames(SigPrior)
    rownames(Sigs) <- rownames(SigPrior)
    Loads <- out_map$Betas
    dimnames(Loads) <- list(colnames(SigPrior), levels(Mutations$sample), colnames(X))
    mu <- out_map$Mu
    rownames(mu) <- colnames(SigPrior)
    colnames(mu) <- colnames(X)

    #---- Store results in a list
    results$map <- out_map
    results$Signatures <- Sigs
    results$Betas <- Loads
    results$Mu <- mu
    results$SigPrior <- SigPrior

  }

  ######################################################
  # Gibbs
  ######################################################
  if (method == "gibbs") {
    cat("Runing Gibbs sampler for the full posterior \n")
    nsamples <- controls$nsamples
    burnin <- controls$burnin

    # Initialize the sampler
    Mu <- Mu_start
    Betas <- Betas_start
    R <- R_start
    # Store the output
    MU <- array(NA, c(nsamples, K, L))
    BETAS <- array(NA, c(nsamples, K, J, L))
    SIGNS <- array(NA, dim = c(nsamples, I, K))
    # Run the Sampling steps
    pb <- txtProgressBar(style = 3)
    for (iter in 1:(nsamples + burnin)) {
      setTxtProgressBar(pb, iter / (nsamples + burnin))
      # Allocate counts
      Counts <- sample_Counts_cpp(
        X = X, R = R,
        Betas = Betas,
        channel_Id = channel_id - 1,
        sample_Id = sample_id - 1
      )
      # Sample the signatures
      Alpha <- Counts$Y_signatures + SigPrior
      R <- sample_signatures(Alpha)
      # Sample the Beta loadings
      #Betas <- sample_Betas(Counts, Mu, X_bar, a)
      Betas <- sample_Betas_scaled(Counts, Mu, X_bar, a)
      # Sample Mu
      #Mu <- sample_Mu(Beta = Betas, epsilon = epsilon, a = a)
      Mu <- sample_Mu_scaled(Betas, X_bar, a, a0, b0)
      if (iter > burnin) {
        SIGNS[iter - burnin, , ] <- R
        BETAS[iter - burnin, , , ] <- Betas
        MU[iter - burnin, , ] <- Mu
        A[iter - burnin] <- a
      }
    }
    close(pb)

    #------- Postprocess the output
    Sigs <- apply(SIGNS, c(2, 3), mean)
    colnames(Sigs) <- colnames(SigPrior)
    rownames(Sigs) <- rownames(SigPrior)
    mu <- apply(MU, c(2, 3), mean)
    rownames(mu) <- colnames(SigPrior)
    colnames(mu) <- colnames(X)
    Loads <- apply(BETAS, c(2, 3, 4), mean)
    dimnames(Loads) <- list(colnames(SigPrior), levels(Mutations$sample), colnames(X))

    #---------- Save it
    results$Signatures <- Sigs
    results$Betas <- Loads
    results$Mu <- mu
    results$gibbs$MU <- MU
    results$SigPrior <- SigPrior
    results$Xbar <- X_bar
  }

  # Calculate Signature-covariate assignment for each mutation
  Sig_Cov_prob <- calculate_Sig_covariate_prob(X, Sigs, Loads, channel_id - 1, sample_id - 1)
  Probs <- Sig_Cov_prob$Probs
  dimnames(Probs)[[2]] <- colnames(Sigs)
  dimnames(Probs)[[3]] <- dimnames(Loads)[[3]]
  MaxProbIds <- as.data.frame(Sig_Cov_prob$MaxProbIds) %>%
    group_by(V1, V2) %>%
    summarise(Fraction = n()/nrow(Mutations)) %>%
    as.data.frame()
  MaxProbIds$V1 <- colnames(Sigs)[MaxProbIds$V1]
  MaxProbIds$V2 <- dimnames(Loads)[[3]][MaxProbIds$V2]
  colnames(MaxProbIds) <- c("Signature", "Covariate", "Fraction")
  results$Probs <- Probs
  results$AttrSummary <- MaxProbIds
  results$Xbar <- X_bar

  class(results) <- "SigPoisProcess"
  return(results)
}


#---- Initialization functions
#' @export
SigPoisProcess.PriorParams <- function(a = 1.01, alpha = 1.01,
                                       epsilon = 0.001,
                                       Betas_scaled = TRUE,
                                       a0 = NULL,
                                       b0 = NULL,
                                       SigPrior = NULL){
  list(a = a, alpha = alpha, epsilon = epsilon,
       Betas_scaled = Betas_scaled,
       a0 = a0, b0 = b0, SigPrior = SigPrior)
}


#' @export
SigPoisProcess.controls <- function(nsamples = 2000,
                                    burnin = 5000,
                                    maxiter = 1e6,
                                    tol = 1e-6,
                                    nrepl = 1,
                                    ncores = 1,
                                    merge_move = TRUE){
  list(nsamples = nsamples, burnin = burnin,
       maxiter = maxiter, tol = tol, nrepl = nrepl,
       ncores = ncores,
       merge_move = merge_move)
}

#' @export
SigPoisProcess.init <- function(R_start = NULL,
                                Betas_start = NULL,
                                Mu_start = NULL){
  list(R_start = R_start, Betas_start = Betas_start,
       Mu_start = Mu_start)
}



