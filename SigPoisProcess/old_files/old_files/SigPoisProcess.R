#' Poisson process factorization for mutaitonal signatures analysis
#'
#' @param Mutations Dataframe containing the mutations
#' @param X Covariates
#' @param X_total Total surface of the Poisson process
#' @param method Which method to use for the estimation. Available options are 'gibbs', 'cavi' and 'map'
#' @param prior_params The hyperparameters of the prior.
#' @param controls Control parameters
#'
#' @returns An object
#' @export
#' @useDynLib SigPoisProcess
#' @import Rcpp
#' @import RcppArmadillo
SigPoisProcess <- function(Mutations,
                           X,
                           X_total,
                           method = "gibbs",
                           K = 10,
                           prior_params = list(),
                           controls = list(),
                           init = list()) {
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
      colnames(SigPrior_new) <- paste0("Sig_new", 1:K)
      rownames(SigPrior_new) <- rownames(SigPrior)
      SigPrior <- cbind(SigPrior, SigPrior_new)
      K <- Kpre + K
    } else {
      K <- Kpre
    }
  } else {
    SigPrior <- matrix(alpha, nrow = I, ncol = K)
    colnames(SigPrior) <- paste0("Sig", 1:K)
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
    # Calculate Signature-covariate assignment for each mutation
    Sig_Cov_prob <- calculate_Sig_covariate_prob(X, Sigs, Loads, channel_id - 1, sample_id - 1)
    #---- Store results in a list
    results$mle <- out_mle
    results$Signatures <- Sigs
    results$Betas <- Loads
    results$Mu <- NULL
    results$SigPrior <- SigPrior
    results$Sig_Cov_prob <- Sig_Cov_prob
    results$Xbar <- X_bar
  }
  ######################################################
  # MAP
  ######################################################
  if (method == "map") {
    cat("Calculating the maximum-a-posteriori for the model \n")
    maxiter <- controls$maxiter
    tol <- controls$tol
    use_acceleration <- controls$nesterov
    merge_move <- controls$merge_move
    # Run the search for the map
    if(use_acceleration == TRUE){
      out_map <- compute_SigPoisProcess_MAP_Nest(R_start, Betas_start, Mu_start, SigPrior,
        a, a0, b0, X, X_bar, channel_id - 1, sample_id - 1, maxiter, tol, Betas_scaled)
    } else {
      out_map <- compute_SigPoisProcess_MAP(R_start = R_start, Betas_start = Betas_start,
                                            Mu_start = Mu_start, SigPrior = SigPrior,
                                            a = a, a0 = a0, b0 = b0, X = X, X_bar = X_bar,
                                            channel_id = channel_id - 1, sample_id = sample_id - 1,
                                            maxiter = maxiter, tol = tol, scaled = Betas_scaled, merge_move = merge_move)
    }

    #--- Postprocess results
    Sigs <- out_map$R
    colnames(Sigs) <- colnames(SigPrior)
    rownames(Sigs) <- rownames(SigPrior)
    Loads <- out_map$Betas
    dimnames(Loads) <- list(colnames(SigPrior), levels(Mutations$sample), colnames(X))
    mu <- out_map$Mu
    rownames(mu) <- colnames(SigPrior)
    colnames(mu) <- colnames(X)
    # Calculate Signature-covariate assignment for each mutation
    Sig_Cov_prob <- calculate_Sig_covariate_prob(X, Sigs, Loads, channel_id - 1, sample_id - 1)
    #---- Store results in a list
    results$map <- out_map
    results$Signatures <- Sigs
    results$Betas <- Loads
    results$Mu <- mu
    results$SigPrior <- SigPrior
    results$Sig_Cov_prob <- Sig_Cov_prob
    results$Xbar <- X_bar
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
    A <- rep(NA, nsamples)
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
      # Sample a
      if(controls$sample_a){
        Mu_all <- array(NA, dim = c(K, J, L))
        for(j in 1:J) Mu_all[, j, ] <- Mu
        a <- sample_Post_a(1, epsilon, Betas, Mu_all, X_bar)
        #print(a)
        #if(a < 0.01){
        #  browser()
        #  }
      }
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
    Loads <- apply(BETAS, c(2, 3, 4), mean)
    rownames(mu) <- colnames(SigPrior)
    #---------- Save it
    results$Signatures <- Sigs
    results$Betas <- Loads
    results$Mu <- mu
    results$gibbs$MU <- MU
    results$gibbs$A <- A
  }

  #------ CAVI
  if (method == "cavi") {
    cat("Mean-field approximation of the posterior via CAVI \n")
    maxiter <- controls$maxiter
    tol <- controls$tol
    #------ Initialize quantities
    # Signatures
    alpha_r <- matrix(rgamma(I * K, SigPrior), nrow = I, ncol = K) + 1e-10 # Nugget for the gamma
    # Loadings parameters
    a_beta <- array(rgamma(K * J * L, a), dim = c(K, J, L))
    b_beta <- array(rgamma(K * J * L, a / epsilon), dim = c(K, J, L))
    # Relevance weights
    a_mu <- matrix(2, nrow = K, ncol = L) # matrix(rgamma(K * L, a0), nrow = K, ncol = L)
    b_mu <- matrix(2, nrow = K, ncol = L) # matrix(rgamma(K * L, b0), nrow = K, ncol = L)
    # Mutation signature-covariate probabilities
    Phi <- array(1 / (K * L), dim = c(N, K, L))
    # Run the CAVI algorithm for the inhomogeneous Poisson factorization
    logX <- log(X + 1e-18)
    out_cavi <- InhomogeneousPoissonNMF_CAVI(Phi, alpha_r, a_beta, b_beta,
      a_mu, b_mu, logX, X_bar,
      channel_id - 1,
      sample_id - 1,
      SigPrior,
      a = a, a0, b0, maxiter = maxiter
    )
    #---------- Postprocess output
    Sigs <- apply(out_cavi$alpha_r, 2, function(x) x / sum(x))
    colnames(Sigs) <- colnames(SigPrior)
    rownames(Sigs) <- rownames(SigPrior)
    Loads <- out_cavi$a_beta / out_cavi$b_beta
    mu <- out_cavi$b_mu / (out_cavi$a_mu - 1)
    rownames(mu) <- colnames(SigPrior)
    #---------- Save it
    results$cavi <- out_cavi
    results$Signatures <- Sigs
    results$Betas <- Loads
    results$Mu <- mu
  }

  return(results)
}


#---- Initialization functions
#' @export
SigPoisProcess.PriorParams <- function(a = 1.01, alpha = 1.01,
                                       epsilon = 0.001,
                                       Betas_scaled = TRUE,
                                       a0 = NULL, b0 = NULL,
                                       SigPrior = NULL){
  list(a = a, alpha = alpha, epsilon = epsilon,
       Betas_scaled = Betas_scaled,
       a0 = a0, b0 = b0, SigPrior = SigPrior)
}


#' @export
SigPoisProcess.controls <- function(nsamples = 2000,
                                    burnin = 5000,
                                    sample_a = FALSE,
                                    maxiter = 1e6,
                                    tol = 1e-6,
                                    nrepl = 1,
                                    ncores = 1,
                                    nesterov = FALSE,
                                    merge_move = TRUE){
  list(nsamples = nsamples, burnin = burnin, sample_a = sample_a,
       maxiter = maxiter, tol = tol, nrepl = nrepl,
       ncores = ncores,
       nesterov = nesterov,
       merge_move = merge_move)
}

#' @export
SigPoisProcess.init <- function(R_start = NULL, Betas_start = NULL, Mu_start = NULL){
  list(R_start = R_start, Betas_start = Betas_start, Mu_start = Mu_start)
}



