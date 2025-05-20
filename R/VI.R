InhomogeneousPoissonNMF_VI <- function(df_trinucl, X, X_total,
                                       K = 10,
                                       alpha = 0.5,
                                       SigPrior = NULL,
                                       #genome = BSgenome.Hsapiens.UCSC.hg19,
                                       a = 1,
                                       epsilon = 0.001,
                                       tol = 1e-7,
                                       maxiter = 100){

  # Find trinucleotide channels for each mutation in maf file
  #df_trinucl <- get_df_trinucl(maf, genome = genome)

  # Useful quantities
  channel_id <- as.numeric(df_trinucl$Channel)
  sample_id <- as.numeric(df_trinucl$sample)

  I <- length(levels(df_trinucl$Channel)) # Assume the usual 96 single base substitution signatures
  J <- max(sample_id)
  L <- ncol(X)
  N <- nrow(df_trinucl)

  # Compressive prior for Mu
  a0 <- a * J + 1
  b0 <- epsilon * a * J

  # Prior for the signatures
  if(!is.null(SigPrior)){
    if(length(setdiff(rownames(SigPrior), levels(df_trinucl$Channel))) > 0){
      stop("Must have rownames(SigPrior) equal to levels(df_trinucl$Channel)")
    }
    Kpre <- ncol(SigPrior)
    SigPrior <- SigPrior[levels(df_trinucl$Channel), ] # reorder
    if(K > 0){
      SigPrior_new <- matrix(alpha, nrow = I, ncol = K)
      colnames(SigPrior_new) <- paste0("SBSnew", 1:K)
      rownames(SigPrior_new) <- levels(df_trinucl$Channel)
      SigPrior <- cbind(SigPrior, SigPrior_new)
      K <- Kpre + K
    }
  } else {
    SigPrior <- matrix(alpha, nrow = I, ncol = K)
    colnames(SigPrior) <- paste0("SBSnew", 1:K)
    rownames(SigPrior) <- levels(df_trinucl$Channel)
  }

  #---------------------Initialize the array of total surface of poisson process
  X_bar <- array(NA, dim = c(K, J, L))
  for(k in 1:K) X_bar[k, , ] <- X_total

  #------------------ Initialize variational parameters at random starting point
  # Signatures
  alpha_r <- matrix(rgamma(I * K, SigPrior), nrow = I , ncol = K) + 1e-10 # Nugget for the gamma
  # Loadings parameters
  a_beta <- array(rgamma(K * J * L, a), dim = c(K, J, L))
  b_beta <- array(rgamma(K * J * L, a / epsilon), dim = c(K, J, L))
  # Relevance weights
  a_mu <- matrix(2, nrow = K, ncol = L)#matrix(rgamma(K * L, a0), nrow = K, ncol = L)
  b_mu <- matrix(2, nrow = K, ncol = L)#matrix(rgamma(K * L, b0), nrow = K, ncol = L)
  # Mutation signature-covariate probabilities
  Phi <- array(1/(K * L), dim = c(N, K, L))
  # Run the CAVI algorithm for the inhomogeneous Poisson factorization
  logX <- log(X + 1e-18)
  out_cavi <- InhomogeneousPoissonNMF_CAVI(Phi, alpha_r, a_beta, b_beta,
                                           a_mu, b_mu, logX, X_bar, channel_id, sample_id, SigPrior,
                                           a = a, a0, b0, maxiter = maxiter)

  # post-process the signatures and the weights
  Sigs <- apply(out_cavi$alpha_r, 2, function(x) x/sum(x))
  colnames(Sigs) <- colnames(SigPrior)
  rownames(Sigs) <- rownames(SigPrior)
  mu <- out_cavi$b_mu/(out_cavi$a_mu - 1)
  rownames(mu) <- colnames(SigPrior)
  return(list("Signatures" = Sigs,
              "Loadings" = out_cavi$a_beta/out_cavi$b_beta,
              "RelWeights" = mu,
              out_cavi = out_cavi))
}
