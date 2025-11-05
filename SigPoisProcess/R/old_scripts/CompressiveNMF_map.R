
#' @export
#' @useDynLib SigPoisProcess
CompressiveNMF_map <- function(X, K, S = NULL, alpha = 1, a = 1,
                               epsilon = 0.001,
                               a0 = a * ncol(X) + 1,
                               b0 = a * epsilon * ncol(X),
                               cutoff = 5 * epsilon,
                               tol = 1e-6,
                               maxiter = 1e6){



  # Check if the matrix X is made of integers
  if(any(round(X) != X)){
    stop("The matrix X has to be made of integers! Consider using round(X).")
  }
  # Matrix dimensions
  I <- nrow(X); J <- ncol(X)

  # Check if the matrix S is correctly specified
  if(!is.null(S)){
    if(!is.matrix(S)){
      stop("The hyperparameter S must be a matrix!")
    }
  }

  # Fix the prior matrix for the signatures
  SignaturePrior <- cbind(S, matrix(alpha, nrow = I, ncol = K))
  colnames(SignaturePrior) <- c(colnames(S), paste0(paste0("Sig_new", c(1:K))))
  rownames(SignaturePrior) <- rownames(X)
  SignaturePrior <- pmax(SignaturePrior, 1e-10) # correction for Dirichlet

  # Random initialization from the prior
  Ktot <- ncol(SignaturePrior)
  R <- CompressiveNMF:::sample_signatures(SignaturePrior)
  Mu <- 1/rgamma(Ktot, a0, b0)
  shape_mat <- matrix(a, nrow = Ktot, ncol = J)
  rate_mat <- as.matrix(a/Mu)[, rep(1, J)]
  Theta <- CompressiveNMF:::sample_weights(shape_mat, rate_mat)
  # Run the code in Rcpp
  res <- compute_CompressiveNMF_MAP(X, R, Theta, Mu, SignaturePrior,
                                    a0, b0, a,
                                    maxiter, tol)
  # Postprocess results
  ids <- which(c(res$Mu) > cutoff)
  R <- res$R[, ids]
  colnames(R) <- colnames(SignaturePrior[, ids]); rownames(R) <- rownames(X)
  Theta <- res$Theta[ids, ]
  rownames(Theta) <- colnames(SignaturePrior[, ids]); colnames(Theta) <- colnames(X)
  Mu <- c(res$Mu)[ids]
  names(Mu) <- colnames(SignaturePrior[, ids])

  return(list(Signatures = R,
              Theta = Theta,
              Mu = Mu,
              mapOutput = res,
              a0 = a0,
              b0 = b0,
              a = a,
              SignaturePrior = SignaturePrior))
}





