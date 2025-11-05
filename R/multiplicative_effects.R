library(rvest)
# Initialization of the SigPoisProcess method using standard NMF +
# Poisson regression with identity link.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
#library(SigPoisProcess)
library(rtracklayer)
library(foreach)
library(doParallel)
# Load simulated data
devtools::document()
data <- readRDS("~/SigPoisProcess/data/temp_simulation.rds.gzip")

set.seed(10)
J <- length(unique(data$gr_Mutations$sample))

# Covariates
Xsignal <- cbind(as.matrix(mcols(data$gr_SignalTrack)))
Xtotal <- array(NA,  dim = c(nrow(Xsignal), J, ncol(Xsignal)),
                dimnames = list(1:nrow(Xsignal),
                                unique(data$gr_Mutations$sample),
                                colnames(Xsignal)))
for(j in 1:J) Xtotal[, j, ] <- Xsignal
Bins <- t(matrix(width(data$gr_SignalTrack))[, rep(1, J)])


p <- ncol(Xsignal)

X <- as.matrix(mcols(data$gr_Mutations))[, -c(1,2)]
channel_id = as.numeric(data$gr_Mutations$channel) - 1
sample_id = as.numeric(data$gr_Mutations$sample) - 1

Xfield <- SigPoisProcess:::build_X_field(X = as.matrix(mcols(data$gr_Mutations))[, -c(1,2)],
                                         I = 96,
                                         J = J,
                                         channel_id = as.numeric(data$gr_Mutations$channel) - 1,
                                         sample_id = as.numeric(data$gr_Mutations$sample) - 1)

# Initialization.
K <- 8
MutMatrix <- getTotalMutations(data$gr_Mutations)
set.seed(10)
StartSol <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = K)

R_start <- StartSol$Signatures
Theta_start <- StartSol$Theta/sum(width(data$gr_SignalTrack))
Betas_start <- matrix(0, nrow = K, ncol = p)#matrix(rnorm(K * p, sd = 0.01), ncol = p)
colnames(Betas_start) <- colnames(Xsignal)

Betas <- Betas_start
R <- R_start
Theta <- Theta_start
# Update for theta
k <- 1
j <- 1
i <- 1

R_new <- Update_Signatures_MLE_loglinear(R, Theta, Betas, X, channel_id, sample_id)

Theta_new <- Update_Theta_MLE_loglinear(R, Theta, Betas, X,
                                        channel_id,
                                        sample_id,
                                        Xtotal, Bins)

Betas <- Update_Theta_Betas_loglinear(R, Theta, Betas, X, channel_id, sample_id, Xtotal, Bins)
plot(solve(Hess), Hes_inv)

plot(Hess[lower.tri(Hess)], t(Hess)[lower.tri(Hess)])

t(Hess) - Hess

Xfield
xx <- compute_expXtB_field(X, I = 96, J, channel_id, sample_id, Betas)
xx[1,1][[1]] * matrix(R[1, ])[, rep(1, 14)] == xx[1,1][[1]] * R[1, ]

pp <- xx[1,1][[1]] * R[1, ] * Theta[, 1]
pp <- apply(pp, 2, function(x) x/sum(x))
min(Xfield[1,1][[1]] * t(pp[1, ])[rep(1, 10), ] - tmp)



Theta_upd <- Theta
tt <- 0
for(i in 1:96){
  tots <- sapply(1:K, function(k) R[i,k] * exp(c(Betas[k, ] %*% Xfield[i, j][[1]])))
  if(!is.matrix(tots)){
    tots <- t(tots)
  }
  tt <- tt + sum(tots[, k] / rowSums(tots * t(Theta[, rep(j, nrow(tots))])))
}

bins <- width(data$gr_SignalTrack)

Theta[k, j] * tt/sum(bins * exp(Xtotal[, j, ] %*% Betas[k, ])) # <---- important update for theta :)

# Update for regression coefficients
Hess <- 0
grad_1 <- 0
for(j in 1:J){
  w <- c(Theta[k, j] * bins * exp(Xtotal[, j, ] %*% Betas[k, ]))
  Hess <- pp + crossprod(Xtotal[, j, ] * w, Xtotal[, j, ])
  grad_1 <- grad_1 + colSums(Xtotal[, j, ] * w)
}
solve(Hess)

grad_2 <- 0
for(j in 1:J){
  for(i in 1:96){
    if(MutMatrix[i, j] > 0){
      tots <- sapply(1:K, function(k) R[i,k] * exp(c(Betas[k, ] %*% Xfield[i, j][[1]])) * Theta[k, j])
      if(!is.matrix(tots)){
        tots <- t(tots)
      }
      probs <- apply(tots, 1, function(x) x/sum(x))[k, ]
      grad_2 <- grad_2 + rowSums(Xfield[i, j][[1]] * t(probs)[rep(1, p), ])
    }
  }
}

beta_k <- Betas[k, ] + solve(Hess) %*% (grad_2 - grad_1)


update_Theta <- function(Theta, R, Betas, Xfield, Xtotal, bins, MutMatrix){
  Theta_upd <- Theta
  for(k in 1:K){
    for(j in 1:J){
      tt <- 0
      for(i in 1:96){
        if(MutMatrix[i, j] > 0){
        tots <- sapply(1:K, function(k) R[i,k] * exp(c(Betas[k, ] %*% Xfield[i, j][[1]])))
        if(!is.matrix(tots)){
          tots <- t(tots)
        }
        tt <- tt + sum(tots[, k] / rowSums(tots * t(Theta[, rep(j, nrow(tots))])))
        }
      }
      Theta_upd[k, j] <- Theta[k, j] * tt/sum(bins * exp(Xtotal[, j, ] %*% Betas[k, ])) # <---- important update for theta :)
    }
  }
  return(Theta_upd)
}

update_Betas <- function(Theta, R, Betas, Xfield, Xtotal, bins, MutMatrix){
  Betas_upd <- Betas
  K <- nrow(Theta)
  J <- ncol(Theta)
  for(k in 1:K){
    Hess <- 0
    grad_1 <- 0
    for(j in 1:J){
      w <- c(Theta[k, j] * bins * exp(Xtotal[, j, ] %*% Betas[k, ]))
      Hess <- Hess + crossprod(Xtotal[, j, ] * w, Xtotal[, j, ])
      grad_1 <- grad_1 + colSums(Xtotal[, j, ] * w)
    }
    #solve(Hess)

    grad_2 <- 0
    for(j in 1:J){
      for(i in 1:96){
        if(MutMatrix[i, j] > 0){
          tots <- sapply(1:K, function(k) R[i,k] * exp(c(Betas[k, ] %*% Xfield[i, j][[1]])) * Theta[k, j])
          if(!is.matrix(tots)){
            tots <- t(tots)
          }
          probs <- apply(tots, 1, function(x) x/sum(x))[k, ]
          grad_2 <- grad_2 + rowSums(Xfield[i, j][[1]] * t(probs)[rep(1, p), ])
        }
      }
    }
    Betas_upd[k, ] <- Betas[k, ] + solve(Hess + 10 * diag(nrow(Hess))) %*% (grad_2 - grad_1)
  }
  return(Betas_upd)
}

Betas_upd <- Betas
Theta_upd <- Theta
for(it in 1:20){
  print(it)
  Theta_upd <- update_Theta(Theta_upd, R, Betas_upd, Xfield, Xtotal, bins, MutMatrix)
  Betas_upd_n <- update_Betas(Theta_upd, R, Betas_upd, Xfield, Xtotal, bins, MutMatrix)
  print(max(abs(Betas_upd_n/Betas_upd - 1)))
  Betas_upd <- Betas_upd_n
}

hist(Betas_upd - Betas)


hist(Theta_upd - Theta)
hist(Betas)
Betas_upd - Betas

round(Betas_upd, 8)
round(Theta_upd, 8)




# Install if needed: install.packages("einsum")
sol <- compute_MAP_loglinear(R_start = R_start,
                             Theta_start = Theta_start,
                             Betas_start = Betas_start,
                             X = X,
                             channel_id = channel_id,
                             sample_id = sample_id,
                             Xtotal = Xtotal,
                             Bins = Bins,
                             maxiter = 50)

sol$Betas
sol$Theta

T1 <- update_Theta(Theta, R, Betas_upd, Xfield, Xtotal, bins, MutMatrix)
T2 <- update_Theta_vectorized(Theta, R, Betas_upd, Xfield, Xtotal, bins, MutMatrix)



