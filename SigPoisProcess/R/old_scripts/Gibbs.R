InhomogeneousPoissonNMF <- function(df_trinucl,
                                    X, X_total,
                                    S = NULL,
                                    betah = 1000,
                                    R_start = NULL,
                                    K = 5,
                                    alpha = 0.5,
                                    a = 1,
                                    epsilon = 0.001,
                                    nsamples = 20,
                                    burnin = 20){

  # Initialize the sampler
  J <- max(df_trinucl$sample_id)
  p <- ncol(X)
  I <- 96
  if(!is.null(R_start)){
    K <- ncol(R_start)
  }

  if(!is.null(S)){
    H <- ncol(S)
    Ktot <- K + H
    prior_names <- colnames(S)
    if(K > 0){
      prior_names <- c(prior_names, paste0("SBSnew", c(1:K)))
    }
  } else {
    Ktot <- K
    prior_names <- paste0("SBSnew", c(1:K))
  }

  K <- Ktot


  # Prepare the hyperparameters in the signatures
  SignaturePrior <- matrix(alpha, nrow = I, ncol = K)
  if(!is.null(S)){
    SignaturePrior[, 1:H] <- t(betah)[rep(1, I), ] * S
  }
  colnames(SignaturePrior) <- prior_names

  # Initialize the signatures
  #Alpha <- matrix(alpha, nrow = I, ncol = K)
  if(is.null(R_start)){
    R <- sample_signatures(SignaturePrior)
  } else {
    R <- R_start
  }

  # Initialize the array of total surface of poisson process for easy summation
  X_bar <- array(NA, dim = c(K, J, p))
  for(k in 1:K) X_bar[k, , ] <- X_total

  # Initialize Mu from the prior
  Mu <- matrix(1/rgamma(K * p, 2 * a * J + 1, epsilon *  a * J ), nrow = K, ncol = p)
  # Initialize the loadings
  Betas <- array(rgamma(K * J * p, a, a / Mu), dim = c(K, J, p))

  # Store the output
  MU <- array(NA, c(nsamples,  K, p))
  BETAS <- array(NA, c(nsamples, K, J, p))
  SIGNS <- array(NA, dim = c(nsamples, I, K))
  # Run the Sampling steps
  pb <- txtProgressBar(style=3)
  for(iter in 1:(nsamples + burnin)){
    setTxtProgressBar(pb, iter/(nsamples + burnin))
    # Allocate counts
    #Counts <- sample_Counts(df_trinucl = df_trinucl, R = R, Beta = Betas, X = X)
    Counts <- sample_Counts_cpp(X = X, R = R, Betas = Betas, channel_Id = df_trinucl$Channel_id, sample_Id = df_trinucl$sample_id)
    # Sample the signatures
    Alpha <- Counts$Y_signatures + SignaturePrior
    R <- sample_signatures(Alpha)
    # Sample the Beta loadings
    Betas <- sample_Betas(Counts, Mu, X_bar, a)
    # Sample Mu
    Mu <- sample_Mu(Beta = Betas, epsilon = epsilon, a = a)
    if(iter > burnin) {
      SIGNS[iter - burnin, ,] <- R
      BETAS[iter - burnin, , ,] <- Betas
      MU[iter - burnin, ,] <- Mu
    }
  }
  close(pb)

  return(list("Signatures" =  SIGNS, "Betas" = BETAS, "Mu" = MU))
}


InhomogeneousPoissonNMF_v2 <- function(df_trinucl, R_start = NULL, X, X_total, K = 5, alpha = 0.5, a = 1,
                                       epsilon = 0.001, I = 96, nsamples = 20, burnin = 20){

  # Initialize the sampler
  J <- max(df_trinucl$sample_id)
  p <- ncol(X)

  if(!is.null(R_start)){
    K <- ncol(R_start)
  }

  # Initialize the signatures
  Alpha <- matrix(alpha, nrow = I, ncol = K)
  if(is.null(R_start)){
    R <- sample_signatures(Alpha)
  } else {
    R <- R_start
  }

  # Initialize the array of total surface of poisson process for easy summation
  X_bar <- array(NA, dim = c(K, J, p))
  for(k in 1:K) X_bar[k, , ] <- X_total

  # Initialize Mu from the prior
  Mu <- 1/rgamma(K, 2 * a * J * p + 1, epsilon *  a * J * p)
  Nu <- 1/rgamma(p, 2 * a, a)
  MuNu <- tcrossprod(Mu, Nu)
  # Initialize the loadings
  Betas <- array(rgamma(K * J * p, a, a / MuNu), dim = c(K, J, p))

  # Store the output
  BETAS <- array(NA, c(nsamples, K, J, p))
  SIGNS <- array(NA, dim = c(nsamples, I, K))
  MU <- array(NA, c(nsamples, K))
  NU <- array(NA, c(nsamples, p))
  # Run the Sampling steps
  pb <- txtProgressBar(style=3)
  for(iter in 1:(nsamples + burnin)){
    setTxtProgressBar(pb, iter/(nsamples + burnin))
    # Allocate counts
    #Counts <- sample_Counts(df_trinucl = df_trinucl, R = R, Beta = Betas, X = X)
    Counts <- sample_Counts_cpp(X = X, R = R, Betas = Betas, channel_Id = df_trinucl$Channel_id, sample_Id = df_trinucl$sample_id)
    # Sample the signatures
    Alpha <- Counts$Y_signatures + alpha
    R <- sample_signatures(Alpha)
    # Sample the Beta loadings
    MuNu <- tcrossprod(Mu, Nu)
    Betas <- sample_Betas(Counts, MuNu, X_bar, a)
    # Sample Mu
    Mu <- sample_Mu_v2(Beta = Betas, Nu = Nu, epsilon = epsilon, a = a)
    # Sample Nu
    Nu <- sample_Nu(Beta = Betas, Mu = Mu, epsilon = epsilon, a = a)
    if(iter > burnin) {
      SIGNS[iter - burnin, ,] <- R
      BETAS[iter - burnin, , ,] <- Betas
      MU[iter - burnin, ] <- Mu
      NU[iter - burnin, ] <- Nu
    }
  }
  close(pb)

  return(list("Signatures" =  SIGNS, "Betas" = BETAS, "Mu" = MU, "Nu" = NU))
}


sample_table <- function(probZW){
  sampled_index <- sample(1:length(probZW), size = 1, prob = c(probZW))
  num_rows <- nrow(probZW)
  z <- (sampled_index - 1) %% num_rows + 1
  w <- (sampled_index - 1) %/% num_rows + 1
  return(c("z" = z, "w" = w))
}

sample_Counts <- function(df_trinucl, R, Beta, X) {
  #Index_table <- matrix(NA, ncol = 2, nrow = nrow(df_trinucl))
  K <- ncol(R)
  J <- dim(Beta)[2]
  p <- ncol(X)
  Y_Beta <- array(0, c(K, J, p))/ncol(X)
  Y_signatures <- matrix(0, ncol = K, nrow = nrow(R))
  for(i in 1:nrow(df_trinucl)) {
    j <- df_trinucl$sample_id[i]
    r_i <- R[df_trinucl$Channel_id[i], ]
    probZW <- outer(r_i, X[i, ]) * Beta[, j, ]
    ZW <- sample_table(probZW)
    #Index_table[i, ] <- ZW
    Y_signatures[df_trinucl$Channel_id[i], ZW[1]] <- Y_signatures[df_trinucl$Channel_id[i], ZW[1]] + 1
    Y_Beta[ZW[1], j, ZW[2]] <- Y_Beta[ZW[1], j, ZW[2]] + 1
  }
  return(list("Y_Beta" = Y_Beta, "Y_signatures" = Y_signatures))
}


sample_signatures <- function(Alpha) {
  # Alpha is a matrix of I x K
  I <- nrow(Alpha)
  Ktot <- ncol(Alpha)
  R <- matrix(rgamma(I * Ktot, Alpha), nrow = I) + 1e-10 # Small nugget to avoid degeneracies of the gamma prior
  R <- apply(R, 2, function(x) x/sum(x))
  colnames(R) <- colnames(Alpha)
  return(R)
}

# Sample the Beta coefficients
sample_Betas <- function(Counts, Mu, X_bar, a) {
  K <- dim(X_bar)[1]
  J <- dim(X_bar)[2]
  p <- dim(X_bar)[3]
  a_Mu <- a / Mu
  shapes <- Counts$Y_Beta + a
  rates <- array(NA, dim = c(K, J, p))
  for(j in 1:J) {
    rates[, j, ] <- a_Mu + X_bar[, j, ]
  }
  Betas <- array(rgamma(K * J * p, shapes, rates), dim = c(K, J, p))
  return(Betas)
}

sample_Betas_scaled <- function(Counts, Mu, X_bar, a) {
  K <- dim(X_bar)[1]
  J <- dim(X_bar)[2]
  L <- dim(X_bar)[3]
  a_Mu <- a / Mu
  shapes <- Counts$Y_Beta + a
  rates <- array(NA, dim = c(K, J, L))
  for(j in 1:J) {
    rates[, j, ] <- (a_Mu + 1) * X_bar[, j, ]
  }
  Betas <- array(rgamma(K * J * L, shapes, rates), dim = c(K, J, L))
  return(Betas)
}

# Sample the compressive weight
sample_Mu <- function(Beta, epsilon, a){
  K <- dim(Beta)[1]
  J <- dim(Beta)[2]
  p <- dim(Beta)[3]
  if(length(epsilon) == 1){
    Mu <- matrix(1/rgamma(K * p, 2 * J + 1, epsilon * J + apply(Beta, c(1, 3), sum)),
                 nrow = K, ncol = p)
  } else {
    Mu <- array(NA, c(K, p))
    for(k in 1:K) Mu[k, ] <- 1/rgamma(p, 2 * J + 1, epsilon[k] * J + apply(Beta, c(1, 3), sum))
  }
  return(Mu)
}

sample_Mu_scaled <- function(Betas, X_bar, a, a0, b0){
  K <- dim(Betas)[1]
  J <- dim(Betas)[2]
  L <- dim(Betas)[3]
  totBetas_Xbar <- apply(Betas * X_bar, c(1, 3), sum)
  Mu <- matrix(1/rgamma(K * L,  a * J + a0, b0 +  a * totBetas_Xbar),
               nrow = K, ncol = L)
  return(Mu)
}


# Sample the compressive weight
sample_Mu_v2 <- function(Beta, Nu, epsilon, a){
  K <- dim(Beta)[1]
  J <- dim(Beta)[2]
  p <- dim(Beta)[3]
  b0 <- rowSums(apply(Beta, c(1, 3), sum) / t(Nu)[rep(1, K), ]) + epsilon * a * J * p
  a0 <- 2 * a * J * p + 1
  Mu <- 1/rgamma(K, a0, b0)
  return(Mu)
}


sample_Nu <- function(Beta, Mu, epsilon, a){
  K <- dim(Beta)[1]
  J <- dim(Beta)[2]
  p <- dim(Beta)[3]
  c0 <- colSums(apply(Beta, c(1, 3), sum) / matrix(Mu)[, rep(1, p)]) + 2 * a
  d0 <- a * J * K + a * 0.01
  Nu <- 1/rgamma(p, c0, d0)
  return(Nu)
}


logPost_a <- function(a, epsilon, Betas, Mu_all, X_bar){
  J <- ncol(Betas)
  K <- nrow(Betas)
  L <- dim(Betas)[3]
  XMuBeta <- X_bar / Mu_all * Betas
  eps_j_Mu <- epsilon * J/Mu_all[, 1 ,]
  hyperprior <- sum(J * log(eps_j_Mu) - eps_j_Mu)
  PriorBetas <- sum(log(XMuBeta) - XMuBeta)
  (2 * (J * K * L) * a + (K * L)) * log(a) - (J * K * L) * lgamma(a) -
    K * L * lgamma(a * J + 1) + a * (hyperprior + PriorBetas) + 99 * log(a) - 100 * a
}

sample_Post_a <- function(nsamples, epsilon, Betas, Mu_all, X_bar) {
  UppU <- -nlminb(1, objective = function(x) -logPost_a(x, epsilon, Betas, Mu_all, X_bar),
                  lower = 1e-6)$objective

  UppV <- -nlminb(1, objective = function(x) -logPost_a(x, epsilon, Betas, Mu_all, X_bar) - 2 * log(x),
                  lower = 1e-6)$objective

  s <- exp(0.5 * (UppV - UppU))
  samples <- rep(NA, nsamples)
  for(n in 1:nsamples){
    accepted <- 0
    while(accepted == 0) {
      u <- runif(1); v <- runif(1)
      x <- s * v / u
      if(log(u) < 0.5 * (logPost_a(x, epsilon, Betas, Mu_all, X_bar) - UppU)){
        samples[n] = x
        accepted = 1
      }
    }
  }
  samples
}


simulation_dir <- "~/CompressiveNMF/output/Compressive_vs_InfiniteFactors/"
data_all <- readRDS("~/CompressiveNMF/output/Compressive_vs_InfiniteFactors/Scenario_100_overdisp_0_Knew_6_theta_100/data.rds.gzip")
data <- data_all[[1]]

X <- CompressiveNMF::Mutations_21breast
library(CompressiveNMF)
set.seed(42)
res <- CompressiveNMF_randa(data$X, K = 15, nsamples = 1000,
                            random_a = TRUE, c0 = 100, d0 = 100,
                            burnin = 5000, cutoff_excluded = 0)
plot_SBS_signature(res$Signatures[, res$RelWeights > 0.05])
plot(res$mcmc_out[[1]]$Mu[, "Sig_new15"])

round(res$RelWeights, 5)
round(res$RelWeights[res$RelWeights>0.05], 5)
plot(res$mcmc_out[[1]]$a, type = "l")


Mu <- res$RelWeights
Theta <- res$Weights

a <- 1.2
logPost_a <-  function(x, Mu, Theta, epsilon, c0, d0){
  K <- length(Mu)
  J <- ncol(Theta)
  K * J * x * log(x) + (K * J * x + K)*log(epsilon*x*J) - K * J * lgamma(x) - K * lgamma(x * J + 1) -
   x * (2 * J * sum(log(Mu)) - sum(log(Theta)) + sum(rowSums(Theta)/Mu) +
          epsilon * J * sum(1/Mu)) + (c0 - 1) * log(x) - x * d0
}

logPost_a(1.2, Mu, Theta, epsilon = 0.001, c0 = 1, d0 = 1)

sample_Post_a <- function(nsamples, Mu, Theta, epsilon, c0, d0,  a_start = 1) {
  UppU <- -nlminb(start = a_start, objective = function(x) -logPost_a(x, Mu, Theta, epsilon, c0, d0),
                  lower = 1e-6)$objective

  UppV <- -nlminb(start = a_start, objective = function(x) -logPost_a(x, Mu, Theta, epsilon, c0, d0) - 2 * log(x),
                  lower = 1e-6)$objective

  s <- exp(0.5 * (UppV - UppU))
  samples <- rep(NA, nsamples)
  for(n in 1:nsamples){
    accepted <- 0
    while(accepted == 0) {
      u <- runif(1); v <- runif(1)
      x <- s * v / u
      if(log(u) < 0.5 * (logPost_a(x, Mu, Theta, epsilon, c0, d0) - UppU)){
        samples[n] = x
        accepted = 1
      }
    }
  }
  samples
}

epsilon <- 0.001
c0 <- 1
d0 <- 1
my_a <- sample_Post_a(5000, Mu, Theta, epsilon, c0, d0)

library(rstan)

stan_code <- "
data {
  int<lower=1> K;
  int<lower=1> J;
  real<lower=0> c0;
  real<lower=0> d0;
  real<lower=0> eps;
  vector<lower=0>[K] mu;
  matrix<lower=0>[K, J] theta;
}
parameters {
  real<lower=0> a;
}
model {
  a ~ gamma(c0, d0);
  mu ~ inv_gamma(a * J + 1, eps * a * J);
  for (k in 1:K)
    for (j in 1:J)
      theta[k,j] ~ gamma(a, a / mu[k]);
}
"

stan_data <- list(
  K = nrow(Theta),
  J = ncol(Theta),
  c0 = 1,
  d0 = 1,
  eps = 0.001,
  mu = Mu,
  theta = Theta
)

fit <- stan(
  model_code = stan_code,
  data = stan_data,
  seed = 123,
  chains = 4,
  iter = 2000,
  warmup = 1000
)

ext <- extract(fit)

print(fit, pars = "a")
# posterior samples of a:
a_draws <- extract(fit, pars = "a")$a
plot(density(a_draws))
lines(density(my_a), col = "red")
qqplot(a_draws, my_a)

rAlpha_ru <- function(nsamples = 1, Mu, Theta, epsilon = 0.001, c0 = 1, d0 = 1){
  x2 <- suppressWarnings(rust::ru(logf = logPost_a, d = 1, n = nsamples,
                                  lower = 1e-12, init = 1,
                                  Mu = Mu, Theta = Theta, epsilon = epsilon,
                                  c0 = c0, d0 = d0))
  return(c(x2$sim_vals))
}



