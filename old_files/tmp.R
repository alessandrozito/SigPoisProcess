# Load the package
library(MCMCpack)
library(TruncatedNormal)

# Set parameters
v <- 25                   # degrees of freedom (> dimension - 1)
Sigma <- matrix(c(1,0.25,0.25,1), 2, 2)  # Scale matrix (must be positive-definite)

# Simulate one sample
IW_sample <- riwish(v, Sigma)
IW_sample

# To generate multiple samples:
n <- 10
IW_samples <- replicate(n, riwish(v, Sigma), simplify="array") # (2 x 2 x n) array

S <- IW_samples[, ,4]
pmvnorm(mu = rep(0, ncol(Sigma)), sigma = S, lb = rep(0, ncol(Sigma)))
tt <-rtmvnorm(1000, rep(0, ncol(Sigma)), sigma = S)
plot(tt)
abline(h = 0, v = 0)
mean(rowSums(tt>0) == 2)



rho <- 0.2
s1 <- 300 * pi/2
s2 <- 0.01 * pi/2
S <- matrix(c(s1, sqrt(s1) * sqrt(s2) * rho,
              sqrt(s1) * sqrt(s2) * rho, s2), 2, 2)
colMeans(rtmvnorm(10000, rep(0, ncol(Sigma)), sigma = S, lb = rep(0, ncol(Sigma))))
sqrt(diag(S)) * sqrt(2/pi)

(1 + rho) * sqrt(diag(S)) / (sqrt(2 * pi) *
  (pmvnorm(mu = rep(0, ncol(Sigma)), sigma = matrix(c(1,rho,rho, 1), 2,2), lb = rep(0, ncol(Sigma)))))

S <- riwish(10, diag(5) * 100)

J <- 100
eps <- 0.01
S <- riwish(J + 5 + 1, diag(10) * J * eps^2 * pi / 2)

replicate(1000, riwish(J + 5 + 1, diag(5) * J * eps^2 * pi / 2))
apply(replicate(1000, riwish(J + 5 + 1, diag(5) * J * eps^2 * pi / 2)), c(1,2), mean)

tt <- rtmvnorm(10000, rep(0, ncol(S)), sigma = S, lb = rep(0, ncol(S)))
crossprod(tt)/10000

sqrt(diag(S)) * sqrt(2/pi)
cov2cor(S)

hist(S - t(S))


riwish(10, matrix(c(1, 0.3, 0.5,1), 2, 2))
