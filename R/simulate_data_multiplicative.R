library(tidyverse)
library(maftools)
library(BSgenome)
library(sigminer)
library(Biostrings)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(foreach)
library(doParallel)
library(corrplot)
devtools::document()

match_MutSign <- function(R_true, R_hat) {
  # Need to make sure that R_hat and R_true are the same
  k_true <- ncol(R_true)
  k_hat <- ncol(R_hat)
  k_tot <- max(c(k_hat, k_true))
  I <- nrow(R_true)
  mat0 <- matrix(100, nrow = I, ncol = abs(k_hat - k_true)) #
  if (k_hat > k_true) {
    colnames(mat0) <- paste0("new_extra", 1:ncol(mat0))
    R_true <- cbind(R_true, mat0)
  } else if (k_hat < k_true) {
    R_hat <- cbind(R_hat, mat0)
  }

  # Match mutational signatures using the Hungarian algorithm
  CosMat <- matrix(1, k_tot, k_tot)
  for (i in 1:k_tot) {
    for (j in 1:k_tot) {
      CosMat[i, j] <- 1 - cosine(R_true[, i], R_hat[, j])
    }
  }
  match <- RcppHungarian::HungarianSolver(CosMat)$pairs[, 2]
  R_hat_matched <- R_hat[, match]

  # Change 100 with 0s
  R_hat_matched[R_hat_matched == 100] <- 0
  R_true[R_true == 100] <- 0

  return(list("R_hat" = R_hat_matched, "R_true" = R_true, "match" = match))
}

generate_SignalTrack <- function(length_genome = 2e6, tilewidth = 100,
                                 rho = 0.99, p = 10, Sigma = diag(p)){
  # Generate the genome and the splits
  grX <- tileGenome(c("chrsim" = length_genome), tilewidth = tilewidth,
                    cut.last.tile.in.chrom = TRUE)
  # Generate the covariates from random starting points
  X <- matrix(NA, nrow = length(grX), ncol = p)
  D <- diag(rep(rho, p))
  X[1, ] <- mvtnorm::rmvnorm(1, sigma = Sigma)
  for(i in 2:nrow(X)){
    eps <- c(mvtnorm::rmvnorm(1, sigma = Sigma))
    X[i, ] <- c(D %*% X[i-1, ]) + c(eps)
  }
  # Make them positive
  #X <- exp(scale(X))
  #X <- -log10(pnorm(scale(X)))
  #X <- apply(X, 2, function(x) (x - min(x))/(max(x) - min(x)))
  #X2 <- apply(X, 2, function(x) (x - min(x))/(max(x) - min(x)))
  colnames(X) <- paste0("Encode", sprintf("%02d", 1:ncol(X)))
  # Append and return
  grX$bin_weight <- tilewidth
  mcols(grX) <- cbind(mcols(grX), scale(X))
  return(grX)
}

generate_Parameters <- function(cosmic_sigs = c("SBS1", "SBS2", "SBS13", "SBS3"),
                                K_new = 4,
                                J = 50,
                                a = 1,
                                p = 5,
                                theta = 200,
                                prob_zero = 0,
                                sd_beta = 1/3,
                                alpha = 0.1){

  Rmat_cos <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37[, cosmic_sigs]
  if(K_new > 0){
    Rmat_random <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, 96)))
    colnames(Rmat_random) <- paste0("SBSnew", 1:K_new)
  }
  Rmat <- cbind(Rmat_cos, Rmat_random)

  # Generate Theta
  K_true <- ncol(Rmat)
  exposures <- rgamma(K_true, theta, 1)
  Theta <- matrix(rgamma(K_true * J, 0.5, 0.5), ncol = J, nrow = K_true)
  Theta <- apply(Theta, 2, function(x) x * exposures)
  rownames(Theta) <- colnames(Rmat)

  # Generate Betas
  Betas <- matrix(rnorm(K_true * p, mean = 0, sd = sd_beta), nrow = p, ncol = K_true)
  Null <- 1 - matrix(rbinom(n = length(Betas), size = 1, prob = prob_zero),
         nrow = nrow(Betas), ncol = ncol(Betas))
  if(any(colSums(Null) == 0)) {
    id <- sample(1:p, size = 1) # add one covariate at random
    Null[id, colSums(Null) == 0] <- 1
  }
  Betas <- Betas * Null # Multiply to generate the mixture
  return(list(R = Rmat, Theta = Theta, Betas = Betas))
}

# p_all <- 10
# p_to_use <- 5
# grX <- generate_SignalTrack(length_genome = 2e6,
#                             tilewidth = 100, rho = 0.995, p = p_all,
#                             Sigma = diag(p_all))
# pars <- generate_Parameters(p = p_to_use, sd_beta = 1/3, theta = 200)
# Xcov <- as.matrix(mcols(grX)[, -1][, 1:p_to_use])
# #
# plot(Xcov[, 1], type = "l")
# lines(Xcov[, 2], col = "red")
# lines(Xcov[, 3], col = "green")
# #
# # #Xcov <- -log10(pnorm(log(Xcov)))
# CumSums <- apply(exp(Xcov %*% pars$Betas), 2, function(x) cumsum(width(grX) * x))
# i <- sample(1:96, size = 1)
# j <- sample(1:40, size = 1)
# Theta_adj <- pars$Theta/sum(width(grX))
# mm <- exp(Xcov %*% pars$Betas)
#
#
# Lambda_Tmax <- c(crossprod(pars$R[i, ], Theta_adj[, j] * c(tail(CumSums, 1))))
# Lambdas <- c(crossprod(pars$R[i, ], Theta_adj[, j] * t(CumSums)))
# N_Tmax <- rpois(1, Lambda_Tmax)
# u <- sort(runif(n = N_Tmax))
# ts <- findInterval(u, Lambdas/Lambda_Tmax) * 100 + sample(1:(tilewidth-1), size = length(u), replace = TRUE)
# id <- ts[1]
#   rg <- RangesX[id]#ranges(grX)[id]
#   t <- sample(start(rg):end(rg), size = 1)

# plot(Lambdas, type = "l")
# Lambda_Tmax
#
#
# LambdaMat <- matrix(NA, nrow = 96, ncol = 50)
# for(i in 1:96){
#  for(j in 1:50){
#    LambdaMat[i, j] <- c(crossprod(pars$R[i, ], Theta_adj[, j] * c(tail(CumSums, 1))))
#  }
# }
# hist(LambdaMat)



sample_mutations <- function(i, j, R, Theta, RangesX, CumSums, tilewidth = 100){

  #Lambda_Tmax <- c(crossprod(R[i, ], Theta[, j] %*% c(tail(CumSums, 1))))
  Lambda_Tmax <- c(crossprod(R[i, ], Theta[, j] * c(tail(CumSums, 1))))
  Lambdas <- c(crossprod(R[i, ], Theta[, j] * t(CumSums)))
  N_Tmax <- rpois(1, Lambda_Tmax)
  #print(N_Tmax)
  if(N_Tmax > 0){
    # Sample from the uniform
    u <- sort(runif(n = N_Tmax))
    ts <- findInterval(u, Lambdas/Lambda_Tmax) * tilewidth + sample(1:(tilewidth - 1), size = length(u), replace = TRUE)
    # Invert to find location
    # ts <- sapply(u, function(uu) {
    #   # ids <- which(Lambdas <= uu * Lambda_Tmax)
    #   # if(length(ids)==0){
    #   #   id <- 1
    #   # } else {
    #   #   id <- max(ids)
    #   # }
    #   # rg <- RangesX[id]#ranges(grX)[id]
    #   # t <- sample(start(rg):end(rg), size = 1)
    #   # t
    #   ids
    # })
  } else {
    ts <- NA
  }
  return(ts)
}


sample_dataset <- function(R, Theta, Betas, RangesX, CumSums, ncores = 20){
  # Model dimensions
  I <- nrow(R)
  J <- ncol(Theta)
  Mutations <- data.frame()

  #registerDoParallel(ncores)
  for(j in 1:J) {
    print(paste0("Simulating Patient_", sprintf("%02d", j)))
    for(i in 1:I) {
      ch <- rownames(R)[i]
      ts <- sample_mutations(i, j, R, Theta, RangesX, CumSums)
      if(!is.na(ts[1])) {
        Mutations <- rbind(Mutations, data.frame(sample = paste0("Patient_", sprintf("%02d", j)),
                                                 pos = ts,
                                                 channel = ch))
        # data.frame(sample = paste0("Patient_", sprintf("%02d", j)),
        #            pos = ts,
        #            channel = ch)
      }
    }
  }
  Mutations$sample <- factor(Mutations$sample, levels = unique(Mutations$sample))
  Mutations$channel <- as.factor(Mutations$channel)
  return(Mutations)
}

sample_dataset_parallel <- function(R, Theta, Betas, RangesX, CumSums, ncores = 20){
  # Model dimensions
  I <- nrow(R)
  J <- ncol(Theta)
  #Mutations <- data.frame()

  registerDoParallel(ncores)
  Mutations <- foreach(j = 1:J, .combine = "rbind") %dopar% {
    print(paste0("Simulating Patient_", sprintf("%02d", j)))
    mut_temp <- data.frame()
    for(i in 1:I) {
      ch <- rownames(R)[i]
      ts <- sample_mutations(i, j, R, Theta, RangesX, CumSums)
      if(!is.na(ts[1])) {
        mut_temp <- rbind(mut_temp, data.frame(sample = paste0("Patient_", sprintf("%02d", j)),
                                                 pos = ts,
                                                 channel = ch))
        # data.frame(sample = paste0("Patient_", sprintf("%02d", j)),
        #            pos = ts,
        #            channel = ch)
      }
    }
    mut_temp
  }
  Mutations$sample <- factor(Mutations$sample, levels = unique(Mutations$sample))
  Mutations$channel <- as.factor(Mutations$channel)
  return(Mutations)
}

generate_MutationData <- function(J = 50,
                                  cosmic_sigs = c("SBS1", "SBS2", "SBS13", "SBS3"),
                                  K_new = 4,
                                  theta = 200,
                                  ngaps = 10,
                                  a = 1,
                                  sd_beta = 1/3,
                                  length_genome = 2e6,
                                  tilewidth = 100,
                                  rho = 0.99,
                                  p_all = 10,
                                  p_to_use = 5,
                                  corr = "indep"){

  # Step 1 - generate a covariance matrix for the
  if (corr == "indep") {
    Sigma = diag(p_all)
  } else if (corr == "onion"){
    corrmat <- clusterGeneration::genPositiveDefMat(p_all, covMethod="onion")$Sigma
    Sigma <- cov2cor(corrmat)
  }

  # Step 2 - Generate the SignalTrack, and calculate the Cumulative Sums
  grX <- generate_SignalTrack(length_genome = length_genome,
                              tilewidth = tilewidth, rho = rho, p = p_all,
                              Sigma = Sigma)

  # Step 3 - Sample the model parameters
  pars_true <- generate_Parameters(cosmic_sigs = cosmic_sigs,
                                   K_new = K_new, J = J,
                                   a = a,
                                   p = p_to_use,
                                   theta = theta, sd_beta = sd_beta)
  # Step 4 - Rescale the Betas
  Xcovs <- mcols(grX)[, -1] # Remove bin_weights
  bin_weight = grX$bin_weight
  CumSums <- apply(exp(as.matrix(Xcovs[, 1:p_to_use]) %*% pars_true$Betas),
                   2, function(x) cumsum(bin_weight * x))
  #X_totals <- CumSums[rep(nrow(CumSums), J), ]
  #X_bar <- array(NA, dim = c(ncol(pars_true$R), J, p_to_use))
  #for(k in 1:nrow(X_bar)) X_bar[k, ,] <- X_totals
  Theta_scaled <- pars_true$Theta / sum(bin_weight)
  # Step 5 - Generate the dataset of mutations
  Mutations <- sample_dataset_parallel(R = pars_true$R,
                                      Betas = pars_true$Betas,
                                      Theta = Theta_scaled,
                                      RangesX = ranges(grX),
                                      CumSums = CumSums)

  # Step 6 - Normalize the covariates to have mean 1, and then merge with covariates
  gr_Mutations <- GRanges(seqnames = "chrsim",
                          IRanges(start = Mutations$pos, end = Mutations$pos),
                          sample = Mutations$sample,
                          channel = Mutations$channel)
  overlaps <- findOverlaps(gr_Mutations, grX)
  mcols(gr_Mutations) <- cbind(mcols(gr_Mutations), mcols(grX[subjectHits(overlaps)]))
  gr_Mutations$bin_weight <- NULL

  data <- list("gr_Mutations" = gr_Mutations,
               "covariate_used" = colnames(mcols(grX)[, -1])[1:p_to_use],
               "R" = pars_true$R,
               "Betas" = pars_true$Betas,
               "Theta" = Theta_scaled,
               "gr_SignalTrack" = grX,
               "simulation_parameters" = list(J = J,
                                              cosmic_sigs = cosmic_sigs,
                                              K_new = K_new,
                                              theta = theta,
                                              a = a,
                                              length_genome = length_genome,
                                              tilewidth = tilewidth,
                                              rho = rho,
                                              p_all = p_all,
                                              p_to_use = p_to_use,
                                              corr = corr))
  return(data)

}

#------- Generate data
set.seed(24, kind = "L'Ecuyer-CMRG")
data <- generate_MutationData(J = 50, sd_beta = 1/2)
saveRDS(data, "~/SigPoisProcess/data/temp_simulation_multiplicative.rds.gzip", compress="gzip")
data <- readRDS("~/SigPoisProcess/data/temp_simulation_multiplicative.rds.gzip")

#---- Plot the number of mutation in each region
length(data$gr_Mutations)
corrplot(cor(as.matrix(mcols(data$gr_SignalTrack))[, -1]))
corrplot(cor(as.matrix(mcols(data$gr_Mutations))[, -c(1,2)]))

#-- Intensity function
grX <- data$gr_SignalTrack
p_all <- 10
p_to_use <- 5
Xcov <- as.matrix(mcols(grX)[, -1][, 1:p_to_use])
ExpTrack <- apply(exp(Xcov %*% data$Betas), 2, function(x) width(grX) * x)
LambdaAll <- rep(0, length(grX))
for(i in 1:96){
 for(j in 1:ncol(data$Theta)){
   LambdaAll <- LambdaAll +  c(crossprod(data$R[i, ], data$Theta[, j] * t(ExpTrack)))
 }
}
plot(LambdaAll, type = "l")

#-- Number of mutation in each bin
df_tmp <- as.data.frame(data$gr_Mutations)
over <- findOverlaps(data$gr_Mutations, data$gr_SignalTrack)
df_tmp$region <- subjectHits(over)
df_all <- df_tmp %>%
  group_by(region)%>%
  summarize(n = n())
MutAll <- rep(0, length(grX))
MutAll[df_all$region] <- df_all$n
plot(MutAll, type = "l")

plot(MutAll, type = "l")
lines(LambdaAll, col = "red")


mat <- getTotalMutations(data$gr_Mutations)

# Model 1 ---------- Try the model with the full set of covariates
set.seed(10)
gr_Mutations <- data$gr_Mutations
SignalTrack <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
bin_weigth <- data$gr_SignalTrack$bin_weight

#--------------
SignalTrack <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
#SignalTrack_red <- SignalTrack#[, -c(1:5)]
#mcols(gr_Mutations_red) <- mcols(gr_Mutations_red)[, c(1, 2, 8:12)]
bin_weigth <- data$gr_SignalTrack$bin_weight

set.seed(10)
resFull <- SigPoisProcess_mult(gr_Mutations = gr_Mutations,
                           SignalTrack = SignalTrack,
                           bin_weight = bin_weigth,
                           K = 10,
                           method = "map",
                           update_Signatures = TRUE,
                           update_Theta = TRUE,
                           update_Betas = TRUE,
                           #init = SigPoisProcess_mult.init(R_start = resSol$Signatures),
                           controls = SigPoisProcess_mult.controls(maxiter = 1000, tol = 1e-6),
                           prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                          c0 =  100, d0 = 1))
saveRDS(resFull, "~/SigPoisProcess/output/Simulation_full_nocorr.rds")

p1 <- CompressiveNMF::plot_SBS_signature(resFull$Signatures)
p2 <- plot_betas(resFull$Betas)
p1 + p2

match <- match_MutSign(R_true =  data$R, R_hat = resFull$Signatures)
plot(resFull$sol$Betas[, match$match][1:5, 1:8], data$Betas)
abline(a = 0, b = 1)

par(mfrow = c(2,1))
plot(resFull$sol$trace[-c(1:3)])
plot(resFull$sol$MUchain[1:resFull$sol$iter, 1][seq(1, resFull$sol$iter, 10)][-c(1:3)], type = "l")
plot(resFull$sol$MUchain[1:resFull$sol$iter, 6][seq(1, resFull$sol$iter, 10)][-c(1:3)], type = "l")
plot(resFull$sol$Sigma2[1:resFull$sol$iter, 6][seq(1, resFull$sol$iter, 10)][-c(1:3)], type = "l")

#----- Total mutation rate
ExpTrackPred <- apply(exp(SignalTrack %*% resFull$Betas), 2, function(x) width(grX) * x)
LambdaPred<- rep(0, length(grX))
for(i in 1:96){
  for(j in 1:ncol(resFull$Theta)){
    LambdaPred <- LambdaPred +  c(crossprod(resFull$Signatures[i, ], resFull$Theta[, j] * t(ExpTrackPred)))
  }
}
plot(MutAll, type = "l", col = "gray")
lines(LambdaPred, col = "red", lty = "dashed")

estim


# Model 2 ---------- Remove all true covariates, throw in the redundant ones
set.seed(10)
gr_Mutations_red <- data$gr_Mutations
SignalTrack_red <- SignalTrack[, -c(1:5)]
mcols(gr_Mutations_red) <- mcols(gr_Mutations_red)[, c(1, 2, 8:12)]
bin_weigth <- data$gr_SignalTrack$bin_weight

set.seed(10)
resNoise <- SigPoisProcess_mult(gr_Mutations = gr_Mutations_red,
                               SignalTrack = SignalTrack_red,
                               bin_weight = bin_weigth,
                               K = 10,
                               method = "map",
                               update_Signatures = TRUE,
                               update_Theta = TRUE,
                               update_Betas = TRUE,
                               #init = SigPoisProcess_mult.init(R_start = resSol$Signatures),
                               controls = SigPoisProcess_mult.controls(maxiter = 1000, tol = 5e-7),
                               prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                              c0 =  100, d0 = 1))

p3 <- CompressiveNMF::plot_SBS_signature(resNoise$Signatures)
p4 <- plot_betas(resNoise$Betas)
p3 + p4

resNoise$Betas
plot_betas(resNoise$Betas[, -c(4, 8, 9)])
rowMeans(resNoise$Betas[, c(4, 8, 9)])

plot(resNoise$sol$trace[-c(1:3)])
colnames(resNoise$Mu) <- colnames(resNoise$Signatures)

SigEst <- resNoise$Signatures[, -c(2,4)]
match <- match_MutSign(R_true =  data$R, R_hat = SigEst)
plot(SigEst[, match$match][, 1:8], data$R)
abline(a = 0, b = 1)

#------ Prediction under bad value
ExpTrackPred <- apply(exp(SignalTrack_red %*% resNoise$Betas), 2, function(x) width(grX) * x)
LambdaPred<- rep(0, length(grX))
for(i in 1:96){
  for(j in 1:ncol(resFull$Theta)){
    LambdaPred <- LambdaPred +  c(crossprod(resNoise$Signatures[i, ], resNoise$Theta[, j] * t(ExpTrackPred)))
  }
}
plot(MutAll, type = "l", col = "gray")
lines(LambdaPred, col = "red", lty = "dashed")
plot(LambdaPred, MutAll)


#----------- Generate now correlated covariates
set.seed(24, kind = "L'Ecuyer-CMRG")
data_corr <- generate_MutationData(J = 50, sd_beta = 1/2, corr = "onion")
saveRDS(data_corr, "~/SigPoisProcess/data/temp_simulation_multiplicative_corr.rds.gzip", compress="gzip")
#data <- readRDS("~/SigPoisProcess/data/temp_simulation_multiplicative.rds.gzip")

length(data_corr$gr_Mutations)
corrplot(cor(as.matrix(mcols(data_corr$gr_SignalTrack))[, -1]))
corrplot(cor(as.matrix(mcols(data_corr$gr_Mutations))[, -c(1,2)]))

df_tmp <- as.data.frame(data_corr$gr_Mutations)
over <- findOverlaps(data_corr$gr_Mutations, data_corr$gr_SignalTrack)
df_tmp$region <- subjectHits(over)
df_all <- df_tmp %>%
  group_by(region)%>%
  summarize(n = n())
MutAll <- rep(0, length(grX))
MutAll[df_all$region] <- df_all$n
plot(MutAll, type = "l")

plot(MutAll, type = "l")
lines(LambdaAll, col = "red")


#--------------
gr_Mutations_corr <- data_corr$gr_Mutations
SignalTrack_corr <- as.matrix(mcols(data_corr$gr_SignalTrack))[, -1]
#SignalTrack_red <- SignalTrack#[, -c(1:5)]
#mcols(gr_Mutations_red) <- mcols(gr_Mutations_red)[, c(1, 2, 8:12)]
bin_weigth <- data_corr$gr_SignalTrack$bin_weight

set.seed(10)
resFull_corr <- SigPoisProcess_mult(gr_Mutations = gr_Mutations_corr,
                               SignalTrack = SignalTrack_corr,
                               bin_weight = bin_weigth,
                               K = 10,
                               method = "map",
                               update_Signatures = TRUE,
                               update_Theta = TRUE,
                               update_Betas = TRUE,
                               #init = SigPoisProcess_mult.init(R_start = resSol$Signatures),
                               controls = SigPoisProcess_mult.controls(maxiter = 1000, tol = 1e-6),
                               prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                              c0 =  100, d0 = 1))


p1 <- CompressiveNMF::plot_SBS_signature(resFull_corr$Signatures)
p2 <- plot_betas(resFull_corr$Betas)
p1 + p2

ExpTrackPred <- apply(exp(SignalTrack_corr %*% resFull_corr$Betas), 2, function(x) width(grX) * x)
LambdaPred<- rep(0, length(grX))
for(i in 1:96){
  for(j in 1:ncol(resFull$Theta)){
    LambdaPred <- LambdaPred +  c(crossprod(resFull_corr$Signatures[i, ], resFull_corr$Theta[, j] * t(ExpTrackPred)))
  }
}
plot(MutAll, type = "l", col = "gray")
lines(LambdaPred, col = "red", lty = "dashed")
plot(LambdaPred, MutAll)


match <- match_MutSign(R_true =  data_corr$R, R_hat = resFull_corr$Signatures)
plot(resFull_corr$sol$Betas[, match$match][1:5, 1:8], data_corr$Betas)
abline(a = 0, b = 1)

plot(resFull_corr$sol$trace[-c(1,4)])


# Model 4 ---------- Remove all true covariates, throw in the redundant ones
set.seed(10)
gr_Mutations_red_corr <- data_corr$gr_Mutations
SignalTrack_red_corr <- SignalTrack_corr[, -c(1:5)]
mcols(gr_Mutations_red_corr) <- mcols(gr_Mutations_red_corr)[, c(1, 2, 8:12)]
bin_weigth <- data$gr_SignalTrack$bin_weight

set.seed(10)
resNoise_corr <- SigPoisProcess_mult(gr_Mutations = gr_Mutations_red_corr,
                                SignalTrack = SignalTrack_red_corr,
                                bin_weight = bin_weigth,
                                K = 10,
                                method = "map",
                                update_Signatures = TRUE,
                                update_Theta = TRUE,
                                update_Betas = TRUE,
                                #init = SigPoisProcess_mult.init(R_start = resSol$Signatures),
                                controls = SigPoisProcess_mult.controls(maxiter = 1000, tol = 1e-6),
                                prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                               c0 =  100, d0 = 1))

p3 <- CompressiveNMF::plot_SBS_signature(resNoise_corr$Signatures)
p4 <- plot_betas(resNoise_corr$Betas)
p3 + p4

resFull_corr$Mu

ExpTrackPred <- apply(exp(SignalTrack_red_corr %*% resNoise_corr$Betas), 2, function(x) width(grX) * x)
LambdaPred<- rep(0, length(grX))
for(i in 1:96){
  for(j in 1:ncol(resFull$Theta)){
    LambdaPred <- LambdaPred +  c(crossprod(resNoise_corr$Signatures[i, ], resNoise_corr$Theta[, j] * t(ExpTrackPred)))
  }
}
plot(MutAll, type = "l", col = "gray")
lines(LambdaPred, col = "red", lty = "dashed")
plot(LambdaPred, MutAll)

plot(resNoise_corr$sol$trace[-c(1:3)])
#-------------------------------------------------------------------------------


set.seed(10)
gr_Mutations <- gr_Mutations_red <- data$gr_Mutations

SignalTrack <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
SignalTrack_red <- SignalTrack[, -c(1:5)]
mcols(gr_Mutations_red) <- mcols(gr_Mutations_red)[, c(1, 2, 8:12)]
bin_weigth <- data$gr_SignalTrack$bin_weight

res2 <- SigPoisProcess_mult(gr_Mutations = gr_Mutations_red,
                           SignalTrack = SignalTrack_red,
                           bin_weight = bin_weigth,
                           K = 12,
                           method = "map",
                           update_Signatures = TRUE,
                           update_Theta = TRUE,
                           update_Betas = TRUE,
                           #init = SigPoisProcess_mult.init(R_start = data$R),
                           controls = SigPoisProcess_mult.controls(maxiter = 1000,
                                                                   tol = 1e-7),
                           prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "none",
                                                                          c0 =  10, d0 = 1))

res <- SigPoisProcess_mult(gr_Mutations = gr_Mutations_red,
                           SignalTrack = SignalTrack_red,
                           bin_weight = bin_weigth,
                           K = 10,
                           method = "map",
                           update_Signatures = TRUE,
                           update_Theta = TRUE,
                           update_Betas = TRUE,
                           #init = SigPoisProcess_mult.init(R_start = resSol$Signatures),
                           controls = SigPoisProcess_mult.controls(maxiter = 1000, tol = 1e-6),
                           prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                          c0 =  100, d0 = 1))

resSol <- CompressiveNMF::CompressiveNMF_map(mat)
CompressiveNMF::plot_SBS_signature(resSol$Signatures)


plot(res$sol$BETASchain[, 2, 1])
plot(res$sol$trace[-c(1:3)])

plot(res$sol$Sigma2, res$sol$Mu)
plot(res$Betas[1:5, ], data$Betas)
abline(a = 0, b = 1)
plot(res$Theta, data$Theta)
abline(a = 0, b = 1)
round(res$Thetas[5, ] * sum(bin_weigth), 4)
res$Betas[, 5]
p1 <- CompressiveNMF::plot_SBS_signature(res$Signatures)
p2 <- plot_betas(res$Betas)
p1 + p2

plot_betas(res$Betas[, -c(1,4,6,8)])
round(res$sol$Lambda, 3)

match <- match_MutSign(R_true =  data$R, R_hat = res$Signatures)
plot(res$sol$Betas[, match$match][1:5, 1:8], data$Betas)
abline(a = 0, b = 1)

j <- 3
plot(res$sol$Betas[, match$match][1:5, j], data$Betas[, j])
abline(a = 0, b = 1)

rowSums(res$sol$Betas^2)

plot(res$Signatures[, match$match[1:8]], data$R)
abline(a = 0, b = 1)

plot(res$Thetas[match$match, ], data$Theta)
abline(a = 0, b = 1)

plot(match$R_hat, match$R_true)
abline(a = 0, b = 1)

plot(res$Thetas[match$match[1:8], ], data$Theta)
abline(a = 0, b = 1)


# Now, I want to predict the mutation rate.
ExpTrackPred <- apply(exp(SignalTrack_red %*% res$Betas), 2, function(x) width(grX) * x)
LambdaPred<- rep(0, length(grX))
for(i in 1:96){
  for(j in 1:ncol(res$Theta)){
    LambdaPred <- LambdaPred +  c(crossprod(res$Signatures[i, ], res$Theta[, j] * t(ExpTrackPred)))
  }
}
plot(LambdaAll, type = "l")
lines(LambdaPred, col = "red")

mm <- SigPoisProcess::estimateTotalCounts.SigPoisProcess_mult(res,
                                                              SignalTrack_red,
                                                              width(grX))
plot(mat, mm)

#----- Let's try to fit the model without the data. It should go to the prior.
gr_Mutations_prior <- gr_Mutations_red[c(1, 4363)]
gr_Mutations_prior$sample <- as.factor(as.numeric(gr_Mutations_prior$sample))
res <- SigPoisProcess_mult(gr_Mutations = gr_Mutations_red,
                           SignalTrack = SignalTrack_red,
                           bin_weight = bin_weigth,
                           K = 10,
                           method = "map",
                           update_Signatures = TRUE,
                           update_Theta = TRUE,
                           update_Betas = TRUE,
                           #init = SigPoisProcess_mult.init(R_start = data$R),
                           controls = SigPoisProcess_mult.controls(maxiter = 500,
                                                                   tol = 1e-6),
                           prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "none",
                                                                          c0 =  100, d0 = 1))

plot(res$sol$trace[-c(1:15)], type = "l")
plot(res$sol$trace)

plot(res$sol$trace_logLik[-c(1:10)], type = "l")
plot(res$sol$trace_logLik, type = "l")
plot(res$sol$trace_logPrior[-c(1:10)], type = "l")
plot(res$sol$trace, type = "l")

plot(round(res$sol$BETASchain[, 2,1], 6), type = "l")

plot(res$sol$BETASchain[-c(1:10), 1,9])

plot(res$sol$MUchain[-c(1:10), 1])

plot(res$sol$THETAchain[-c(1:10), 1,1])
plot(res$sol$THETAchain[-c(1:10), 10,2], type = "l")

plot(res$sol$SIGSchain[-c(1:10), 10,2])


plot_betas(res$Betas)
res$Mu
res$Thetas* sum(bin_weigth)
CompressiveNMF::plot_SBS_signature(res$Signatures)
res$Mu
plot_betas(res$Betas)

# How well do we retrieve the mutation matrix?
res3 <- res
res3$Signatures <- res3$Signatures[, c(res3$Mu) > 0.01]
res3$Thetas <- res3$Thetas[c(res3$Mu) > 0.01, ]
res3$Betas <- res3$Betas[, c(res3$Mu) > 0.01]
Mat_est <- estimateTotalCounts.SigPoisProcess_mult(object = res3,
                                        bin_weight = bin_weigth,
                                        SignalTrack = SignalTrack_red)
MutMatrix <- SigPoisProcess::getTotalMutations(gr_Mutations_red)
plot(Mat_est, MutMatrix)

gr_Mutations$Encode01
X <- as.matrix(mcols(gr_Mutations)[, -c(1,2)])
bigProd <- exp(X %*% res$Betas) * t(res$Thetas)[gr_Mutations$sample, ] * res$Signatures[gr_Mutations$sample, ]
Probs <- bigProd/rowSums(bigProd)
mm <- apply(Probs, 1, which.max)
table(mm)

res$Mu
hist(Probs[, 9], breaks = 100)
round(quantile(Probs[, 8], 0.99), 4)


out <- CompressiveNMF::CompressiveNMF_map(MutMatrix)
plot_96SBSsig(out$Signatures)

R_start <- out$Signatures
Theta_start <- out$Theta/sum(bin_weigth)


Means <- colMeans(SignalTrack_red)
sds <- apply(SignalTrack_red, 2, sd)
SignalTrack_std <- SignalTrack_red - t(Means)[rep(1, nrow(SignalTrack_red)),]
SignalTrack_std <- SignalTrack_std / t(sds)[rep(1, nrow(SignalTrack_red)),]

Xnorm <- as.matrix(mcols(gr_Mutations)[, 3:12])
Xnorm <- Xnorm - t(Means)[rep(1, nrow(Xnorm)), ]
Xnorm <- Xnorm / t(sds)[rep(1, nrow(Xnorm)), ]


gr_Mutations_red <- gr_Mutations
mcols(gr_Mutations_red) <- cbind(mcols(gr_Mutations)[, 1:2], Xnorm)

set.seed(73)
res <- SigPoisProcess_mult(gr_Mutations = gr_Mutations_red,
                           SignalTrack = SignalTrack_std,
                           bin_weight = bin_weigth,
                           K = 8,
                           method = "map",
                           update_Signatures = FALSE, update_Theta = FALSE,
                           init = SigPoisProcess_mult.init(R_start = data$R, Theta_start = data$Theta),
                           controls = SigPoisProcess_mult.controls(maxiter = 100,
                                                                   tol = 1e-7),
                           prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                          c0 =  5 + 1, d0 = 0.01 * 5))

res$sol$Sigma2
res$sol$Mu
res$sol$Betas

sqrt(sum(res$sol$Betas^2))
round(res$sol$Betas, 5)

plot(res$sol$trace[-c(1:5)])
plot_96SBSsig(res$Signatures)
round()
dim(SignalTrack)
mutMatrix

match <- match_MutSign(R_true =  data$R, R_hat = res$Signatures)
plot(res$sol$Betas[, match$match][1:5, ], data$Betas)
abline(a = 0, b = 1)

match <- match_MutSign(R_true =  data$R, R_hat = R_start)
plot(R_start[, match$match], data$R)
abline(a = 0, b = 1)

res$Mu

plot(res$sol$Theta[match$match, ], data$Theta)
abline(a = 0, b = 1)

match <- match_MutSign(R_true = data$R, R_hat = out$Signatures)
plot(Theta_start[match$match, ], data$Theta)
abline(a = 0, b = 1)

plot(R_start[, match$match], data$R)
abline(a = 0, b = 1)

plot(MutMatrix, R_start %*% Theta_start)

# How well do we retrieve the mutation rate?
res$Signatures
res$Thetas
res$Betas

Xest <- estimateTotalCounts.SigPoisProcess_mult(res, SignalTrack, bin_weight)
plot(Xest, MutMatrix, xlim = c(0, 100), ylim = c(0, 100))

i <- sample(1:96, size = 1)
j <- sample(1:30, size = 1)
Track_est <- estimateMutationRate.SigPoisProcess_mult(i, j, R = res$Signatures,
                                                      Theta =res$Thetas, Betas = res$Betas,
                                                      SignalTrack, bin_weigth)

Track_true <- estimateMutationRate.SigPoisProcess_mult(i, j, R = data$R,
                                                      Theta = data$Theta,
                                                      Betas = data$Betas,
                                                      SignalTrack[, 1:5], bin_weigth)
plot(Track_true, type = "l")
lines(Track_est, col = "red")




#Xcov <- -log10(pnorm(log(Xcov)))
CumSums <- apply(exp(Xcov %*% pars$Betas),
                 2, function(x) cumsum(width(grX) * x))
i <- sample(1:96, size = 1)
j <- sample(1:50, size = 1)

mm <- exp(Xcov %*% pars$Betas)
Lambda_Tmax <- c(crossprod(pars$R[i, ], Theta_adj[, j] * c(tail(CumSums, 1))))
Lambdas <- c(crossprod(pars$R[i, ], Theta_adj[, j] * t(CumSums)))
plot(Lambdas, type = "l")
Lambda_Tmax


#
# LambdaMat <- matrix(NA, nrow = 96, ncol = 50)
# for(i in 1:96){
#  for(j in 1:50){
#    LambdaMat[i, j] <- c(crossprod(pars$R[i, ], Theta_adj[, j] * c(tail(CumSums, 1))))
#  }
# }
# hist(LambdaMat)





X <- as.matrix(mcols(data$gr_Mutations)[, -c(1,2)])
SignalTrack <- as.matrix(mcols(data$gr_SignalTrack))
bin_weight <- width(data$gr_SignalTrack)

# Starting solution
MutMatrix <- getTotalMutations(data$gr_Mutations)
StartSol <- CompressiveNMF::CompressiveNMF_map(MutMatrix)
R_start <- StartSol$Signatures
Theta_start <- StartSol$Theta/sum(width(data$gr_SignalTrack))
Betas_start <- matrix(0, nrow = ncol(X), ncol = ncol(R_start))#matrix(rnorm(K * p, sd = 0.01), ncol = p)

channel_id = as.numeric(data$gr_Mutations$channel) - 1
sample_id = as.numeric(data$gr_Mutations$sample) - 1

set.seed(42)
K <- 10
R_start <- sample_signatures(matrix(1, 96, K))
Theta_start <- matrix(rgamma(K * 50, 1,1)/200000, nrow = K)
Betas_start <- matrix(0, nrow = ncol(X), ncol = ncol(R_start))
#Betas_start <- matrix(rnorm(p * K, sd = 2), nrow = p)

SigPrior <- matrix(1.1, nrow = 96, ncol = K)

sol <- PPF_loglinear(R_start = R_start,
                     Theta_start = Theta_start,
                     Betas_start = Betas_start,
                     X = X,
                     SignalTrack = SignalTrack,
                     bin_weight = bin_weight,
                     channel_id = channel_id,
                     sample_id = sample_id,
                     SigPrior = SigPrior,
                     a = 1.1, a0 = 1.1 * 50 + 1, b0 = 1.1 * 50 * 0.01,
                     update_R = TRUE,
                     update_Theta = TRUE,
                     update_Betas = TRUE,
                     method = "map",
                     maxiter = 1000,
                     n_iter_betas = 2, tol = 1e-13)

sol$Mu
library(sigminer)
match_id <- match_MutSign(data$R, sol$R)$match
-1911822.90429
-1976667.46997
rownames(sol$R) <- rownames(data$R)
#colnames(sol$R)[match_id] <- paste0("Sig", 1:ncol(sol$R))
colnames(sol$R) <- paste0("Sig", 1:ncol(sol$R))
plot_96SBSsig(sol$R)

plot(sol$R[, match_id[1:8]], data$R); abline(a = 0, b = 1)
plot(sol$Theta[match_id[1:8], ], data$Theta); abline(a = 0, b = 1)
plot(sol$Betas[, match_id[1:8]], t(data$Betas)); abline(a = 0, b = 1)

sol$Betas[, -match_id[1:8]]
sol$Theta[-match_id[1:8], ]

plot(sol$R[, match_id[-c(9,10)]], data$R); abline(a = 0, b = 1)
plot(sol$Theta[match_id[-c(9,10)], ], data$Theta); abline(a = 0, b = 1)
plot(sol$Betas[, match_id[-c(9,10)]], t(data$Betas)); abline(a = 0, b = 1)

sol$Betas[, c(5, 11)]
sol$Betas[, 11] <- 0
sol$R[, 11] <- 1/96
sol$Theta[5, ] <- sol$Theta[5, ] + sol$Theta[11, ]
sol$Theta[11, ] <- 0

plot(sol$Theta[5, ], data$Theta[1, ])

data$Betas[1, ]
plot(sol$Theta[5, ], data$Theta[1, ]); abline(a = 0, b = 1)

tots <- colSums(bin_weight * exp(as.matrix(SignalTrack %*% sol$Betas)))

plot(sol$R %*% (sol$Theta * matrix(tots)[, rep(1, 50)]), MutMatrix)

for(i in 1:96){
  for(j in 1:50)
}
tots * sol$Theta[, 1]

CumSums <- apply(exp(as.matrix(SignalTrack %*% sol$Betas)),
                 2, function(x) cumsum(width(grX) * x))
Lambda_Tmax <- c(crossprod(R[i, ], Theta[, j] * c(tail(CumSums, 1))))


set.seed(10)
Betas <- matrix(rnorm(p * 8, sd = 1), nrow = 10)#t(matrix(0, nrow = nrow(data$Theta), ncol = p))
Betas <- Update_Beta(X = X,
                     SignaTrack = SignalTrack,
                     bin_weight = bin_weight,
                     R = data$R,
                     Theta = data$Theta,
                     Betas = Betas,
                     channel_id = channel_id,
                     sample_id = sample_id,
                     maxiter = 10)

plot(Betas, t(data$Betas)); abline(a = 0, b = 1)
solve(Hess)
rowSums(Theta_start)

R2 <- R
tmp <- Update_R(X, R2, Theta, Betas, channel_id, sample_id)
R2 <- tmp
plot(tmp, R)


TT <- Update_Theta(X = X,
             SignaTrack = SignalTrack,
             bin_weight = bin_weight,
             R = data$R,
             Theta = TT,#data$Theta,
             Betas = Betas,
             channel_id = channel_id,
             sample_id = sample_id)

Betas_start <- matrix(0, nrow = p, ncol = 8)
res <- PPF_loglinear(R_start = R_start,
                     Theta_start = Theta_start,
                     Betas_start = Betas_start,
                     X = X,
                     SignalTrack = SignalTrack,
                     bin_weight = bin_weight,
                     channel_id = channel_id,
                     sample_id = sample_id,
                     update_R = TRUE,
                     update_Theta = TRUE,
                     update_Betas = TRUE)

plot(res$Theta, data$Theta)
plot(res$R, R_start)
plot(res$Betas, data$Betas)

gr_Mutations <- data$gr_Mutations

MutMatrix <- getTotalMutations(data$gr_Mutations)
StartSol <- CompressiveNMF::CompressiveNMF_map(MutMatrix)
#plot_96SBSsig(StartSol$Signatures)

J <- 50
Xsignal <- cbind(as.matrix(mcols(data$gr_SignalTrack)))
scoreTracks <- array(NA,  dim = c(nrow(Xsignal), J, ncol(Xsignal)),
                dimnames = list(1:nrow(Xsignal),
                                unique(data$gr_Mutations$sample),
                                colnames(Xsignal)))
for(j in 1:J) scoreTracks[, j, ] <- Xsignal
Bins <- matrix(width(data$gr_SignalTrack))[, rep(1, J)]

p <- ncol(Xsignal)
X <- as.matrix(mcols(data$gr_Mutations)[, -c(1,2)])
channel_id = as.numeric(data$gr_Mutations$channel) - 1
sample_id = as.numeric(data$gr_Mutations$sample) - 1

R_start <- StartSol$Signatures
Theta_start <- StartSol$Theta/sum(width(data$gr_SignalTrack))
Betas_start <- matrix(0, nrow = nrow(Theta_start), ncol = p)#matrix(rnorm(K * p, sd = 0.01), ncol = p)
colnames(Betas_start) <- colnames(Xsignal)

tp <- matrix(rgamma(J * K, (colSums(MutMatrix))), nrow = K, byrow = TRUE)
plot(colSums(MutMatrix), colSums(tp))

plot(MutMatrix, R_start %*% tp)

nnls::nnls(b = MutMatrix[,2], A = R_start)$x

set.seed(10)

R_start <- t(LaplacesDemon::rdirichlet(n = 8, alpha = rep(1, 96)))
Theta_start <- matrix(rgamma(8 * 50, 10, 0.5), nrow = 8)/sum(width(data$gr_SignalTrack))
sol <- compute_MAP_loglinear(R_start = R_start,
                             Theta_start = Theta_start,
                             Betas_start = Betas_start,
                             X = X,
                             channel_id = channel_id,
                             sample_id = sample_id,
                             Xtotal = Xtotal,
                             Bins = Bins,
                             maxiter = 50)

rownames(sol$R) <- rownames(data$R)
colnames(sol$R) <- paste0("Sig", 1:ncol(sol$R))
plot_96SBSsig(sol$R)
match_id <- match_MutSign(data$R, sol$R)$match

abline(a=0, b = 1)
plot(sol$R, R_start)
plot(sol$Theta, Theta_start)
rowname(sol$R) <- dimnames(R_start)
plot_96SBSsig(sol$R)
plot_96SBSsig(R_start)


set.seed(10)
R_start <- t(LaplacesDemon::rdirichlet(n = 12, alpha = rep(1, 96)))
Theta_start <- matrix(rgamma(12 * 50, 10, 0.5), nrow = 12)/sum(width(data$gr_SignalTrack))
Betas_start <- matrix(0, nrow = 12, ncol = 10)
sol2 <- compute_MAP_loglinear(R_start = R_start,
                             Theta_start = Theta_start,
                             Betas_start = Betas_start,
                             X = X,
                             channel_id = channel_id,
                             sample_id = sample_id,
                             Xtotal = Xtotal,
                             Bins = Bins,
                             maxiter = 20)

rownames(sol2$R) <- rownames(data$R)
colnames(sol2$R) <- paste0("Sig", 1:ncol(sol2$R))
plot_96SBSsig(sol2$R)


match_id <- match_MutSign(data$R, sol2$R)$match

abline(a=0, b = 1)
plot(sol$R, R_start)
plot(sol$Theta, Theta_start)

plot(sol2$Betas[12, ], data$Betas[1, ])


rowname(sol$R) <- dimnames(R_start)
plot_96SBSsig(sol$R)
plot_96SBSsig(R_start)


library(sigminer)

match_id <- match_MutSign(data$R, sol2$R)$match
plot(sol$Betas[match_id, ], data$Betas); abline(a = 0, b = 1)
plot(sol$Theta[match_id, ], data$Theta); abline(a = 0, b = 1)
plot(sol$Theta[match_id, ], data$Theta); abline(a = 0, b = 1)

solve(sol$Hess) - solve(.5 *sol$Hess + .5 *t(sol$Hess))
sol$Hess
solve(sol$Hess, sol$grad2 - sol$grad1)


#####
plot_96SBSsig(signatures = Sigs)
library(sigminer)
match_id <- match_MutSign(data$R, Sigs)$match
Betas_true <- rbind(data$Betas, matrix(0, nrow = 2, ncol = p))
Thetas_true <- rbind(data$Theta, matrix(0, nrow = 2, ncol = J))
plot(Betas[match_id, ], Betas_true); abline(a = 0, b = 1)
plot(Thetas[match_id, ], Thetas_true); abline(a = 0, b = 1)


tmp <- 0#sample_Categorical(t(Sigs))
for(j in 1:100000){
  tmp <- tmp +  sample_Categorical(t(Sigs))
}
plot(tmp/100000, t(Sigs))

plot_96SBSsig(sol$R)
sol$Betas
sol$Theta

R <- R_start
Theta <- Theta_start
Betas <- Betas_start
x <- Update_Theta_Betas_loglinear(sol$R,
                                  sol$Theta,
                                  Betas = sol$Betas,
                                  X = X,
                                  channel_id = channel_id,
                                  sample_id = sample_id,
                                  Xtotal,
                                  Bins)
x$Hess
compute_expXtB_field(X, 96, J, channel_id, sample_id, sol$Betas)[2, 1]


rr <- Update_Signatures_MLE_loglinear2(R_start, Theta_start, Betas_start, X, channel_id, sample_id)
tt <- Update_Theta_MLE_loglinear(sol$R, Theta_start, Betas_start, X, channel_id,
                           sample_id, Xtotal, Bins)

Betas_start <- matrix(0, nrow = nrow(Theta_start), ncol = p)
Betas <- Update_Theta_Betas_loglinear(data$R, data$Theta, Betas_start, X,
                             channel_id, sample_id, Xtotal, Bins, iter = 50)

Betas <- Update_Theta_Betas_loglinear(data$R, data$Theta, Betas_start, X,
                                      channel_id, sample_id, Xtotal, Bins, iter = 30)
plot(Betas, data$Betas)
plot(Betas[3, ], data$Betas[3, ])
abline(a = 0, b = 1)

Xfield <- build_X_field2(X, 96, J, channel_id, sample_id)
Betas <- data$Betas#matrix(0, nrow = nrow(Theta_start), ncol = p)
for(i in 1:5){
  print(i)
  plot(Betas, data$Betas); abline(a = 0, b = 1)
  Betas <- update_Betas(data$Theta, data$R, Betas, Xfield, Xtotal,
                        bins = Bins[1, ], MutMatrix = MutMatrix)
}



plot(sol2$Ej, tt)
plot(sol$Theta, Theta_start)
plot(sol$Theta, Theta_start)
plot(sol$R, R_start)
plot(sol$Betas, data$Betas)

plot(data$Theta, Theta_start)

# log10p <- import("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bigWig/ENCFF356KVD.bigWig")
# fc <- import("~/ENCFF136OER.bigWig")
# log10p2 <- log10p[1:1e6]
# fc2 <- fc[1:1e6]
# plot(log2(fc2$score))
# overlaps <- findOverlaps(fc2, log10p2)
# sc1 <- log10p2$score[subjectHits(overlaps)]
# sc2 <- fc2$score[queryHits(overlaps)]
# plot(log2(sc2 + 1e-12), type = "l")
# quantile(sc2)
#
# plot(sc2, sc1)
# points(sc2, sc1)
#
# hist(sc2[sc2>0], breaks = 40)
#
# plot(log2(sc2 + 1e-12), type = "l")
#
# plot(rep(fc2$score, width(fc2)), type = "l")
#
# overlapPairs <- findOverlapPairs(fc2, log10p2)
# overlapPairs$log10p <- log10p2$score[subjectHits(overlaps)]
# overlapPairs$fc <- fc2$score[queryHits(fc2)]
#
#
#
# methyl <- import("~/SigPoisProcess/data/ENCODE_pvalues/Methylation_data/Geo_methylation_bigWig/GSM5652369_Gastric-antrum-Endocrine-Z00000438.bigwig")
# overlaps <- findOverlaps(fc, log10p)
# overlapPairs <- findOverlapPairs(fc, log10p)
# overlapPairs$log10p <- log10p$score[subjectHits(log10p)]
# overlapPairs$fc <- fc$score[queryHits(fc)]
#
# score <- fc$score[1:1000]
# score2 <- log10p$score[1:1000]
#
# plot(score, score2)
#
# sum(width(log10p))
# sum(width(fc))
#
# hist(fc$score)
# findOver
#
# length(fc)
#
# tmp <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Breast-Cancer_female_H3K9me3.bigWig")
# hist(tmp$score, breaks=50)
# hist(scale(qnorm(1- 10^-tmp$score)), breaks=50)
#
# pvals <- pmax(10^-tmp$score, 0.9999999)
# pp <- 1 - qnorm(pvals + 1e-18)
# quantile(pp, 0.999999)
#
#
#
#
df_tmp <- as.data.frame(data$gr_Mutations)
over <- findOverlaps(data$gr_SignalTrack, data$gr_Mutations)
df_tmp$region <- queryHits(over)
df_all <- df_tmp %>%
  group_by(region)%>%
  summarize(n = n())
plot(df_all$region, df_all$n)



sol <- function(a){
  sqrt(sqrt(a^2 + 10 * a + 1) + a - 1)/sqrt(6)
}

sol(0.00000)

