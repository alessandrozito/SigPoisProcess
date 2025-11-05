library(tidyverse)
library(maftools)
library(BSgenome)
library(Biostrings)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(SigPoisProcess)

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
  X <- -log10(pnorm(X)) + 1e-12
  colnames(X) <- paste0("Encode", sprintf("%02d", 1:ncol(X)))
  # Append and return
  mcols(grX) <- X
  return(grX)
}


generate_Parameters <- function(cosmic_sigs = c("SBS1", "SBS2", "SBS13", "SBS3"),
                                K_new = 4,
                                J = 50, a = 1, p = 5,
                                theta = 100,
                                ngaps = 10,
                                alpha = 0.25){

  Rmat_cos <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37[, cosmic_sigs]
  if(K_new > 0){
    Rmat_random <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, 96)))
    colnames(Rmat_random) <- paste0("SBSnew", 1:K_new)
  }
  Rmat <- cbind(Rmat_cos, Rmat_random)

  # Generate the weights
  K_true <- ncol(Rmat)
  Mu_true <- matrix(rgamma(ncol(Rmat) * p, theta, 1), nrow = K_true)
  # Fill it with random gaps
  if(ngaps > 0){
    id <- sample(1:length(Mu_true), ngaps)
    Mu_true[id] <- 1e-5
  }
  # Simulate Betas
  Betas_true <- array(NA, dim = c(K_true, J, p))
  for(j in 1:J) Betas_true[, j, ] <- matrix(rgamma(K_true * p, a, a/Mu_true),
                                            nrow = K_true)
  return(list(R = Rmat, Betas = Betas_true, Mu = Mu_true))
}


sample_mutations <- function(i, j, R, Betas, RangesX, CumSums){
  Lambda_Tmax <- c(crossprod(R[i, ], Betas[, j, ] %*% c(tail(CumSums, 1))))
  Lambdas <- c(crossprod(R[i, ], Betas[, j, ] %*% t(CumSums)))
  N_Tmax <- rpois(1, Lambda_Tmax)
  if(N_Tmax > 0){
    # Sample from the uniform
    u <- sort(runif(n = N_Tmax))
    # Invert to find location
    ts <- sapply(u, function(uu) {
      ids <- which(Lambdas <= uu * Lambda_Tmax)
      if(length(ids)==0){
        id <- 1
      } else {
        id <- max(ids)
      }
      rg <- RangesX[id]#ranges(grX)[id]
      t <- sample(start(rg):end(rg), size = 1)
      t
    })
  } else {
    ts <- NA
  }
  return(ts)
}


sample_dataset <- function(R, Betas, RangesX, CumSums){
  # Model dimensions
  I <- nrow(R)
  J <- ncol(Betas)
  Mutations <- data.frame()
  for(j in 1:J) {
    print(paste0("Simulating Patient_", sprintf("%02d", j)))
    for(i in 1:I) {
      ch <- rownames(R)[i]
      ts <- sample_mutations(i, j, R, Betas, RangesX, CumSums)
      if(!is.na(ts[1])) {
        Mutations <- rbind(Mutations, data.frame(sample = paste0("Patient_", j),
                                                 pos = ts,
                                                 channel = ch))
      }
    }
  }
  Mutations$sample <- factor(Mutations$sample, levels = unique(Mutations$sample))
  Mutations$channel <- as.factor(Mutations$channel)
  return(Mutations)
}


generate_MutationData <- function(J = 50,
                                  cosmic_sigs = c("SBS1", "SBS2", "SBS13", "SBS3"),
                                  K_new = 4,
                                  theta = 100,
                                  ngaps = 10,
                                  a = 1,
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
  CumSums <- apply(mcols(grX)[, 1:p_to_use], 2, function(x) cumsum(width(grX) * x))

  # Step 3 - Sample the model parameters
  pars_true <- generate_Parameters(cosmic_sigs = cosmic_sigs,
                                   K_new = K_new, J = J,
                                   a = a, p = p_to_use, theta = theta,
                                   ngaps = ngaps)
  # Step 4 - Rescale the Betas
  X_totals <- CumSums[rep(nrow(CumSums), J), ]
  X_bar <- array(NA, dim = c(ncol(pars_true$R), J, p_to_use))
  for(k in 1:nrow(X_bar)) X_bar[k, ,] <- X_totals
  Betas_scaled <- pars_true$Betas / X_bar

  # Step 5 - Generate the dataset of mutations
  Mutations <- sample_dataset(R = pars_true$R, Betas = Betas_scaled, RangesX = ranges(grX), CumSums = CumSums)

  # Step 6 - Normalize the covariates to have mean 1, and then merge with covariates
  gr_Mutations <- GRanges(seqnames = "chrsim",
                   IRanges(start = Mutations$pos, end = Mutations$pos),
                   sample = Mutations$sample,
                   channel = Mutations$channel)
  grX_norm <- grX
  mcols(grX_norm) <- apply(mcols(grX_norm), 2, function(x) x/mean(x))
  overlaps <- findOverlaps(gr_Mutations, grX_norm)
  mcols(gr_Mutations) <- cbind(mcols(gr_Mutations), mcols(grX_norm[subjectHits(overlaps)]))

  # Step 7 - Find the areas
  totals <- t(crossprod(as.matrix(mcols(grX)), width(grX)))
  df_areas = data.frame("sample" = levels(gr_Mutations$sample),
                           totals[rep(1, J), ])
  totals_norm <- t(crossprod(as.matrix(mcols(grX_norm)), width(grX_norm)))
  df_areas_norm = data.frame("sample" = levels(gr_Mutations$sample),
                             totals_norm[rep(1, J), ])
  # Step 8 - Save all results
  data <- list("gr_Mutations" = gr_Mutations,
               "df_areas" = df_areas,
               "df_areas_norm" = df_areas_norm,
               "covariate_used" = colnames(mcols(grX))[1:p_to_use],
               "R" = pars_true$R,
               "Betas" = Betas_scaled,
               "Mu" = pars_true$Mu,
               "gr_SignalTrack" = grX,
               "simulation_parameters" = list(J = J,
                                              cosmic_sigs = cosmic_sigs,
                                              K_new = K_new,
                                              theta = theta,
                                              ngaps = ngaps,
                                              a = a,
                                              length_genome = length_genome,
                                              tilewidth = tilewidth,
                                              rho = rho,
                                              p_all = p_all,
                                              p_to_use = p_to_use,
                                              ngaps = ngaps,
                                              corr = corr))
  return(data)

}






