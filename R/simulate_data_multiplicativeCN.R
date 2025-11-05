################################################################################
# This file runs the simulation in Section 5 to reproduce the following:
# --- Table 2: Perfomance comparison of the models under two scenarios
# --- Figure 2: Track on the signatures
# --- Table 1.2: xxx
################################################################################

#--- Load the librraries
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
library(knitr)

# Load the R package sigPoisProcess.
devtools::document()

open_rds_file <- function(file){
  if(file.exists(file)){
    out <- readRDS(file)
  } else {
    out <- NULL
  }
  return(out)
}

create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}


################################################################################
# Part 1 - Functions to simulate the data.
################################################################################

generate_SignalTrack <- function(length_genome = 2e6,
                                 tilewidth = 100,
                                 rho = 0.99, p = 10,
                                 Sigma = diag(p)){
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
  colnames(X) <- paste0("Encode", sprintf("%02d", 1:ncol(X)))
  # Append and return
  grX$bin_weight <- tilewidth
  mcols(grX) <- cbind(mcols(grX), scale(X))
  return(grX)
}

generate_CopyTrack <- function(J = 50,
                               mu_copy = 1,
                               size_copy = 1,
                               length_genome = 2e6,
                               tilewidth = 100,
                               mean_segments = 20){
  grX <- tileGenome(c("chrsim" = length_genome), tilewidth = tilewidth,
                    cut.last.tile.in.chrom = TRUE)
  n_tiles <- length(grX)
  cnmat <- matrix(NA, nrow = n_tiles, ncol = J)
  for (j in seq_len(J)) {
    # number of CNA segments
    n_seg <- rpois(1, lambda = mean_segments)
    # segment boundaries
    segs <- sort(sample(1:n_tiles, n_seg, replace = FALSE))
    segs <- c(1, segs, n_tiles+1)
    profile <- numeric(n_tiles)
    for (k in seq_len(length(segs)-1)) {
      s <- segs[k]
      e <- segs[k+1] - 1
      if (e < s) next
      # CN state for this segment
      cn <- 2 + rnbinom(1, mu = mu_copy, size = size_copy)
      profile[s:e] <- cn
    }
    cnmat[, j] <- profile
  }
  colnames(cnmat) <- paste0("Patient_", sprintf("%02d", 1:J))
  mcols(grX) <- as.data.frame(cnmat)
  return(grX)
}

generate_Parameters <- function(cosmic_sigs = c("SBS1", "SBS2", "SBS13", "SBS3"),
                                K_new = 4,
                                J = 50,
                                a = 1,
                                p_all = 10,
                                p_to_use = 5,
                                theta = 200,
                                prob_zero = 0,
                                sd_beta = 1/2,
                                alpha = 0.1){
  p <- p_to_use
  Rmat_cos <- SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37[, cosmic_sigs]
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
  # Pad Betas with zeros to simulate redundant covariates
  if(p_all > p_to_use) {
    Betas <- rbind(Betas, matrix(0, nrow = p_all - p_to_use, ncol = K_true))
  }

  return(list(R = Rmat, Theta = Theta, Betas = Betas))
}


sample_mutations <- function(i, j, R, Theta, RangesX, ExpBetas, CopyTrack,
                             tilewidth = 100){
  # CopyTrack
  CumSums <- apply(ExpBetas, 2, function(x) cumsum(0.5 * tilewidth * CopyTrack[, j] * x))
  # Total expected number of
  Lambda_Tmax <- c(crossprod(R[i, ], Theta[, j] * c(tail(CumSums, 1))))
  Lambdas <- c(crossprod(R[i, ], Theta[, j] * t(CumSums)))
  # Sample the points
  N_Tmax <- rpois(1, Lambda_Tmax)
  # Distribute the points
  if(N_Tmax > 0){
    # Sample from the uniform
    u <- sort(runif(n = N_Tmax))
    ts <- findInterval(u, Lambdas/Lambda_Tmax) * tilewidth + sample(1:(tilewidth - 1), size = length(u), replace = TRUE)
  } else {
    ts <- NA
  }
  return(ts)
}

sample_dataset_parallelCN <- function(R, Theta, Betas, RangesX, ExpBetas,
                                      CopyTrack, tilewidth, ncores){
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
      ts <- sample_mutations(i, j, R, Theta, RangesX, ExpBetas, CopyTrack, tilewidth)
      if(!is.na(ts[1])) {
        mut_temp <- rbind(mut_temp, data.frame(sample = paste0("Patient_", sprintf("%02d", j)),
                                               pos = ts,
                                               channel = ch))
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
                                  theta = 100,
                                  ngaps = 10,
                                  a = 1,
                                  mu_copy = 1,
                                  size_copy = 10,
                                  sd_beta = 1/2,
                                  length_genome = 2e6,
                                  tilewidth = 100,
                                  rho = 0.99,
                                  p_all = 10,
                                  p_to_use = 5,
                                  corr = "indep",
                                  ncores = 1){

  # Step 1 - generate a covariance matrix for the
  if (corr == "indep") {
    Sigma <- diag(p_all)
  } else if (corr == "onion"){
    corrmat <- clusterGeneration::genPositiveDefMat(p_all, covMethod="onion")$Sigma
    Sigma <- cov2cor(corrmat)
  }

  # Step 2 - Generate the SignalTrack
  grX <- generate_SignalTrack(length_genome = length_genome,
                              tilewidth = tilewidth, rho = rho, p = p_all,
                              Sigma = Sigma)

  # Step 3 - Generate the CopyTrack
  gr_CopyTrack <- generate_CopyTrack(J = J, mu_copy = mu_copy, size_copy = size_copy,
                                     length_genome = length_genome,
                                     tilewidth = tilewidth)

  # Step 3 - Sample the model parameters
  pars_true <- generate_Parameters(cosmic_sigs = cosmic_sigs,
                                   K_new = K_new, J = J,
                                   a = a,
                                   p_all = p_all,
                                   p_to_use = p_to_use,
                                   theta = theta,
                                   sd_beta = sd_beta)

  # Step 4 - Rescale the tracks and the hyperparameters
  Xcovs <- mcols(grX)[, -1] # Remove bin_weights
  bin_weight = grX$bin_weight
  ExpBetas <- exp(as.matrix(Xcovs) %*% pars_true$Betas)
  CopyTrack <- as.matrix(mcols(gr_CopyTrack))

  #CopyTots <- colSums(apply(CopyTrack/2, 2, function(x) x * bin_weight))
  colnames(pars_true$Theta) <- names(CopyTrack)
  Theta_scaled <- pars_true$Theta / sum(bin_weight)

  # Step 5 - Generate the dataset of mutations
  Mutations <- sample_dataset_parallelCN(R = pars_true$R,
                                          Betas = pars_true$Betas,
                                          Theta = Theta_scaled,
                                          RangesX = ranges(grX),
                                          ExpBetas = ExpBetas,
                                          CopyTrack = CopyTrack,
                                         tilewidth = tilewidth)

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
               "gr_CopyTrack" = gr_CopyTrack,
               "simulation_parameters" = list(J = J,
                                              cosmic_sigs = cosmic_sigs,
                                              K_new = K_new,
                                              theta = theta,
                                              a = a,
                                              mu_copy = mu_copy,
                                              size_copy = size_copy,
                                              length_genome = length_genome,
                                              tilewidth = tilewidth,
                                              rho = rho,
                                              p_all = p_all,
                                              p_to_use = p_to_use,
                                              corr = corr))
  return(data)

}

aggregate_SignalTrack_CopyTrack <- function(data, lenght_bins = 200) {

  df <- as.data.frame(data$gr_SignalTrack)
  df_copy <- as.data.frame(data$gr_CopyTrack)

  #------- Aggregate the SignalTrack and CopyTrack
  df$bin_start <- df_copy$bin_start <- floor((df$start - 1) / lenght_bins) * lenght_bins + 1
  df$bin_end <- df_copy$bin_end <- df$bin_start + lenght_bins - 1

  # Summarize SignalTrack
  df_agg_sum <- df %>%
    group_by(seqnames, bin_start, bin_end) %>%
    summarize(across(where(is.numeric), mean))
  gr_agg <- GRanges(seqnames = df_agg_sum$seqnames,
                    ranges = IRanges(start=df_agg_sum$bin_start, end=df_agg_sum$bin_end))
  mcols(gr_agg) <- df_agg_sum[, -(1:6)]
  gr_agg$bin_weight <- lenght_bins

  # Summarize CopyTrack
  df_copy_sum <- df_copy %>%
    group_by(seqnames, bin_start, bin_end) %>%
    summarize(across(where(is.numeric), mean))
  gr_copy <- GRanges(seqnames = df_copy_sum$seqnames,
                     ranges = IRanges(start=df_copy_sum$bin_start, end=df_copy_sum$bin_end))
  mcols(gr_copy) <- df_copy_sum[, -(1:6)]

  # Aggregate SignalTrack using mutate
  df_agg_tot <- df %>%
    group_by(seqnames, bin_start, bin_end) %>%
    mutate(across(where(is.numeric), mean))
  gr_agg_tot <- GRanges(seqnames = df_agg_tot$seqnames,
                        ranges = IRanges(start=df_agg_tot$bin_start, end=df_agg_tot$bin_end))
  mcols(gr_agg_tot) <- df_agg_tot[, -(1:6)] %>%
    ungroup() %>%
    dplyr::select(-bin_start, -bin_end)

  # Aggregate CopyTrack using mutate
  df_copy_tot <- df_copy %>%
    group_by(seqnames, bin_start, bin_end) %>%
    mutate(across(where(is.numeric), mean))
  gr_copy_tot <- GRanges(seqnames = df_copy_tot$seqnames,
                         ranges = IRanges(start=df_copy_tot$bin_start, end=df_copy_tot$bin_end))
  mcols(gr_copy_tot) <- df_copy_tot[, -(1:5)] %>%
    ungroup() %>%
    dplyr::select(-bin_start, -bin_end)

  # Find overlaps with mutation positions
  gr_Mutations_agg <- data$gr_Mutations
  overlaps <- findOverlaps(gr_Mutations_agg, gr_agg)
  mcols(gr_Mutations_agg) <- cbind(mcols(gr_Mutations_agg)[, c(1,2)], mcols(gr_agg[subjectHits(overlaps)]))
  gr_Mutations_agg$bin_weight <- NULL

  # Return quantities
  return(list(gr_Mutations = gr_Mutations_agg,
              SignalTrack = as.matrix(mcols(gr_agg))[, -1],
              SignalTrackRed = as.matrix(mcols(gr_agg_tot)),
              CopyTrack = as.matrix(mcols(gr_copy))/2 * lenght_bins,
              CopyTrackRed = as.matrix(mcols(gr_copy_tot))/2 * data$simulation_parameters$tilewidth))
}

################################################################################
# Part 0 - Functions to evaluate performance
################################################################################
Compute_sensitivity_precision <- function(R_hat, R_true, cos_cutoff = 0.9){
  # Sensitivity: proportion of ground truth signatures that were estimated with
  #              sufficiently high similarity
  sig_sens <- sapply(1:ncol(R_true), function(i)
    max(sapply(1:ncol(R_hat), function(j) sigminer::cosine(R_true[, i], R_hat[, j]))))
  # Precision: proportion of estimated signatures that were sufficiently similar to
  #            ground truth signatures
  sig_prec <- sapply(1:ncol(R_hat), function(i)
    max(sapply(1:ncol(R_true), function(j) sigminer::cosine(R_true[, j], R_hat[, i]))))
  return(c("Sensitivity" = mean(sig_sens > cos_cutoff), "Precision" = mean(sig_prec > cos_cutoff)))
}

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

compute_RMSE_Theta <- function(Theta_true, Theta_hat, match){
  # Need to make sure that Theta_hat and Theta_true are the same
  k_true <- nrow(Theta_true)
  k_hat <- nrow(Theta_hat)
  k_tot <- max(c(k_hat, k_true))
  J <- ncol(Theta_true)
  mat0 <- matrix(0, nrow = abs(k_hat - k_true), ncol = J)
  if (k_hat > k_true) {
    rownames(mat0) <- paste0("new_extra", 1:nrow(mat0))
    Theta_true <- rbind(Theta_true, mat0)
  } else if (k_hat < k_true) {
    Theta_hat <- rbind(Theta_hat, mat0)
  }
  Theta_hat <- Theta_hat[match, ]
  return(sqrt(mean((Theta_true - Theta_hat)^2)))
}

compute_RMSE_Betas <- function(Beta_true, Beta_hat, match) {
  p_true <- nrow(Beta_true)
  k_true <- ncol(Beta_true)
  p_hat <- nrow(Beta_hat)
  k_hat <- ncol(Beta_hat)
  p_max <- max(p_true, p_hat)
  k_max <- max(k_true, k_hat)
  pad_matrix <- function(mat, nrows, ncols) {
    out <- matrix(0, nrows, ncols)
    out[1:nrow(mat), 1:ncol(mat)] <- mat
    out
  }
  Beta_true_pad <- pad_matrix(Beta_true, p_max, k_max)
  Beta_hat_pad  <- pad_matrix(Beta_hat,  p_max, k_max)
  Beta_hat_matched <- Beta_hat_pad[, match, drop = FALSE]
  Beta_true_matched <- Beta_true_pad[, seq_along(match), drop = FALSE]
  stopifnot(all(dim(Beta_true_matched) == dim(Beta_hat_matched)))
  sqrt(mean((Beta_true_matched - Beta_hat_matched)^2))
}

getModelEstimates <- function(res, model = "SigPoisProcess", SignalTrack, CopyTrack){
  # Cutoff values for signatures
  if (model == "CompNMF") {
    cutoff <- 0
  } else {
    cutoff <- 5 * res$PriorParams$a * res$PriorParams$epsilon
  }
  filter_mu <- res$Mu > 5 * cutoff
  filter_sigs <- c(sigminer::cosine(res$Signatures, matrix(rep(1/96, 96))) < 0.975)
  # Cutoff values for signatures
  R_hat <- res$Signatures[, filter_mu & filter_sigs]
  if (model == "CompNMF") {
    theta_baseline <- res$Theta[filter_mu & filter_sigs, ]
    theta_baseline <- theta_baseline/(t(colSums(CopyTrack))[rep(1, nrow(theta_baseline)), ])
    Theta_total <- res$Theta[filter_mu & filter_sigs, ]
    Beta_hat <- matrix(0, nrow = ncol(SignalTrack), ncol = nrow(theta_baseline))
    Lambda_hat <- Reconstruct_Lambda(SignalTrack = SignalTrack, CopyTrack = CopyTrack,
                                     R = R_hat, Theta = theta_baseline, Betas = Beta_hat)
  } else {
    Beta_hat <- res$Betas[, filter_mu & filter_sigs]
    theta_baseline <- res$Thetas[filter_mu & filter_sigs, ]
    Theta_total <- theta_baseline * crossprod(exp(SignalTrack[, rownames(Beta_hat)] %*% Beta_hat), CopyTrack)
    Lambda_hat <- Reconstruct_Lambda(SignalTrack = SignalTrack[, rownames(Beta_hat)], CopyTrack = CopyTrack,
                                     R = R_hat, Theta = theta_baseline, Betas = Beta_hat)
  }
  return(list("R_hat" = R_hat, "Theta_total" = Theta_total,
              "theta_baseline" = theta_baseline, "Beta_hat" = Beta_hat,
              "Lambda_hat" = Lambda_hat))
}

get_counts_from_data <- function(data){
  df_tmp <- as.data.frame(data$gr_Mutations)
  over <- findOverlaps(data$gr_Mutations, data$gr_SignalTrack)
  df_tmp$region <- subjectHits(over)
  df_all <- df_tmp %>%
    group_by(region, sample) %>%
    summarize(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = n, values_fill = 0)  %>%
    complete(region = 1:length(data$gr_SignalTrack), fill = list(n = 0))
  mat_counts <- df_all %>% dplyr::select(-region) %>% as.matrix()
  mat_counts <- mat_counts[, levels(data$gr_Mutations$sample)]
  mat_counts[is.na(mat_counts)] <- 0
  return(mat_counts)
}

get_PosteriorEffectiveSize <- function(chain, burnin = 1500){
  if(is.null(dim(chain))){
    coda::effectiveSize(chain[-c(1:burnin)])
  } else if(length(dim(chain)) == 2) {
    apply(chain[-c(1:burnin), ], 2, function(x) coda::effectiveSize(x))
  } else if (length(dim(chain)) == 3){
    apply(chain[-c(1:burnin), ,], c(2,3), function(x) coda::effectiveSize(x))
  }
}
#
# res <- readRDS("~/SigPoisProcess/output/Simulation/Scenario_A_indep/Simulation_05/output_mcmc_FullModel.rds.gzip")
# res$controls$burnin

compute_RMSE_paramters <- function(res,
                                   data,
                                   Lambda_true,
                                   Theta_total_true,
                                   R_hat,
                                   theta_baseline,
                                   Theta_total,
                                   Beta_hat,
                                   Lambda_hat) {
  Sens_prec <- Compute_sensitivity_precision(R_hat = R_hat, R_true = data$R)
  # Match signatures with the real ones
  MatchedSigs <- match_MutSign(R_true = data$R, R_hat = R_hat)
  # Calculate the RMSE for each quantities
  rmse_sigs <- sqrt(mean((MatchedSigs$R_hat - MatchedSigs$R_true)^2))
  rmse_theta <- compute_RMSE_Theta(Theta_total_true, Theta_total, MatchedSigs$match)
  rmse_Betas <- compute_RMSE_Betas(data$Betas, Beta_hat, MatchedSigs$match)
  rmse_lambda <- sqrt(mean(rowMeans((Lambda_hat - Lambda_true))^2))
  # RMSE for the counts
  counts_true <- get_counts_from_data(data)
  rmse_counts <- sqrt(mean(rowMeans((Lambda_hat - counts_true))^2)) #sqrt(mean((rowSums(Lambda_hat) - counts_true)^2))
  if(!is.null(res$sol)){
    iter <- res$sol$iter
    sampling_details <- c(
      iter              = iter,
      effectiveBetas    = NA,
      effectiveSigs     = NA,
      effectiveTheta    = NA,
      effectiveMu       = NA,
      effectiveSigma2   = NA,
      effectiveLogPost  = NA,
      effectiveLogLik   = NA,
      effectiveLogPrior = NA)
  } else if (!is.null(res$mapOutput)){
    iter <- res$mapOutput$iter
    sampling_details <- c(
      iter              = iter,
      effectiveBetas    = NA,
      effectiveSigs     = NA,
      effectiveTheta    = NA,
      effectiveMu       = NA,
      effectiveSigma2   = NA,
      effectiveLogPost  = NA,
      effectiveLogLik   = NA,
      effectiveLogPrior = NA)
  } else {
    # Calculate average effective sample size of all quantities R, Theta, Beta, Mu, sigma2
    burnin <- res$controls$burnin
    keep <- colnames(R_hat)
    #browser()
    # Effective sample sizes for model parameter
    effectiveBetas <- mean(get_PosteriorEffectiveSize(res$chain$BETASchain[, , keep], burnin = burnin))
    effectiveSigs <- mean(get_PosteriorEffectiveSize(res$chain$SIGSchain[, , keep], burnin = burnin))
    effectiveTheta <- mean(get_PosteriorEffectiveSize(res$chain$THETAchain[, keep, ], burnin = burnin))
    effectiveMu <- mean(get_PosteriorEffectiveSize(res$chain$MUchain[, keep], burnin = burnin))
    effectiveSigma2 <- mean(get_PosteriorEffectiveSize(res$chain$SIGMA2chain[, keep], burnin = burnin))
    # Effective sample sizes for logposterior
    effectiveLogPost <- unname(get_PosteriorEffectiveSize(c(res$chain$logPostchain), burnin = burnin))
    effectiveLogLik <- unname(get_PosteriorEffectiveSize(c(res$chain$logLikchain), burnin = burnin))
    effectiveLogPrior <- unname(get_PosteriorEffectiveSize(c(res$chain$logPriorchain), burnin = burnin))
    iter <- NA
    sampling_details <- c(
      iter              = iter,
      effectiveBetas    = effectiveBetas,
      effectiveSigs     = effectiveSigs,
      effectiveTheta    = effectiveTheta,
      effectiveMu       = effectiveMu,
      effectiveSigma2   = effectiveSigma2,
      effectiveLogPost  = effectiveLogPost,
      effectiveLogLik   = effectiveLogLik,
      effectiveLogPrior = effectiveLogPrior)
  }
  return(c("Kest" = ncol(R_hat),
           Sens_prec,
           "rmse_sig" = rmse_sigs, "rmse_theta" = rmse_theta,
           "rmse_Betas" = rmse_Betas, "rmse_lambda" = rmse_lambda,
           "rmse_counts" = rmse_counts,
           "time" = as.numeric(res$time, units = "mins"),
            sampling_details))
}

postProcessOutput <- function(out_dir){
  # STEP 1 - LOAD THE DATA
  data <- readRDS(paste0(out_dir, "data.rds.gzip"))
  tilewidth <- data$simulation_parameters$tilewidth
  J <- data$simulation_parameters$J

  # SignalTrack and copy number track
  SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
  CopyTrackSim <- tilewidth * as.matrix(mcols(data$gr_CopyTrack)) / 2

  # Calculate the true intensity function
  Lambda_true <- Reconstruct_Lambda(SignalTrack = SignalTrackSim,
                                    CopyTrack = CopyTrackSim,
                                    R = data$R, Theta = data$Theta, Betas = data$Betas)

  Theta_total_true <- data$Theta * crossprod(exp(SignalTrackSim %*% data$Betas), CopyTrackSim)
  # Do not adjust for copy number
  CopyTrackNull <- matrix(tilewidth, nrow = nrow(SignalTrackSim), ncol = J)
  colnames(CopyTrackNull) <- colnames(CopyTrackSim)

  results <- data.frame()
  ###########################
  # Model 1 - CompressiveNMF MAP for the compressive NMF
  ###########################
  outCompNMFBase <- open_rds_file(paste0(out_dir, "output_CompNMFBase.rds.gzip"))
  if(!is.null(outCompNMFBase)){
    print("outCompNMFBase")
    estimates <- getModelEstimates(res = outCompNMFBase,
                                   model = "CompNMF",
                                   SignalTrack = SignalTrackSim,
                                   CopyTrack = CopyTrackNull)
    results_tmp <- compute_RMSE_paramters(res = outCompNMFBase,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results,
                     data.frame("model" = "CompNMFBase",
                                t(results_tmp)))
  }

  ###########################
  # Model 2 - SigPoisProcess with only the copy number
  ###########################
  out_map_CopyOnly <- open_rds_file(paste0(out_dir, "output_map_CopyOnly.rds.gzip"))
  if(!is.null(out_map_CopyOnly)){
    print("out_map_CopyOnly")
    estimates <- getModelEstimates(res = out_map_CopyOnly,
                                   SignalTrack = SignalTrackSim,
                                   CopyTrack = CopyTrackSim)
    results_tmp <- compute_RMSE_paramters(res = out_map_CopyOnly,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results,
                     data.frame("model" = "map_CopyOnly",
                                t(results_tmp)))
  }

  ###########################
  # Model 3 - SigPoisProcess with only the covariates, but no copy number
  ###########################
  out_map_NoCopies  <- open_rds_file(paste0(out_dir, "output_map_NoCopies.rds.gzip"))
  out_mcmc_NoCopies <- open_rds_file(paste0(out_dir, "output_mcmc_NoCopies.rds.gzip"))
  # ------------ 3.1 - Maximum-a-posteriori
  if(!is.null(out_map_NoCopies)){
    print("out_map_NoCopies")
    estimates <- getModelEstimates(res = out_map_NoCopies,
                                   SignalTrack = SignalTrackSim,
                                   CopyTrack = CopyTrackNull)
    results_tmp <- compute_RMSE_paramters(res = out_map_NoCopies,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results, data.frame("model" = "map_NoCopies", t(results_tmp)))
  }
  # ------------ 3.2 - Bayes
  if(!is.null(out_mcmc_NoCopies)){
    print("out_mcmc_NoCopies")
    estimates <- getModelEstimates(res = out_mcmc_NoCopies,
                                   SignalTrack = SignalTrackSim,
                                   CopyTrack =CopyTrackNull)
    results_tmp <- compute_RMSE_paramters(res = out_mcmc_NoCopies,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results, data.frame("model" = "mcmc_NoCopies", t(results_tmp)))
  }

  ###########################
  # Model 4 - SigPoisProcess with only correct covariates
  ###########################
  SignalTrackSim_true <- SignalTrackSim[, 1:data$simulation_parameters$p_to_use]
  gr_Mutations_true <- data$gr_Mutations
  mcols(gr_Mutations_true) <- mcols(gr_Mutations_true)[, c("sample", "channel", colnames(SignalTrackSim_true))]
  out_map_TrueCovs  <- open_rds_file(paste0(out_dir, "output_map_TrueCovs.rds.gzip"))
  out_mcmc_TrueCovs <- open_rds_file(paste0(out_dir, "output_mcmc_TrueCovs.rds.gzip"))
  # ------------ 4.1 - Maximum-a-posteriori
  if(!is.null(out_map_TrueCovs)){
    print("out_map_TrueCovs")
    estimates <- getModelEstimates(res = out_map_TrueCovs,
                                   SignalTrack = SignalTrackSim_true,
                                   CopyTrack =CopyTrackSim)
    results_tmp <- compute_RMSE_paramters(res = out_map_TrueCovs,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results, data.frame("model" = "map_TrueCovs", t(results_tmp)))
  }
  # ------------ 4.2 - Bayes
  if(!is.null(out_mcmc_TrueCovs)){
    print("out_mcmc_TrueCovs")
    estimates <- getModelEstimates(res = out_mcmc_TrueCovs,
                                   SignalTrack = SignalTrackSim_true,
                                   CopyTrack =CopyTrackSim)
    results_tmp <- compute_RMSE_paramters(res = out_mcmc_TrueCovs,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results, data.frame("model" = "mcmc_TrueCovs", t(results_tmp)))
  }

  ###########################
  # Model 5 - SigPoisProcess with also redundant covariates
  ###########################
  out_map_Full <- open_rds_file(paste0(out_dir, "output_map_FullModel.rds.gzip"))
  out_mcmc_Full <- open_rds_file(paste0(out_dir, "output_mcmc_FullModel.rds.gzip"))
  # ------------ 5.1 - Maximum-a-posteriori
  if(!is.null(out_map_Full)){
    print("out_map_Full")
    estimates <- getModelEstimates(res = out_map_Full,
                                   SignalTrack = SignalTrackSim,
                                   CopyTrack = CopyTrackSim)
    results_tmp <- compute_RMSE_paramters(res = out_map_Full,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results, data.frame("model" = "map_Full", t(results_tmp)))
  }
  # ------------ 5.2 - Bayes
  if(!is.null(out_mcmc_Full)){
    print("out_mcmc_Full")
    estimates <- getModelEstimates(res = out_mcmc_Full,
                                   SignalTrack = SignalTrackSim,
                                   CopyTrack = CopyTrackSim)
    results_tmp <- compute_RMSE_paramters(res = out_mcmc_Full,
                                          data = data,
                                          Lambda_true = Lambda_true,
                                          Theta_total_true = Theta_total_true,
                                          R_hat = estimates$R_hat,
                                          theta_baseline = estimates$theta_baseline,
                                          Theta_total = estimates$Theta_total,
                                          Beta_hat = estimates$Beta_hat,
                                          Lambda_hat =  estimates$Lambda_hat)
    results <- rbind(results, data.frame("model" = "mcmc_Full", t(results_tmp)))
  }

  if(nrow(results) >0) {
    #Save Name
    results$Simulation <- basename(out_dir)
    results$Scenario <- basename(dirname(out_dir))
    results$n <- length(data$gr_Mutations)
  }

  return(results)
}

postProcessOutputAgg <- function(out_dir, dims_aggreg = c(200, 500)){
  # STEP 1 - LOAD THE DATA
  data <- readRDS(paste0(out_dir, "data.rds.gzip"))
  tilewidth <- data$simulation_parameters$tilewidth
  J <- data$simulation_parameters$J

  # SignalTrack and copy number track
  SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
  CopyTrackSim <- tilewidth * as.matrix(mcols(data$gr_CopyTrack)) / 2

  # Calculate the true intensity function
  Lambda_true <- Reconstruct_Lambda(SignalTrack = SignalTrackSim,
                                    CopyTrack = CopyTrackSim,
                                    R = data$R, Theta = data$Theta, Betas = data$Betas)

  Theta_total_true <- data$Theta * crossprod(exp(SignalTrackSim %*% data$Betas), CopyTrackSim)

  results <- data.frame()
  for(t in 1:length(dims_aggreg)){
    # Aggregate the signal
    data_agg <- aggregate_SignalTrack_CopyTrack(data, lenght_bins = dims_aggreg[t])
    # Load data
    out_map  <- open_rds_file(paste0(out_dir, "output_map_FullModel_agg", dims_aggreg[t],".rds.gzip"))
    out_mcmc <- open_rds_file(paste0(out_dir, "output_mcmc_FullModel_agg", dims_aggreg[t],".rds.gzip"))
    # ------------ 4.1 - Maximum-a-posteriori
    if(!is.null(out_map)){
      print(paste0("out_map ", dims_aggreg[t]))
      estimates <- getModelEstimates(res = out_map,
                                     SignalTrack = data_agg$SignalTrackRed,
                                     CopyTrack = data_agg$CopyTrackRed)
      results_tmp <- compute_RMSE_paramters(res = out_map,
                                            data = data,
                                            Lambda_true = Lambda_true,
                                            Theta_total_true = Theta_total_true,
                                            R_hat = estimates$R_hat,
                                            theta_baseline = estimates$theta_baseline,
                                            Theta_total = estimates$Theta_total,
                                            Beta_hat = estimates$Beta_hat,
                                            Lambda_hat =  estimates$Lambda_hat)
      results <- rbind(results, data.frame("model" = paste0("map_TrueCovs", dims_aggreg[t]), t(results_tmp)))
    }
    # ------------ 4.2 - Bayes
    if(!is.null(out_mcmc)){
      print(paste0("out_mcmc ", dims_aggreg[t]))
      estimates <- getModelEstimates(res = out_mcmc,
                                     SignalTrack = data_agg$SignalTrackRed,
                                     CopyTrack =data_agg$CopyTrackRed)
      results_tmp <- compute_RMSE_paramters(res = out_mcmc,
                                            data = data,
                                            Lambda_true = Lambda_true,
                                            Theta_total_true = Theta_total_true,
                                            R_hat = estimates$R_hat,
                                            theta_baseline = estimates$theta_baseline,
                                            Theta_total = estimates$Theta_total,
                                            Beta_hat = estimates$Beta_hat,
                                            Lambda_hat =  estimates$Lambda_hat)
      results <- rbind(results, data.frame("model" = paste0("mcmc_TrueCovs", dims_aggreg[t]), t(results_tmp)))
    }
  }
  if(nrow(results) >0) {
    #Save Name
    results$Simulation <- basename(out_dir)
    results$Scenario <- basename(dirname(out_dir))
    results$n <- length(data$gr_Mutations)
  }
  return(results)
}


################################################################################
# Part 2 - Simulate all the datasets in the two scenarios
################################################################################

# Set-up of the simuation
J <- 40             # N. of patients in each dataset
n_datasets <- 20 #40
corr_scenarios <- c("Scenario_A_indep", "Scenario_B_corr")
out_dir <- "~/SigPoisProcess/output/Simulation/"
ncores_simulation <- 20

# Genetate the data
regenerate_data <- FALSE # <---- Set to TRUE if one needs to re-run it
if(regenerate_data) {
  # Set seed for reproducibility
  set.seed(10, kind = "L'Ecuyer-CMRG")
  for(sc in seq_along(corr_scenarios)){
    create_directory(paste0(out_dir, corr_scenarios[sc]))
    corr_type <- if(corr_scenarios[sc] == "Scenario_A_indep") "indep" else "onion"
    for(i in 1:n_datasets) {
      print(paste0("<--- Scenario ", corr_scenarios[sc], " - Simulation ", i, "--->"))
      # Simulate data
      data <- generate_MutationData(J = J, ncores = ncores_simulation,
                                    corr = corr_type)
      # Store the output
      dir_data <- paste0(out_dir, corr_scenarios[sc], "/Simulation_", sprintf("%02d", i), "/")
      create_directory(dir_data)
      file_name <- paste0(dir_data, "data.rds.gzip")
      saveRDS(data, file = file_name, compress = "gzip")
    }
  }
}

################################################################################
# Part 3 - Run all models
################################################################################

# I run the model sequentially using this function
run_models <- function(out_dir,
                       K = 15,
                       a = 1.01, alpha = 1.01, epsilon = 0.001,
                       c0 = 100, d0 = 1,
                       tol = 1e-7,
                       maxiter = 4000,
                       nsamples = 3000,
                       burnin = 1500) {
  # Load the data
  #out_dir <- "~/SigPoisProcess/output/Simulation/Scenario_A_indep/Simulation_01/"
  data <- readRDS(paste0(out_dir, "data.rds.gzip"))
  J <- length(unique(data$gr_Mutations$sample))

  #------- Preliminary useful quantities
  tilewidth <- data$simulation_parameters$tilewidth
  # SignalTrack and copy number track
  SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
  CopyTrackSim <- tilewidth * as.matrix(mcols(data$gr_CopyTrack)) / 2
  # Do not adjust for copy number
  CopyTrackNull <- matrix(tilewidth, nrow = nrow(SignalTrackSim), ncol = J)
  colnames(CopyTrackNull) <- colnames(CopyTrackSim)
  # Aggregate mutation matrix
  MutMatrix <- SigPoisProcess::getTotalMutations(data$gr_Mutations)
  # Controls for SigPoisProcess
  controls <- SigPoisProcess_mult.controls(maxiter = maxiter,
                                           tol = tol,
                                           nsamples = nsamples,
                                           burnin = burnin)
  # Prior Parameters for sigpoisprocess
  prior_params <- SigPoisProcess_mult.PriorParams(c0 = c0, d0 = d0,
                                                  epsilon = epsilon, a = a,
                                                  alpha = alpha)
  Betas_start <- matrix(0, nrow = ncol(SignalTrackSim), ncol = K)

  ###########################
  # Model 1 - CompressiveNMF MAP for the compressive NMF
  ###########################
  start_time <- Sys.time()
  outCompNMFBase <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = K, a = a,
                                                       alpha = alpha, epsilon = epsilon, tol = tol)
  outCompNMFBase$time <- Sys.time() - start_time
  saveRDS(outCompNMFBase, file = paste0(out_dir, "output_CompNMFBase.rds.gzip"), compress = "gzip")

  ###########################
  # Model 2 - SigPoisProcess with only the copy number
  ###########################
  start_time <- Sys.time()
  out_map_CopyOnly <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
                                            SignalTrack = SignalTrackSim,
                                            CopyTrack = CopyTrackSim,
                                            K = K,
                                            method = "map",
                                            update_Signatures = TRUE,
                                            update_Theta = TRUE,
                                            update_Betas = FALSE,
                                            controls = controls,
                                            prior_params = prior_params,
                                            init = SigPoisProcess_mult.init(Betas_start = Betas_start))
  out_map_CopyOnly$time <- Sys.time() - start_time
  saveRDS(out_map_CopyOnly, file = paste0(out_dir, "output_map_CopyOnly.rds.gzip"), compress = "gzip")

  ###########################
  # Model 3 - SigPoisProcess with only the covariates, but no copy number
  ###########################

  # ------------ 3.1 - Maximum-a-posteriori
  start_time <- Sys.time()
  out_map_NoCopies <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
                                            SignalTrack = SignalTrackSim,
                                            CopyTrack = CopyTrackNull,
                                            K = K,
                                            method = "map",
                                            update_Signatures = TRUE,
                                            update_Theta = TRUE,
                                            update_Betas = TRUE,
                                            controls = controls,
                                            prior_params = prior_params)
  out_map_NoCopies$time <- Sys.time() - start_time
  saveRDS(out_map_NoCopies, file = paste0(out_dir, "output_map_NoCopies.rds.gzip"), compress = "gzip")

  # # ------------ 3.2 - Bayes
  # start_time <- Sys.time()
  # out_mcmc_NoCopies <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
  #                                            SignalTrack = SignalTrackSim,
  #                                            CopyTrack = CopyTrackNull,
  #                                            K = K,
  #                                            method = "mcmc",
  #                                            update_Signatures = TRUE,
  #                                            update_Theta = TRUE,
  #                                            update_Betas = TRUE,
  #                                            controls = controls,
  #                                            prior_params = prior_params)
  # out_mcmc_NoCopies$time <- Sys.time() - start_time
  # saveRDS(out_mcmc_NoCopies, file = paste0(out_dir, "output_mcmc_NoCopies.rds.gzip"), compress = "gzip")

  ###########################
  # Model 4 - SigPoisProcess with only correct covariates
  ###########################
  SignalTrackSim_true <- SignalTrackSim[, 1:data$simulation_parameters$p_to_use]
  gr_Mutations_true <- data$gr_Mutations
  mcols(gr_Mutations_true) <- mcols(gr_Mutations_true)[, c("sample", "channel", colnames(SignalTrackSim_true))]
  # ------------ 4.1 - Maximum-a-posteriori
  start_time <- Sys.time()
  out_map_TrueCovs <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_true,
                                            SignalTrack = SignalTrackSim_true,
                                            CopyTrack = CopyTrackSim,
                                            K = K,
                                            method = "map",
                                            update_Signatures = TRUE,
                                            update_Theta = TRUE,
                                            update_Betas = TRUE,
                                            controls = controls,
                                            prior_params = prior_params)
  out_map_TrueCovs$time <- Sys.time() - start_time
  saveRDS(out_map_TrueCovs, file = paste0(out_dir, "output_map_TrueCovs.rds.gzip"), compress = "gzip")

  # # ------------ 4.2 - Bayes
  # init <- SigPoisProcess_mult.init(R_start = out_map_TrueCovs$Signatures,
  #                                  Theta_start = out_map_TrueCovs$Thetas,
  #                                  Betas_start = out_map_TrueCovs$Betas,
  #                                  Mu_start = out_map_TrueCovs$Mu,
  #                                  Sigma2_start = out_map_TrueCovs$Sigma2)
  # start_time <- Sys.time()
  # out_mcmc_TrueCovs <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_true,
  #                                            SignalTrack = SignalTrackSim_true,
  #                                            CopyTrack = CopyTrackSim,
  #                                            K = K,
  #                                            method = "mcmc",
  #                                            update_Signatures = TRUE,
  #                                            update_Theta = TRUE,
  #                                            update_Betas = TRUE,
  #                                            init = init,
  #                                            controls = controls,
  #                                            prior_params = prior_params)
  # out_mcmc_TrueCovs$time <- Sys.time() - start_time
  # saveRDS(out_mcmc_TrueCovs, file = paste0(out_dir, "output_mcmc_TrueCovs.rds.gzip"), compress = "gzip")

  ###########################
  # Model 5 - SigPoisProcess with also redundant covariates
  ###########################

  # ------------ 5.1 - Maximum-a-posteriori
  start_time <- Sys.time()
  out_map_Full <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
                                        SignalTrack = SignalTrackSim,
                                        CopyTrack = CopyTrackSim,
                                        K = K,
                                        method = "map",
                                        update_Signatures = TRUE,
                                        update_Theta = TRUE,
                                        update_Betas = TRUE,
                                        controls = controls,
                                        prior_params = prior_params)
  out_map_Full$time <- Sys.time() - start_time
  saveRDS(out_map_Full, file = paste0(out_dir, "output_map_FullModel.rds.gzip"), compress = "gzip")

  # ------------ 5.2 - Bayes
  init <- SigPoisProcess_mult.init(R_start = out_map_Full$Signatures,
                                   Theta_start = out_map_Full$Thetas,
                                   Betas_start = out_map_Full$Betas,
                                   Mu_start = out_map_Full$Mu,
                                   Sigma2_start = out_map_Full$Sigma2)
  start_time <- Sys.time()
  out_mcmc_Full <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
                                         SignalTrack = SignalTrackSim,
                                         CopyTrack = CopyTrackSim,
                                         K = K,
                                         method = "mcmc",
                                         update_Signatures = TRUE,
                                         update_Theta = TRUE,
                                         update_Betas = TRUE,
                                         init = init,
                                         controls = controls,
                                         prior_params = prior_params)
  out_mcmc_Full$time <- Sys.time() - start_time
  saveRDS(out_mcmc_Full, file = paste0(out_dir, "output_mcmc_FullModel.rds.gzip"), compress = "gzip")

}


run_models_aggreg <- function(out_dir,
                              dims_aggreg = c(200, 500),
                              K = 15,
                              a = 1.01, alpha = 1.01, epsilon = 0.001,
                              c0 = 100, d0 = 1,
                              tol = 1e-7,
                              maxiter = 4000,
                              nsamples = 3000,
                              burnin = 1500) {
  # Load the data
  #out_dir <- "~/SigPoisProcess/output/Simulation/Scenario_A_indep/Simulation_01/"
  data <- readRDS(paste0(out_dir, "data.rds.gzip"))
  J <- length(unique(data$gr_Mutations$sample))
  controls <- SigPoisProcess_mult.controls(maxiter = maxiter,
                                           tol = tol,
                                           nsamples = nsamples,
                                           burnin = burnin)
  # Prior Parameters for sigpoisprocess
  prior_params <- SigPoisProcess_mult.PriorParams(c0 = c0, d0 = d0,
                                                  epsilon = epsilon, a = a,
                                                  alpha = alpha)

  #------- Preliminary useful quantities
  tilewidth <- data$simulation_parameters$tilewidth
  for(t in 1:length(dims_aggreg)){
    print(paste0("Run MAP ", basename(out_dir)," ", basename(dirname(out_dir)), " width ", dims_aggreg[t]))
    # Aggregate the signal
    data_agg <- aggregate_SignalTrack_CopyTrack(data, lenght_bins = dims_aggreg[t])
    # ------------ 5.1 - Maximum-a-posteriori
    start_time <- Sys.time()
    out_map_Full <- SigPoisProcess_multCN(gr_Mutations = data_agg$gr_Mutations,
                                          SignalTrack = data_agg$SignalTrack,
                                          CopyTrack = data_agg$CopyTrack,
                                          K = K,
                                          method = "map",
                                          verbose = FALSE,
                                          controls = controls,
                                          prior_params = prior_params)
    out_map_Full$time <- Sys.time() - start_time
    saveRDS(out_map_Full, file = paste0(out_dir, "output_map_FullModel_agg", dims_aggreg[t],".rds.gzip"), compress = "gzip")
    print(paste0("Run BAYES ", basename(out_dir)," ", basename(dirname(out_dir)), " width ", dims_aggreg[t]))
    # ------------ 5.2 - Bayes
    init <- SigPoisProcess_mult.init(R_start = out_map_Full$Signatures,
                                     Theta_start = out_map_Full$Thetas,
                                     Betas_start = out_map_Full$Betas,
                                     Mu_start = out_map_Full$Mu,
                                     Sigma2_start = out_map_Full$Sigma2)
    start_time <- Sys.time()
    out_mcmc_Full <- SigPoisProcess_multCN(gr_Mutations = data_agg$gr_Mutations,
                                           SignalTrack = data_agg$SignalTrack,
                                           CopyTrack = data_agg$CopyTrack,
                                           K = K,
                                           method = "mcmc",
                                           verbose = FALSE,
                                           init = init,
                                           controls = controls,
                                           prior_params = prior_params)
    out_mcmc_Full$time <- Sys.time() - start_time
    saveRDS(out_mcmc_Full, file = paste0(out_dir, "output_mcmc_FullModel_agg", dims_aggreg[t],".rds.gzip"), compress = "gzip")
  }
}


#--- Run all models, using the standard aggregation
ncores <- 20
dims_aggreg <- c(200, 500)
rerun <- FALSE # <---- Set to TRUE if one needs to re-run it
if(rerun) {
  # Set seed for reproducibility
  set.seed(10, kind = "L'Ecuyer-CMRG")
  for(sc in seq_along(corr_scenarios)){
    corr_type <- if(corr_scenarios[sc] == "Scenario_A_indep") "indep" else "onion"
    # Estimate models in parallel
    registerDoParallel(ncores)
    res <- foreach(i = 1:n_datasets) %dopar% {
      print(paste0("<--- Scenario ", corr_scenarios[sc], " - Simulation ", i, "--->"))
      out_dir <- paste0("~/SigPoisProcess/output/Simulation/",
                        corr_scenarios[sc], "/Simulation_", sprintf("%02d", i), "/")
      # Run the models on the standard aggregation
      run_models(out_dir)
      # Run the models for the total aggregation
      run_models_aggreg(out_dir, dims_aggreg = dims_aggreg)
    }
  }
}


################################################################################
# Part 4 - Post-process the output
################################################################################

# I need to compute:
# 1) Number of signatures estimated
# 2) RMSE for the signatures
# 3) RMSE for Theta-total
# 3) RMSE for Betas
# 4) RMSE for the Lambda across the track.


reprocess_output <- FALSE
if(reprocess_output){
  out_dir_base <- "~/SigPoisProcess/output/Simulation/"
  results_all <- data.frame()
  for(sc in seq_along(corr_scenarios)){
    for(i in 1:n_datasets) {
      print(paste0("<--- Scenario ", corr_scenarios[sc], " - Simulation ", i, "--->"))
      out_dir <- paste0(out_dir_base, corr_scenarios[sc], "/Simulation_", sprintf("%02d", i), "/")
      results_temp <- postProcessOutput(out_dir = out_dir)
      results_temp_agg <- postProcessOutputAgg(out_dir = out_dir, dims_aggreg = dims_aggreg)
      results_all <- rbind(results_all, results_temp, results_temp_agg)
    }
  }
  #write_tsv(results_all, "~/SigPoisProcess/output/Simulation/simulation_results2.tsv")
}

results_all <- read_tsv("~/SigPoisProcess/output/Simulation/simulation_results2.tsv")


#--- Re-load the results
results_all <- results_all %>%
  mutate(F1 = 2 * Precision * Sensitivity / (Precision + Sensitivity))

# Compute and scale means and sds, as you wish (sig*1e3, beta*1e2)
results_summary <- results_all %>%
  mutate(model2 = case_when(
    model == "CompNMFBase" ~ "1. Baseline NMF",
    model == "map_CopyOnly" ~ "2. No covariates",
    model == "map_NoCopies" ~ "3. All covariates, no copy",
    model == "map_TrueCovs" ~ "4. True covariates, true copy",
    model == "map_Full" ~ "5. All covariates, true copy",
    model == "mcmc_Full" ~ "6. MCMC100",
    model == "mcmc_TrueCovs200" ~ "7. MCMC200",
    model == "mcmc_TrueCovs500" ~ "8. MCMC500",
    TRUE ~ NA
  )) %>%
  filter(!is.na(model2)) %>%
  group_by(model2, Scenario) %>%
  summarise(across(where(is.numeric),
                   list(mean = ~if (cur_column() == "rmse_sig")  mean(.x, na.rm = TRUE) * 1e3
                     else if (cur_column() == "rmse_Betas") mean(.x, na.rm = TRUE) * 1e2
                     else                                    mean(.x, na.rm = TRUE),
                     sd   = ~if      (cur_column() == "rmse_sig")  sd(.x, na.rm = TRUE) * 1e3
                     else if (cur_column() == "rmse_Betas") sd(.x, na.rm = TRUE) * 1e2
                     else   sd(.x, na.rm = TRUE)), .names = "{.col}_{.fn}"), .groups = "drop")

# 2. Construct interleaved table: first the means, then the sds in parentheses
num_vars <- names(results_all)[sapply(results_all, is.numeric)]
# Remove any 'n' column or any numeric col you don't want to report (optional).
num_vars <- setdiff(num_vars, c("n", "rmse_lambda", "time", "iter", "Precision", "Sensitivity"))
num_vars <- num_vars[c(1, 6, 2, 3, 4, 5)]

# Get column names for means and sds
mean_cols <- paste0(num_vars, "_mean")
sd_cols   <- paste0(num_vars, "_sd")

# Optional: rename display variables for table
display_vars <- sub("^rmse_", "", num_vars)

# 3. Build mean and sd (parentheses) data.frames
means_df <- results_summary %>%
  dplyr::select(model2, Scenario, all_of(mean_cols))
colnames(means_df) <- c("Model", "Scenario", display_vars)
means_df[display_vars] <- lapply(means_df[display_vars], function(x) sprintf("%.2f", x))

sds_df <- results_summary %>%
  dplyr::select(model2, Scenario, all_of(sd_cols))
colnames(sds_df) <- c("Model", "Scenario", display_vars)
sds_df$Model <- ""  # remove model name for sd row
# Format sd: as (0.32), not plain 0.32
sds_df[display_vars] <- lapply(sds_df[display_vars], function(x) sprintf("~{\\scriptsize (%.2f)}", x))

# 4. Interleave means and sds rows
table_long <- lapply(seq_len(nrow(means_df)), function(i) {
  rbind(means_df[i,], sds_df[i,])
}) %>% do.call(rbind, .)

# Optional: filter (e.g. only one scenario)
names_cols <- c("\\textsc{model}", "$\\hat{K}$", "$F_1$", "$R$", "$\\Theta([0, T])$", "$B$", "$N_b$")
table_longA <- table_long %>% filter(Scenario == "Scenario_A_indep") %>% dplyr::select(-Scenario)
table_longB <- table_long %>% filter(Scenario == "Scenario_B_corr") %>% dplyr::select(-Scenario)
colnames(table_longA) <- colnames(table_longB) <- names_cols

# 5. Ready to output
kable(cbind(table_longA, table_longB[,-1]),
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      align = c("l", rep("c", ncol(table_long)-1)),
      caption = "Means and standard deviations (in parentheses) for each model. SD rows unlabelled.")

#------------------------------------------------------
# Make a plot now
#------------------------------------------------------
results_all2 <- results_all %>%
  mutate(model2 = case_when(
    model == "CompNMFBase" ~ "4. MAP, CompNMF",
    model == "map_TrueCovs" ~ "1. MAP, true x",
    model == "map_Full" ~ "2. MAP, all x",
    model == "mcmc_Full" ~ "3. MCMC, all x",
    TRUE ~ NA
  )) %>%
  filter(!is.na(model2))

results_all2 %>%
  mutate(Lambda = rmse_lambda,
         Signatures = rmse_sig,
         Theta = rmse_theta,
         Betas = case_when(model2 == "4. MAP, CompNMF" ~ NA,
                           TRUE~rmse_Betas),
         Scenario2 = case_when(Scenario == "Scenario_A_indep" ~ "A - indep",
                               TRUE ~ "B - corr")) %>%
  dplyr::select(model2, Scenario2, Lambda:Betas) %>%
  gather(key = "key", value = "value", -model2, -Scenario2) %>%
  mutate(key = factor(key, levels = c("Signatures", "Theta", "Betas", "Lambda"))) %>%
  ggplot()+
  geom_boxplot(aes(x = Scenario2, y = log(value), fill = model2, color = model2), alpha = 0.6) +
  facet_wrap(~key, scales = "free", nrow = 1) +
  scale_fill_manual(name = "Model", values = c("#CD2626", "#000D8B", "#89BBF6", "antiquewhite")) +
  scale_color_manual(name = "Model", values = c("#CD2626", "#000D8B", "#89BBF6", "antiquewhite3")) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab("Scenario") +
  ylab("log RMSE")
ggsave("../figures/Simualations_results.pdf", width = 9.70, height=2.34)


#------------------------------------------------------
# Make a plot now
#------------------------------------------------------
results_all3 <- results_all %>%
  mutate(model2 = case_when(
    model == "map_TrueCovs" ~ "0. MAP, true x",
    model == "CompNMFBase" ~ "1. MAP, CompNMF",
    model == "map_CopyOnly" ~ "2. MAP, copy only",
    model == "map_Full" ~ "3. MAP, all x",
    model == "mcmc_Full" ~ "4. MCMC, all x",
    model == "mcmc_TrueCovs200" ~ "5. MCMC, all x, Delta = 200",
    model == "mcmc_TrueCovs500" ~ "6. MCMC, all x, Delta = 500",
    TRUE ~ NA
  )) %>%
  filter(!is.na(model2))

colors <- c("#CD2626","antiquewhite3","#FF8C00", "#000D8B", "#89BBF6", "lightgreen", "forestgreen")
results_all3 %>%
  mutate(K = Kest,
         Lambda = log(rmse_lambda),
         Signatures = log(rmse_sig),
         Theta = log(rmse_theta),
         Betas = case_when(model2 == "1. MAP, CompNMF" ~ NA,
                           model2 == "2. MAP, copy only" ~ NA,
                           TRUE~log(rmse_Betas)),
         Scenario2 = case_when(Scenario == "Scenario_A_indep" ~ "A - indep",
                               TRUE ~ "B - corr")) %>%
  dplyr::select(model2, Scenario2, K, F1, Lambda:Betas) %>%
  gather(key = "key", value = "value", -model2, -Scenario2) %>%
  mutate(key = factor(key, levels = c("K", "F1", "Lambda",
                                      "Signatures", "Theta", "Betas"))) %>%
  ggplot()+
  geom_boxplot(aes(x = Scenario2, y = value, fill = model2, color = model2), alpha = 0.6) +
  facet_wrap(~key, scales = "free", nrow = 2) +
  scale_fill_manual(name = "Model", values = colors) +
  scale_color_manual(name = "Model", values = colors) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab("Scenario")+
  theme(axis.title.y = element_blank())
ggsave("../figures/Simualations_results_Supplement.pdf",
       width = 9.04, height=5.20)

# Finally, let's calculate average effective sample sizes and iterations
# Now, Average effective sample sizes and number of iterations

results_summary$effectiveSigs_mean


results_summary <- results_all %>%
  mutate(model2 = case_when(
    model == "map_TrueCovs" ~ "0. MAP, true x",
    model == "CompNMFBase" ~ "1. MAP, CompNMF",
    model == "map_CopyOnly" ~ "2. MAP, copy only",
    model == "map_Full" ~ "3. MAP, all x",
    model == "mcmc_Full" ~ "4. MCMC, all x",
    model == "mcmc_TrueCovs200" ~ "5. MCMC, all x",
    model == "mcmc_TrueCovs500" ~ "6. MCMC, all x",
    TRUE ~ NA
  )) %>%
  filter(!is.na(model2)) %>%
  group_by(model2, Scenario) %>%
  summarise(across(where(is.numeric),
                   list(mean = ~if (cur_column() == "rmse_sig")  mean(.x, na.rm = TRUE) * 1e3
                        else if (cur_column() == "rmse_Betas") mean(.x, na.rm = TRUE) * 1e2
                        else                                    mean(.x, na.rm = TRUE),
                        sd   = ~if      (cur_column() == "rmse_sig")  sd(.x, na.rm = TRUE) * 1e3
                        else if (cur_column() == "rmse_Betas") sd(.x, na.rm = TRUE) * 1e2
                        else   sd(.x, na.rm = TRUE)), .names = "{.col}_{.fn}"), .groups = "drop")%>%
  dplyr::select(model2, Scenario, time_mean:effectiveLogPost_sd) %>%
  as.data.frame()

# Summary for MAP (time and n.iterations)

results_summary %>%
  filter(grepl("MAP", model2)) %>%
  dplyr::select(model2:iter_sd)


results_summary %>%
  filter(grepl("MCMC", model2)) %>%
  dplyr::select(model2, time_mean, time_sd, effectiveBetas_mean:effectiveLogPost_sd)

# Summary for MCMC (time and eff sample sizes)



results_summary %>%
  mutate(model2 = case_when(
    model == "map_TrueCovs" ~ "0. MAP, true x",
    model == "CompNMFBase" ~ "1. MAP, CompNMF",
    model == "map_CopyOnly" ~ "2. MAP, copy only",
    model == "map_Full" ~ "3. MAP, all x",
    model == "mcmc_Full" ~ "4. MCMC, all x",
    model == "mcmc_TrueCovs200" ~ "5. MCMC, all x",
    model == "mcmc_TrueCovs500" ~ "6. MCMC, all x",
    TRUE ~ NA
  )) %>%
  filter(!is.na(model2)) %>%
  dplyr::select(model2, Scenario, time_mean:effectiveLogPost_sd) %>%
  as.data.frame()

results_summary %>%
  group_by(Scenario) %>%
  summarise(m = mean(n_mean)) %>%
  pull(m)








results_all3 <- results_all2 %>%
  mutate(Lambda = rmse_lambda,
         Signatures = rmse_sig,
         Theta = rmse_theta,
         Betas = rmse_Betas,
         Scenario = case_when(Scenario == "Scenario_A_indep" ~ "A: indep",
                               TRUE ~ "B: corr"))

# Plot for Beta
pBeta <- ggplot(results_all3) +
  geom_boxplot(aes(x = Scenario, y = log(Betas), fill = model2, color = model2),
               alpha = 0.5) +
  facet_grid(~"logRMSE Beta") +
  scale_fill_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  scale_color_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  theme_bw() + theme(axis.title = element_blank())

# Plot for Theta # <----- Consider taking the log??
pTheta <- ggplot(results_all3) +
  geom_boxplot(aes(x = Scenario, y = log(Theta), fill = model2, color = model2),
               alpha = 0.5) +
  facet_grid(~"logRMSE Theta") +
  scale_fill_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  scale_color_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  theme_bw() + theme(axis.title = element_blank())

# Signatures
pSigs <- ggplot(results_all3) +
  geom_boxplot(aes(x = Scenario, y = log(Signatures), fill = model2, color = model2),
               alpha = 0.5) +
  facet_grid(~"logRMSE R") +
  scale_fill_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  scale_color_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  theme_bw() + theme(axis.title = element_blank())

# F1
pF1 <- ggplot(results_all3) +
  geom_boxplot(aes(x = Scenario, y = F1, fill = model2, color = model2),
               alpha = 0.5) +
  facet_grid(~"F1 score") +
  scale_fill_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  scale_color_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  theme_bw() + theme(axis.title = element_blank())

# Kest
pK <- ggplot(results_all3) +
  geom_hline(yintercept = 8, linetype = "dashed") +
  geom_boxplot(aes(x = Scenario, y = Kest, fill = model2, color = model2),
               alpha = 0.5) +
  facet_grid(~"Estimated K") +
  scale_fill_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  scale_color_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  theme_bw() +
  ylim(c(7, 10)) + theme(axis.title = element_blank())

# Lambda
pLambda <- ggplot(results_all3) +
  geom_boxplot(aes(x = Scenario, y = log(Lambda), fill = model2, color = model2),
               alpha = 0.5) +
  facet_grid(~"logRMSE Lambda") +
  scale_fill_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  scale_color_manual(name = "Model", values = c("red", "blue", "darkorange", "forestgreen")) +
  theme_bw() + theme(axis.title = element_blank())

(pK + pF1 + pLambda) / (pSigs + pTheta + pBeta) + plot_layout(guides = 'collect')



# n <- get_counts_from_data(data)
#
#
#
# # STEP 1 - LOAD THE DATA
# data <- readRDS(paste0(out_dir, "data.rds.gzip"))
#
# data <- sim
# tilewidth <- data$simulation_parameters$tilewidth
# J <- data$simulation_parameters$J
#
# data <- generate_MutationData(J = 40, mu_copy = 0.5, size_copy = 20,
#                               ncores = ncores_simulation,
#                               corr = "indep")
#
# # SignalTrack and copy number track
# SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
# CopyTrackSim <- tilewidth * as.matrix(mcols(data$gr_CopyTrack)) / 2
#
# # Calculate the true intensity function
# Lambda_true <- Reconstruct_Lambda(SignalTrack = SignalTrackSim,
#                                   CopyTrack = CopyTrackSim,
#                                   R = data$R, Theta = data$Theta, Betas = data$Betas)
# # Do not adjust for copy number
# CopyTrackNull <- matrix(tilewidth, nrow = nrow(SignalTrackSim), ncol = 40)
# colnames(CopyTrackNull) <- colnames(CopyTrackSim)
#
# Lambda_null <- Reconstruct_Lambda(SignalTrack = SignalTrackSim,
#                                   CopyTrack = CopyTrackNull,
#                                   R = data$R, Theta = data$Theta, Betas = data$Betas)
#
# plot(rowSums(Lambda_true), type = "l")
# lines(rowSums(Lambda_null), type = "l", col = "red")
# plot(rowSums(CopyTrackSim/100))
# plot(rowSums(CopyTrackSim/100))
# table(CopyTrackSim/100)
#
# plot(Lambda_null[, 1], Lambda_true[, 1])
#
# j <- 33
# plot(scale(Lambda_true[, j]), type = "l")
# lines(scale(Lambda_null[, j]), type = "l", col = "red")
# lines(scale(CopyTrackSim[, j]), type = "l", col = "blue")
#
# plot(Lambda_true[, j], Lambda_null[, j])
# abline(a=0, b = 1)
#
# points(CopyTrackSim[, 10]/100 - 1)
# plot(CopyTrackSim[, 10]/100 - 1)
# plot(CopyTrackSim[, 10],)
# boxplot( Lambda_true[, j] ~ CopyTrackSim[, j])
#
# data <- readRDS(paste0(out_dir, "data.rds.gzip"))
# CopyTrack <- data$simulation_parameters$tilewidth * as.matrix(mcols(data$gr_CopyTrack)) / 2
# SignalTrack <- as.matrix(mcols(data$gr_SignalTrack)[, -1])
#
#
#
#
#
# dd <- generate_CopyTrack2()
# plot(dd$Patient_01)
# plot(dd$Patient_20)
# plot(dd$Patient_04)
#
# apply(1:MatchedSigs$R_hat - MatchedSigs$R_true)
#
#
# length(data$gr_Mutations)
# plot(MatchedSigs$R_hat[, 10], MatchedSigs$R_true[, 10])
#
# CompressiveNMF::plot_SBS_signature(MatchedSigs$R_hat)
#
#
#


# res <- readRDS(paste0(out_dir, "output_mcmc_TrueCovs.rds.gzip"))
#
# # Exclude the signatures that are redundant in the factorization
# cutoff <- 5 * res$PriorParams$a * res$PriorParams$epsilon
# filter_mu <- res$Mu > 5 * cutoff
# filter_sigs <- c(sigminer::cosine(res$Signatures, matrix(rep(1/96, 96))) < 0.975)
# round(diag(crossprod(res$Betas)), 5)
#
# # Estimate quantities in each model
# # Estimated quantities after filtering
# R_hat <- res$Signatures[, filter_mu & filter_sigs]
# Beta_hat <- res$Betas[, filter_mu & filter_sigs]
# theta_baseline <- res$Thetas[filter_mu & filter_sigs, ]
# Theta_total <- theta_baseline * crossprod(exp(SignalTrack[, rownames(Beta_hat)] %*% Beta_hat), CopyTrack)
#
# Lambda_hat <- Reconstruct_Lambda(SignalTrack = SignalTrack[, rownames(Beta_hat)], CopyTrack = CopyTrack,
#                                  R = R_hat, Theta = theta_baseline, Betas = Beta_hat)

# data$gr_Mutations
#
#
# #------- Generate data
# set.seed(24, kind = "L'Ecuyer-CMRG")
# data <- generate_MutationData(J = 40, sd_beta = 1/2)
# saveRDS(data, "~/SigPoisProcess/data/temp_simulation_multiplicativeCN.rds.gzip", compress="gzip")
# data <- readRDS("~/SigPoisProcess/data/temp_simulation_multiplicativeCN.rds.gzip")
#
# #---- Plot the number of mutation in each region
# length(data$gr_Mutations)
# corrplot(cor(as.matrix(mcols(data$gr_SignalTrack))[, -1]))
# corrplot(cor(as.matrix(mcols(data$gr_Mutations))[, -c(1,2)]))
#
# #-- Intensity function
# grX <- data$gr_SignalTrack
# CopyTrack <- as.matrix(mcols(data$gr_CopyTrack))
# Xcov <- as.matrix(mcols(grX)[, -1])
# ExpTrack <- exp(Xcov %*% data$Betas)
# LambdaAll <- rep(0, length(grX))
# for(j in 1:ncol(data$Theta)){
#   print(j)
#   ExpTrack_cj <- apply(ExpTrack, 2, function(x) width(grX) * 0.5 * CopyTrack[, j] * x)
#   for(i in 1:96){
#     LambdaAll <- LambdaAll +  c(crossprod(data$R[i, ], data$Theta[, j] * t(ExpTrack_cj)))
#   }
# }
# plot(LambdaAll, type = "l")
#
# #-- Number of mutation in each bin
# df_tmp <- as.data.frame(data$gr_Mutations)
# over <- findOverlaps(data$gr_Mutations, data$gr_SignalTrack)
# df_tmp$region <- subjectHits(over)
# df_all <- df_tmp %>%
#   group_by(region)%>%
#   summarize(n = n())
# MutAll <- rep(0, length(grX))
# MutAll[df_all$region] <- df_all$n
# plot(MutAll, type = "l")
#
# plot(MutAll, type = "l")
# lines(LambdaAll, col = "red")
#
#
# mat <- getTotalMutations(data$gr_Mutations)
#
#
# # Model 0 - Baseline Poisson NMF
#
# # Model 1 - Copy number only
#
# # Model 2 - Junk Covariates only
#
# # Model 3 - Covariate 1
#
# # Model 4 - Covariate 1, 2
#
# # Model 5 - Covariate 1, 2, 3
#
# # Model 6 - Covariate 1, 2, 3, 4
#
# # Model 7 - Covariate 1, 2, 3, 4, 5
#
#
# set.seed(10)
# gr_MutationsSim <- data$gr_Mutations
# SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
# CopyTrackSim <- 100 * as.matrix(mcols(data$gr_CopyTrack)) / 2
#
# sample_channel <- c("sample", "channel")
# covariate_names <- colnames(SignalTrackSim)
# junk_covs <- covariate_names[6:10]
# true_covs <- covariate_names[1:5]
#
#
# # Select the model with only junk covariates
# cols_selected <- c(sample_channel, true_covs[1:5], junk_covs)
# gr_Mutations_temp <- data$gr_Mutations
# mcols(gr_Mutations_temp) <- mcols(gr_Mutations_temp)[cols_selected]
# SignalTrack_temp <- SignalTrackSim[, c(true_covs[1:5], junk_covs)]
#
# set.seed(10)
# out <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_temp,
#                                  SignalTrack = SignalTrack_temp,
#                                  CopyTrack = CopyTrackSim,
#                                  K = 12,
#                                  method = "map",
#                                  update_Signatures = TRUE,
#                                  update_Theta = TRUE,
#                                  update_Betas = TRUE,
#                                  controls = SigPoisProcess_mult.controls(maxiter = 2000, tol = 1e-6),
#                                  prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
#                                                                                 c0 =  100, d0 = 1))
#
#
# p1 <- CompressiveNMF::plot_SBS_signature(out$Signatures)
# p2 <- plot_betas(out$Betas)
# p1 + p2
#
# # do MCMC now
# nsamples <- 1000
# burnin <- 500
# controlsSim <- SigPoisProcess_mult.controls(nsamples = nsamples, burnin = burnin)
# # Set inital quantities
# initSim <- SigPoisProcess_mult.init(R_start = out$Signatures,
#                                     Theta_start = out$Thetas,
#                                     Betas_start = out$Betas,
#                                     Mu_start = out$Mu,
#                                     Sigma2_start = out$Sigma2)
#
#
# set.seed(42)
# #start_MCMC <- Sys.time()
# outMCMCSim <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_temp,
#                                     SignalTrack = SignalTrack_temp,
#                                     CopyTrack = CopyTrackSim,
#                                     K =12,
#                                     method = "mcmc",
#                                     init = initSim,
#                                     controls = controlsSim)
#
# initSim <- SigPoisProcess_mult.init(R_start = data$R)
# outMCMCSim3 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_temp,
#                                     SignalTrack = SignalTrack_temp,
#                                     CopyTrack = CopyTrackSim,
#                                     K =8,
#                                     method = "mcmc",
#                                     controls = controlsSim,
#                                     init = initSim,
#                                     prior_params = SigPoisProcess_mult.PriorParams(c0 = 100, d0 = 1))
#
# get_PosteriorEffectiveSize(outMCMCSim3$chain$BETASchain[-c(1:burnin), ,])
# ci_beta <- get_PosteriorCI(outMCMCSim3$chain$BETASchain[-c(1:burnin), ,])
# ci_sigs <- get_PosteriorCI(outMCMCSim3$chain$SIGSchain[-c(1:burnin), ,])
# ci_theta <- get_PosteriorCI(outMCMCSim3$chain$THETAchain[-c(1:burnin), ,])
#
# p_betas <- plot_betasCI(ci_beta$mean, ci_beta$lowCI, ci_beta$highCI)
# p_sigs <- CompressiveNMF:::plot_SBS_signature(ci_sigs$mean,
#                                               lowCI = ci_sigs$lowCI,
#                                               highCI = ci_sigs$highCI)
# p_sigs + p_betas
# plot(outMCMCSim3$chain$BETASchain[, "Encode08", "SigN01"])
#
# plot(outMCMCSim$chain$BETASchain[, "Encode06", "SigN07"], type = "l")
#
# hist(outMCMCSim$chain$BETASchain[, "Encode06", "SigN07"], breaks = 40)
#
# #---- Try the fixed example with wrong covariates
# set.seed(10)
# out1 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_temp,
#                              SignalTrack = SignalTrack_temp,
#                              CopyTrack = CopyTrackSim,
#                              K = 8,
#                              method = "map",
#                              update_Signatures = FALSE,
#                              update_Theta = TRUE,
#                              update_Betas = TRUE,
#                              init = SigPoisProcess_mult.init(R_start = data$R),
#                              controls = SigPoisProcess_mult.controls(maxiter = 500, tol = 1e-6),
#                              prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
#                                                                             c0 =  100, d0 = 1))
#
# set.seed(10)
# out2 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_temp,
#                               SignalTrack = SignalTrack_temp,
#                               CopyTrack = CopyTrackSim,
#                               K = 8,
#                               method = "map",
#                               update_Signatures = FALSE,
#                               update_Theta = TRUE,
#                               update_Betas = FALSE,
#                               init = SigPoisProcess_mult.init(R_start = data$R, Betas_start = matrix(0, 5, 8)),
#                               controls = SigPoisProcess_mult.controls(maxiter = 500, tol = 1e-6),
#                               prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
#                                                                              c0 =  100, d0 = 1))
#
# logs <- cbind("logPost"= c(out1$sol$trace[5], out2$sol$trace[5]),
#       "logLik"= c(out1$sol$trace_logLik[5], out2$sol$trace_logLik[5]),
#       "logPrior"= c(out1$sol$trace_logPrior[5], out2$sol$trace_logPrior[5]))
# rownames(logs) <- c("Beta Estimated", "Beta set to zero")
# logs
# rbind("logPost Beta Estimated"= out1$sol$trace[5],
#       "logPost Beta set to zero"= out2$sol$trace[5])
# out2$sol$trace[5]
#
# out1$sol$trace_logLik[5]
# out2$sol$trace_logLik[5]
#
# p1 <- CompressiveNMF::plot_SBS_signature(out1$Signatures)
# p2 <- plot_betas(out1$Betas)
# p1 + p2
#
# pp1 <- out1$Thetas * t(crossprod(CopyTrackSim, exp(SignalTrack_temp %*% out1$Betas)))
# pp2 <- out2$Thetas * t(crossprod(CopyTrackSim, exp(SignalTrack_temp %*% out2$Betas)))
# pp <- data$Theta * t(crossprod(CopyTrackSim, exp(SignalTrackSim[, 1:5] %*% data$Betas)))
#
# plot(pp1, pp, xlab = "Estimated Total Loading", ylab = "True total loading")
# abline(a = 0, b = 1)
# plot(pp, pp2)
#
# Lambda1 <- Reconstruct_Lambda(SignalTrack_temp, CopyTrackSim, R = out1$Signatures,
#                    Theta = out1$Thetas, Betas = out1$Betas)
# Lambda2 <- Reconstruct_Lambda(SignalTrack_temp, CopyTrackSim, R = out2$Signatures,
#                               Theta = out2$Thetas, Betas = out2$Betas)
# Lambda <- Reconstruct_Lambda(SignalTrackSim[, 1:5], CopyTrackSim, R = data$R,
#                               Theta = data$Theta, Betas = data$Betas)
#
# plot(rowSums(Lambda), col = "gray", type = "l")
# lines(rowSums(Lambda1), col = "red", type = "l")
# lines(rowSums(Lambda2), col = "blue")
#
# sqrt(mean((rowSums(Lambda) - rowSums(Lambda1))^2))
# sqrt(mean((rowSums(Lambda) - rowSums(Lambda2))^2))
#
# # Model 1 ---------- Try the model with the full set of covariates
# set.seed(10)
# gr_MutationsSim <- data$gr_Mutations
# SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
# CopyTrackSim <- 100 * as.matrix(mcols(data$gr_CopyTrack)) / 2
#
# #--------------
# set.seed(10)
# resFull <- SigPoisProcess_multCN(gr_Mutations = gr_MutationsSim,
#                                SignalTrack = SignalTrackSim,
#                                CopyTrack = CopyTrackSim,
#                                K = 12,
#                                method = "map",
#                                update_Signatures = TRUE,
#                                update_Theta = TRUE,
#                                update_Betas = TRUE,
#                                #init = SigPoisProcess_mult.init(R_start = resSol$Signatures),
#                                controls = SigPoisProcess_mult.controls(maxiter = 1000, tol = 1e-6),
#                                prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
#                                                                               c0 =  100, d0 = 1))
#
#
# p1 <- CompressiveNMF::plot_SBS_signature(resFull$Signatures)
# p2 <- plot_betas(resFull$Betas)
# p1 + p2
#
# #----- Run now the mcmc using map as as starting point
# nsamples <- 2000
# burnin <- 100
# # Set inital quantities
# initSim <- SigPoisProcess_mult.init(R_start = resFull$Signatures,
#                                  Theta_start = resFull$Thetas,
#                                  Betas_start = resFull$Betas,
#                                  Mu_start = resFull$Mu,
#                                  Sigma2_start = resFull$Sigma2)
# controlsSim <- SigPoisProcess_mult.controls(nsamples = nsamples, burnin = burnin)
#
# set.seed(42)
# #start_MCMC <- Sys.time()
# outMCMCSim <- SigPoisProcess_multCN(gr_Mutations = gr_MutationsSim,
#                                  SignalTrack = SignalTrackSim,
#                                  CopyTrack = CopyTrackSim,
#                                  K =12,
#                                  method = "mcmc",
#                                  init = initSim,
#                                  controls = controlsSim)
# #end_MCMC <- Sys.time()
# #outMCMC$time <- end_MCMC - start_MCMC
# plot(outMCMCSim$chain$BETASchain[, 8, 2], type = "l")
# effectiveSize(outMCMCSim$chain$MUchain[-c(1:burnin), ])
# plot(outMCMCSim$chain$MUchain[, 10], type = "l")
# plot(outMCMCSim$chain$logPostchain, type = "l")
#
#
#
# get_PosteriorEffectiveSize(outMCMCSim$chain$BETASchain[-c(1:burnin), ,])
# ci_beta <- get_PosteriorCI(outMCMCSim$chain$BETASchain[-c(1:burnin), ,])
# ci_sigs <- get_PosteriorCI(outMCMCSim$chain$SIGSchain[-c(1:burnin), ,])
# ci_theta <- get_PosteriorCI(outMCMCSim$chain$THETAchain[-c(1:burnin), ,])
#
# p_betas <- plot_betasCI(ci_beta$mean, ci_beta$lowCI, ci_beta$highCI)
# p_sigs <- CompressiveNMF:::plot_SBS_signature(ci_sigs$mean,
#                                               lowCI = ci_sigs$lowCI,
#                                               highCI = ci_sigs$highCI)
#
# p_sigs + p_betas
#
#
#
#
# # Now, let's do the Bayesian
# match <- match_MutSign(R_true = data$R, R_hat = resFull$Signatures)
# plot(resFull$sol$Betas[, match$match][1:5, 1:8], data$Betas)
# abline(a = 0, b = 1)
#
# mutMatrix <- getTotalMutations(gr_MutationsSim)
# MutEst <- Reconstruct_CountMatrix(SignalTrackSim, CopyTrackSim,
#                                   R = resFull$Signatures,
#                                   Theta = resFull$Thetas,
#                                   Betas = resFull$Betas)
#
# plot(MutEst, mutMatrix)
# abline(a = 0, b = 1)
#
#
# # Try the Bayesian solution now.
# nsamples <- 500
# burn_in <- 10
# R_start <- resFull$Signatures
# Theta_start <- resFull$Thetas
# Betas_start <- resFull$Betas
# Mu_start <- resFull$Mu
# Sigma2_start <- resFull$Sigma2
# SigPrior <- resFull$PriorParams$SigPrior
# X <- as.data.frame(GenomicRanges::mcols(gr_MutationsSim)) %>%
#   dplyr::select(where(is.numeric)) %>%
#   as.matrix()
# # Extract mutation data
# Mutations <- as.data.frame(GenomicRanges::mcols(gr_MutationsSim)) %>%
#   dplyr::select(sample, channel)
# Mutations$sample <- as.factor(Mutations$sample) #, levels = unique(Mutations$sample))
# Mutations$channel <- as.factor(Mutations$channel)
#
# # Extract channel ID and mutation ID
# channel_id <- as.numeric(Mutations$channel)
# sample_id <- as.numeric(Mutations$sample)
#
# outBayes2 <- .PoissonProcess_BayesCN(R_start,
#                                    Theta_start,
#                                    Betas_start,
#                                    Mu_start,
#                                    Sigma2_start,
#                                    X,
#                                    SignalTrackSim,
#                                    CopyTrackSim,
#                                    nsamples,
#                                    burn_in,
#                                    channel_id - 1,
#                                    sample_id - 1,
#                                    SigPrior,
#                                    update_R = TRUE,
#                                    update_Theta = TRUE,
#                                    update_Betas = TRUE,
#                                    a = resFull$PriorParams$a,
#                                    a0 = resFull$PriorParams$a0,
#                                    b0 = resFull$PriorParams$b0,
#                                    c0 = 100,
#                                    d0 = 1)
#
# plot(outBayes2$logPostchain)
# plot(outBayes2$logLikchain)
# plot(outBayes2$logPriorchain)
#
#
# library(coda)
# plot(outBayes$THETAchain[, 2, 15], type = "l")
# plot(outBayes$MUchain[, 12], type = "l")
# effectiveSize(outBayes$MUchain)
#
# round(colMeans(outBayes$MUchain), 5)
# plot(hist(outBayes$BETASchain[, 7, 1], breaks = 40))
# abline(v = 0)
# effectiveSize(outBayes$BETASchain[-c(1:500), 5, 12])
# plot(outBayes$BETASchain[, 1, 10], type = "l")
# plot(outBayes$SIGMA2chain[, 12], type = "l")
# mean(outBayes$BETASchain[, 5, 12])
#
# round(colMeans(outBayes$SIGMA2chain), 5)
# plot(outBayes$SIGMA2chain[, 2])
# plot(outBayes$THETAchain[, 1, 2])
#
# chain <- outBayes$BETASchain
# dimnames(chain) <- list(NULL, rownames(Betas_start), colnames(Betas_start))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
