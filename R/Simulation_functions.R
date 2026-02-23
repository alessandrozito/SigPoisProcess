# This file contains all the functions to simulate the data in Section 4 of the
# manuscript

#--- Useful functions to open files and create directories
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
# Part 1 - Functions to simulate the data
################################################################################



###############################################################
#--------- Function to generate the SignalTrack
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

###############################################################
#--------- Function to generate the number of copies for each patient
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

###############################################################
#----- Function to generate the parameters (R, Theta, and Beta)
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

###############################################################
#----- Function to sample the mutations for a single patient j and a single mutation type i
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

###############################################################
#----- Function to sample the dataset in parallel
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


###############################################################
#----- Main function to simulate the whole mutation data
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


###############################################################
#----- Function to aggregate the signal of the covairates and the copy numbers
#      over the genome
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
# Part 2 - Functions to evaluate the output and compute performance
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



