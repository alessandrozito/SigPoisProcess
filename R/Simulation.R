################################################################################
# This file runs the simulation in Section 4 to reproduce the following:
# --- Figure 2: Results of the simulation
# --- Figure S1: Supplemental results of the simulation
# --- Table S1 and S2: time, iterations, and effective sample sizes
#                      for the simulation
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

# Load the SigPoisProcess package
library(SigPoisProcess)

#--- Source all the functions to run the simulation and evaluate performance
source("~/SigPoisProcess/R/Simulation_functions.R")

################################################################################
# Part 1 - Functions to run the simulation and postprocess the output
################################################################################

# Function to run the models, using the aggregation as simulated in the data
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
  controls <- SigPoisProcess.controls(maxiter = maxiter,
                                      tol = tol,
                                      nsamples = nsamples,
                                      burnin = burnin)
  # Prior Parameters for sigpoisprocess
  prior_params <- SigPoisProcess.PriorParams(c0 = c0, d0 = d0,
                                             epsilon = epsilon, a = a,
                                             alpha = alpha)
  Betas_start <- matrix(0, nrow = ncol(SignalTrackSim), ncol = K)

  ###########################
  # CompressiveNMF MAP, no covariates (MAP, CompNMF)
  ###########################
  start_time <- Sys.time()
  outCompNMFBase <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = K, a = a,
                                                       alpha = alpha, epsilon = epsilon, tol = tol)
  outCompNMFBase$time <- Sys.time() - start_time
  saveRDS(outCompNMFBase, file = paste0(out_dir, "output_CompNMFBase.rds.gzip"), compress = "gzip")

  ###########################
  # SigPoisProcess with only the copy number (MAP no x)
  ###########################
  start_time <- Sys.time()
  out_map_CopyOnly <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
                                     SignalTrack = SignalTrackSim,
                                     CopyTrack = CopyTrackSim,
                                     K = K,
                                     method = "map",
                                     update_Signatures = TRUE,
                                     update_Theta = TRUE,
                                     update_Betas = FALSE,
                                     controls = controls,
                                     prior_params = prior_params,
                                     init = SigPoisProcess.init(Betas_start = Betas_start))
  out_map_CopyOnly$time <- Sys.time() - start_time
  saveRDS(out_map_CopyOnly, file = paste0(out_dir, "output_map_CopyOnly.rds.gzip"), compress = "gzip")

  ###########################
  # SigPoisProcess with only the covariates, but no copy number. This is not
  # included in the main document
  ###########################

  # ------------ 3.1 - Maximum-a-posteriori
  start_time <- Sys.time()
  out_map_NoCopies <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
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

  ###########################
  # SigPoisProcess with only correct covariates (MAP, true x)
  ###########################
  SignalTrackSim_true <- SignalTrackSim[, 1:data$simulation_parameters$p_to_use]
  gr_Mutations_true <- data$gr_Mutations
  mcols(gr_Mutations_true) <- mcols(gr_Mutations_true)[, c("sample", "channel", colnames(SignalTrackSim_true))]
  # ------------ 4.1 - Maximum-a-posteriori
  start_time <- Sys.time()
  out_map_TrueCovs <- SigPoisProcess(gr_Mutations = gr_Mutations_true,
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

  ###########################
  # SigPoisProcess with also redundant covariates (MAP, all x and MCMC, all x)
  ###########################

  # ------------ 5.1 - Maximum-a-posteriori
  start_time <- Sys.time()
  out_map_Full <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
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
  init <- SigPoisProcess.init(R_start = out_map_Full$Signatures,
                              Theta_start = out_map_Full$Thetas,
                              Betas_start = out_map_Full$Betas,
                              Mu_start = out_map_Full$Mu,
                              Sigma2_start = out_map_Full$Sigma2)
  start_time <- Sys.time()
  out_mcmc_Full <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
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

# Function to run the models, using a coarser genomic aggregation for the covariates
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
  controls <- SigPoisProcess.controls(maxiter = maxiter,
                                      tol = tol,
                                      nsamples = nsamples,
                                      burnin = burnin)
  # Prior Parameters for sigpoisprocess
  prior_params <- SigPoisProcess.PriorParams(c0 = c0, d0 = d0,
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
    out_map_Full <- SigPoisProcess(gr_Mutations = data_agg$gr_Mutations,
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
    init <- SigPoisProcess.init(R_start = out_map_Full$Signatures,
                                Theta_start = out_map_Full$Thetas,
                                Betas_start = out_map_Full$Betas,
                                Mu_start = out_map_Full$Mu,
                                Sigma2_start = out_map_Full$Sigma2)
    start_time <- Sys.time()
    out_mcmc_Full <- SigPoisProcess(gr_Mutations = data_agg$gr_Mutations,
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
# Part 2 - Simulate all the data
################################################################################

# Set-up of the simuation
J <- 40             # N. of patients in each dataset
n_datasets <- 20    # N. of simulated datasets in each scenarios
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

#--- Run all models
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

reprocess_output <- FALSE #<--- Set to true if you want to re-run everything
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
  write_tsv(results_all, "~/SigPoisProcess/output/Simulation/simulation_results2.tsv")
}



################################################################################
# Part 5 - Draw the figures
################################################################################

# Load the results of the simulation (postprocessed)
results_all <- read_tsv("~/SigPoisProcess/output/Simulation/simulation_results2.tsv")

#--- Calculate F1 score
results_all <- results_all %>%
  mutate(F1 = 2 * Precision * Sensitivity / (Precision + Sensitivity))

#-------------------------------------------------------- Figure 2
results_all2 <- results_all %>%
  mutate(model2 = case_when(
    model == "CompNMFBase" ~ "4. MAP, CompNMF",
    model == "map_TrueCovs" ~ "1. MAP, true x",
    model == "map_Full" ~ "2. MAP, all x",
    model == "mcmc_Full" ~ "3. MCMC, all x",
    TRUE ~ NA
  )) %>%
  filter(!is.na(model2))

plot_Figure2 <- results_all2 %>%
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
ggsave(plot = plot_Figure2, filename = "~/SigPoisProcess/figures/Simulation_Figure2.pdf", width = 9.70, height=2.34)

#-------------------------------------------------------- Figure S2
results_all3 <- results_all %>%
  mutate(model2 = case_when(
    model == "map_TrueCovs" ~ "M0. MAP, true x",
    model == "CompNMFBase" ~ "M1. MAP, CompNMF",
    model == "map_CopyOnly" ~ "M2. MAP, no x",
    model == "map_Full" ~ "M3. MAP, all x",
    model == "mcmc_Full" ~ "M4. MCMC, all x",
    model == "mcmc_TrueCovs200" ~ "M5. MCMC, all x, Delta = 200",
    model == "mcmc_TrueCovs500" ~ "M6. MCMC, all x, Delta = 500",
    TRUE ~ NA
  )) %>%
  filter(!is.na(model2))

colors <- c("#CD2626","antiquewhite3","#FF8C00", "#000D8B", "#89BBF6", "lightgreen", "forestgreen")
plot_FigureS2 <- results_all3 %>%
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
ggsave(plot = plot_FigureS2, filename = "~/SigPoisProcess/figures/Simulation_FigureS2.pdf",
       width = 9.04, height=5.20)

#--------------------------------------------------------
# Time and effective sample sizes
results_summary <- results_all %>%
  mutate(model2 = case_when(
    model == "map_TrueCovs" ~ "M0. MAP, true x",
    model == "CompNMFBase" ~ "M1. MAP, CompNMF",
    model == "map_CopyOnly" ~ "M2. MAP, copy only",
    model == "map_Full" ~ "M3. MAP, all x",
    model == "mcmc_Full" ~ "M4. MCMC, all x",
    model == "mcmc_TrueCovs200" ~ "M5. MCMC, all x",
    model == "mcmc_TrueCovs500" ~ "M6. MCMC, all x",
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

#------------ Table S1 - computational time and iterations for the MAP
results_summary %>%
  filter(grepl("MAP", model2)) %>%
  dplyr::select(model2:iter_sd)

#------------ Table S2 - computational time and effective sample size for MCMC
results_summary %>%
  filter(grepl("MCMC", model2)) %>%
  dplyr::select(model2, time_mean, time_sd, effectiveBetas_mean:effectiveLogPost_sd)



