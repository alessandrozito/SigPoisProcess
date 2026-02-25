# Set the number of threads
library(RhpcBLASctl)
blas_set_num_threads(12)
omp_set_num_threads(1)
################################################################################
# This file runs the pre-specified signatures application of Section 5.
# The final files are:
#    -  resultsMCMC_denovo.rds.gzip: MCMC posterior summaries
#    -  postBurnin_MCMC_denovo.rds.gzip: MCMC chain post burnin
# They will be used to extract figures and outputs.
################################################################################

#---- Load packages
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(foreach)
library(patchwork)
library(doParallel)

# Load SigPoisProcess
library(SigPoisProcess)

################################################################################
# Step 1 - Load the data
################################################################################

# We source a separate script to load the data
source("~/SigPoisProcess/R/Load_ICGC_BreastAdenoCA_data.R")
# Load the cosmic signatures (v3.4)
Sigs_to_use <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8", "SBS20",
                 "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")
CosmicSigs <- SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37[, Sigs_to_use]

################################################################################
# Step 2 - Run the model
################################################################################

# BEWARE: THIS TAKE SEVERAL HOURS
out_dir <- "~/SigPoisProcess/output/Application/Breast_refit/"

rerun_MCMC_fixedSigs <- FALSE
if(rerun_MCMC_fixedSigs) {
  #------ MCMC
  nsamples_all <- 10000
  nsamples_current <- 0
  nsamples <- 100
  # Run the model storing the output every 100 iterations
  file_name <- paste0(out_dir, "MCMCSolution_SemiSupervised_noMap.rds.gzip")
  print("Start from random")
  set.seed(10)
  while(nsamples_current < nsamples_all) {
    print(paste0("<---- N. samples: ", nsamples_current, " ----->"))
    # Initialize either from the MAP or from the last previous iteration
    if(!file.exists(file_name)) {
      init <- SigPoisProcess.init(R_start = CosmicSigs)
      controls <- SigPoisProcess.controls(nsamples = nsamples, burnin = 10)
    } else {
      # Load MCMC chain
      outMCMC <- readRDS(file_name)
      id <- dim(outMCMC$chain$SIGSchain)[1]
      init <- SigPoisProcess.init(R_start = outMCMC$chain$SIGSchain[id, , ],
                                  Theta_start = outMCMC$chain$THETAchain[id, , ],
                                  Betas_start = outMCMC$chain$BETASchain[id, , ],
                                  Mu_start = outMCMC$chain$MUchain[id, ],
                                  Sigma2_start = outMCMC$chain$SIGMA2chain[id, ])
      controls <- SigPoisProcess.controls(nsamples = nsamples, burnin = 10) # Ignore the burnin part here. We simply update the chains
    }

    # Run the model for nsamples iterations
    out_temp <-  SigPoisProcess(gr_Mutations = data$gr_Mutations,
                                SignalTrack = data$SignalTrack,
                                CopyTrack = data$CopyTrack,
                                K = ncol(init$R_start),
                                verbose = TRUE,
                                method = "mcmc",
                                init = init,
                                controls = controls,
                                update_Signatures = FALSE)

    # Update the output
    if(!file.exists(file_name)){
      outMCMC <- out_temp
    } else {
      # Parameters
      outMCMC$chain$SIGSchain <- abind::abind(outMCMC$chain$SIGSchain, out_temp$chain$SIGSchain, along = 1)
      outMCMC$chain$THETAchain <- abind::abind(outMCMC$chain$THETAchain, out_temp$chain$THETAchain, along = 1)
      outMCMC$chain$BETASchain <- abind::abind(outMCMC$chain$BETASchain, out_temp$chain$BETASchain, along = 1)
      outMCMC$chain$MUchain <- rbind(outMCMC$chain$MUchain, out_temp$chain$MUchain)
      outMCMC$chain$SIGMA2chain <- rbind(outMCMC$chain$SIGMA2chain, out_temp$chain$SIGMA2chain)
      # Trace
      outMCMC$chain$logPostchain <- c(outMCMC$chain$logPostchain, out_temp$chain$logPostchain)
      outMCMC$chain$logLikchain <- c(outMCMC$chain$logLikchain, out_temp$chain$logLikchain)
      outMCMC$chain$logPriorchain <- c(outMCMC$chain$logPriorchain, out_temp$chain$logPriorchain)
    }

    # Update the sample tracker
    nsamples_current <- nsamples_current + nsamples

    # Save the model as it is running
    saveRDS(outMCMC, file = file_name, compress = "gzip")
    # Print the logposterior
    png(filename= paste0(out_dir, "logposterior_noMap.png"))
    plot(outMCMC$chain$logPostchain, type = "l")
    dev.off()
  }
}

################################################################################
# Step 2 - Run the model
################################################################################
rerun_postprocessing <- FALSE # <--- set to TRUE
if(rerun_postprocessing) {
  burnin <- c(1:2000)
  # Load the output
  outMCMCFixed <- readRDS(paste0(out_dir, "MCMCSolution_SemiSupervised_noMap.rds.gzip"))

  #----- Calculate posterior means and 95% credible intervals
  # Betas and Thetas
  ThetaCisFix <- get_PosteriorCI(outMCMCFixed$chain$THETAchain[-burnin, , ])
  BetasCisFix <- get_PosteriorCI(outMCMCFixed$chain$BETASchain[-burnin, , ])
  colnames(BetasCisFix$mean) <- colnames(BetasCisFix$lowCI) <- colnames(BetasCisFix$highCI) <- colnames(CosmicSigs)
  rownames(ThetaCisFix$mean) <- rownames(ThetaCisFix$lowCI) <- rownames(ThetaCisFix$highCI) <- colnames(CosmicSigs)

  # Mu and sigma
  sigma2 <- colMeans(outMCMCFixed$chain$SIGMA2chain[-burnin, ])
  mu <- colMeans(outMCMCFixed$chain$MUchain[-burnin, ])

  MuCis <- list(mean = mu,
                lowCI = apply(outMCMCFixed$chain$MUchain[-burnin, ], 2, function(x) quantile(x, 0.025)),
                highCI = apply(outMCMCFixed$chain$MUchain[-burnin, ], 2, function(x) quantile(x, 0.975)))
  sigma2Cis <- list(mean = sigma2,
                    lowCI = apply(outMCMCFixed$chain$SIGMA2chain[-burnin, ], 2, function(x) quantile(x, 0.025)),
                    highCI = apply(outMCMCFixed$chain$SIGMA2chain[-burnin, ], 2, function(x) quantile(x, 0.975)))
  names(MuCis$mean) <-  names(MuCis$lowCI) <- names(MuCis$highCI) <- colnames(CosmicSigs)
  names(sigma2Cis$mean) <-  names(sigma2Cis$lowCI) <- names(sigma2Cis$highCI) <- colnames(CosmicSigs)

  # Save all in a single object
  resultsMCMC_refit <- list("Signatures" = CosmicSigs,
                             "Baselines" = ThetaCisFix,
                             "Betas" = BetasCisFix,
                             "Mu" = MuCis,
                             "sigma2" = sigma2Cis,
                             "logPost" = outMCMCFixed$chain$logPostchain,
                             "logPrior" = outMCMCFixed$chain$logPriorchain,
                             "logLik" = outMCMCFixed$chain$logLikchain)
  # Save the output for plotting
  saveRDS(object = resultsMCMC_refit, file = "~/SigPoisProcess/output/Application/Breast_refit/resultsMCMC_refit.rds.gzip", compress = "gzip")

  ###############################################################################
  # Save the post burnin chain
  postBurnChain_MCMC_refit <- list("Signatures" = CosmicSigs,
                                    "Baselines" = outMCMCFixed$chain$THETAchain[-burnin, , ],
                                    "Betas" = outMCMCFixed$chain$BETASchain[-burnin, , ],
                                    "Mu" = outMCMCFixed$chain$MUchain[-burnin, ],
                                    "sigma2" = outMCMCFixed$chain$SIGMA2chain[-burnin, ],
                                    "logPost" = outMCMCFixed$chain$logPostchain[-burnin],
                                    "logPrior" = outMCMCFixed$chain$logPriorchain[-burnin],
                                    "logLik" = outMCMCFixed$chain$logLikchain[-burnin])
  saveRDS(object = postBurnChain_MCMC_refit, file = "~/SigPoisProcess/output/Application/Breast_refit/postBurnChain_MCMC_refit.rds.gzip", compress = "gzip")

}


