# Set the number of threads
library(RhpcBLASctl)
blas_set_num_threads(12)
omp_set_num_threads(1)
################################################################################
# This file runs the denovo application of Section 5.
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
# Step 0 - useful functions
################################################################################
create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

#---- effective sample sizes
get_PosteriorEffectiveSize <- function(chain, burnin = 1500){
  if(is.null(dim(chain))){
    coda::effectiveSize(chain[-c(1:burnin)])
  } else if(length(dim(chain)) == 2) {
    apply(chain[-c(1:burnin), ], 2, function(x) coda::effectiveSize(x))
  } else if (length(dim(chain)) == 3){
    apply(chain[-c(1:burnin), ,], c(2,3), function(x) coda::effectiveSize(x))
  }
}


#---- Starting values for the model parameters
generate_SigPoisProcess_init <- function(data,
                                         K = 12,
                                         prior_params = SigPoisProcess.PriorParams(),
                                         nnls = FALSE,
                                         sigStartNNLS = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8",
                                                          "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")){
  # Useful quantities
  J <- ncol(data$CopyTrack)
  I <- length(unique(data$gr_Mutations$channel))
  p <- ncol(data$SignalTrack)
  prior_params <- do.call("SigPoisProcess.PriorParams", prior_params)
  a <- prior_params$a
  alpha <- prior_params$alpha
  epsilon <- prior_params$epsilon
  if (is.null(prior_params$a0) | is.null(prior_params$b0)) {
    a0 <- a * J + 1
    b0 <- epsilon * a * J
  } else {
    a0 <- prior_params$a0
    b0 <- prior_params$b0
  }
  c0 <- prior_params$c0
  d0 <- prior_params$d0

  # Extract the matrix of aggregated mutations
  MutMatrix <- getTotalMutations(data$gr_Mutations)
  if(any(colnames(MutMatrix) != colnames(data$CopyTrack))){
    stop("colnames(MutMatrix) != colnames(data$CopyTrack)")
  }
  #----------------------------------------------- Prior for the signatures
  SigPrior <- prior_params$SigPrior
  if (!is.null(SigPrior)) {
    if (!is.matrix(SigPrior)) {
      stop("SigPrior must be a matrix")
    }
    Kpre <- ncol(SigPrior)
    if (K > 0) {
      SigPrior_new <- matrix(alpha, nrow = I, ncol = K)
      colnames(SigPrior_new) <- paste0("SigN", sprintf("%02d", 1:K))
      rownames(SigPrior_new) <- rownames(SigPrior)
      SigPrior <- cbind(SigPrior, SigPrior_new)
      K <- Kpre + K
    } else {
      K <- Kpre
    }
  } else {
    SigPrior <- matrix(alpha, nrow = I, ncol = K)
    #colnames(SigPrior) <- paste0("SigN", sprintf("%02d", 1:K))
    #rownames(SigPrior) <- levels(Mutations$channel)
  }

  #----------------------------------------------- Initial values
  Area_patients <- colSums(data$CopyTrack)
  if(isFALSE(nnls)){
    # Sample Signatures
    R_start <- sample_signatures(SigPrior)
    # Sample baseline activities
    Theta_start <- matrix(rgamma(J * K, colSums(MutMatrix)), nrow = K, byrow = TRUE)#/sum(bin_weight)
    Theta_start <- Theta_start / t(Area_patients)[rep(1, K), ]
    Mu_start <- rowMeans(Theta_start)
    Sigma2_start <- rep(d0/(c0 + 1), K)
    Betas_start <- matrix(rnorm(p * K, mean = 0, sd = 1e-7), nrow = p, ncol = K)
  } else {
    K <- length(sigStartNNLS)
    Cosmic_Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
    R_start <- Cosmic_Sigs[, sigStartNNLS]
    Theta_start <- sapply(1:J, function(j) nnls::nnls(R_start, MutMatrix[, j])$x + 0.001)
    Theta_start <- Theta_start / t(Area_patients)[rep(1, K), ]
    Mu_start <- rowMeans(Theta_start)
    Sigma2_start <- rep(d0/(c0 + 1), K)
    Betas_start <- matrix(rnorm(p * K, mean = 0, sd = 1e-7), nrow = p, ncol = K)
  }

  return(list("R_start" = R_start,
              "Theta_start" = Theta_start,
              "Betas_start" = Betas_start,
              "Mu_start" = Mu_start,
              "Sigma2_start" = Sigma2_start))
}
#------

################################################################################
# Step 1 - Load the data
################################################################################

# We source a separate script to load the data
source("~/SigPoisProcess/R/Load_ICGC_BreastAdenoCA_data.R")

################################################################
# Step 2 - Run the de-novo analysis
################################################################

# Control parameters for the MAP
maxiter <- 5000
controls <- SigPoisProcess.controls(maxiter = maxiter, tol = 1e-7)
# Number of signatures
K <- 12

n_start_point <- 3
# Generate the random starting point
regenerate_start_points <- FALSE
if(regenerate_start_points) {
  set.seed(10)
  for(i in 1:n_start_point){
    out_dir <- paste0("~/SigPoisProcess/output/Application/Breast_DeNovo/MapSolutions/Replicate", sprintf("%02d", i))
    create_directory(out_dir)
    init <- generate_SigPoisProcess_init(data = data, K = K)
    saveRDS(init, paste0(out_dir, "/init.rds"))
  }
}

#--------------------------------------------------------- MAXIMUM A POSTERIORI
#---- RUN the MAP now from the chosen starting point
nrepl <- 3 # 1, 2, and 3 <------ Change the value based on the replicate. 3 has the highest log posterior
out_dir <- paste0("~/SigPoisProcess/output/Application/Breast_DeNovo/MapSolutions/Replicate", sprintf("%02d", nrepl), "/")

rerun_MAP <- FALSE
if(rerun_MAP) {
  # Load initial value
  init <- readRDS(paste0(out_dir, "init.rds"))
  # Find the MAP
  print(paste0("Replicate ", nrepl))
  start_time <- Sys.time()
  outMAP <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
                                  SignalTrack = data$SignalTrack,
                                  CopyTrack = data$CopyTrack,
                                  K = ncol(init$R_start),
                                  verbose = TRUE,
                                  init = init,
                                  controls = controls)
  end_time <- Sys.time()
  outMAP$time <- end_time - start_time

  # Save it
  saveRDS(outMAP,
          file = paste0(out_dir, "MAPSolution.rds.gzip"),
          compress = "gzip")
}

#--------------------------------------------------------- MCMC FROM MAP
# NOTE: THIS FUNCTION TAKES SEVERAL HOURS. BE MINDFUL IF RE-RUNNING IT.
# WE SAVE EVERY 100 ITERATIONS TO AVOID ISSUES.

#---- RUN the MCMC now from the MAP solution
rerun_MCMC <- FALSE
if(rerun_MCMC){
  # Load the MAP
  outMAP <- readRDS(paste0(out_dir, "MAPSolution.rds.gzip"))
  #------ MCMC
  nsamples_all <- 10600
  nsamples_current <- 0
  nsamples <- 100
  # Run the model storing the output every 100 iterations
  file_name <- paste0(out_dir, "MCMCSolution.rds.gzip")
  set.seed(10)
  while(nsamples_current < nsamples_all) {
    print(paste0("<---- N. samples: ", nsamples_current, " ----->"))

    # Initialize either from the MAP or from the last previous iteration
    if(!file.exists(file_name)) {
      init <- SigPoisProcess.init(R_start = outMAP$Signatures,
                                       Theta_start = outMAP$Thetas,
                                       Betas_start = outMAP$Betas,
                                       Mu_start = outMAP$Mu,
                                       Sigma2_start = outMAP$Sigma2)
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
      controls <- SigPoisProcess.controls(nsamples = nsamples, burnin = 10)
    }

    # Run the model for nsamples iterations
    out_temp <-  SigPoisProcess(gr_Mutations = data$gr_Mutations,
                                       SignalTrack = data$SignalTrack,
                                       CopyTrack = data$CopyTrack,
                                       K = ncol(init$R_start),
                                       verbose = TRUE,
                                       method = "mcmc",
                                       init = init,
                                       controls = controls)

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
    png(filename= paste0(out_dir, "logposterior.png"))
    plot(outMCMC$chain$logPostchain, type = "l")
    dev.off()
  }
}

################################################################
# Step 3 - Save the estimates and the posterior credible intervals
################################################################
rerun_postprocessing <- FALSE # <--- set to TRUE
if(rerun_postprocessing) {
  #------ Load the MCMC output
  outMCMC <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate03/MCMCSolution.rds.gzip")

  ###############################################################################
  # Calculate mean and credible intervals
  burnin <- c(1:7600) # There are 7600 initial samples that I drew for burnin. Keep the last 3000

  sigma2 <- colMeans(outMCMC$chain$SIGMA2chain[-burnin, ])
  mu <- colMeans(outMCMC$chain$MUchain[-burnin, ])

  # Calculate posterior means and 95% credible intervals
  SigsCis <- get_PosteriorCI(outMCMC$chain$SIGSchain[-burnin, , ])
  ThetaCis <- get_PosteriorCI(outMCMC$chain$THETAchain[-burnin, , ])
  BetasCis <- get_PosteriorCI(outMCMC$chain$BETASchain[-burnin, , ])
  MuCis <- list(mean = mu,
                lowCI = apply(outMCMC$chain$MUchain[-burnin, ], 2, function(x) quantile(x, 0.025)),
                highCI = apply(outMCMC$chain$MUchain[-burnin, ], 2, function(x) quantile(x, 0.975)))
  sigma2Cis <- list(mean = sigma2,
                    lowCI = apply(outMCMC$chain$SIGMA2chain[-burnin, ], 2, function(x) quantile(x, 0.025)),
                    highCI = apply(outMCMC$chain$SIGMA2chain[-burnin, ], 2, function(x) quantile(x, 0.975)))

  # Save all in a single object
  resultsMCMC_denovo <- list("Signatures" = SigsCis,
                             "Baselines" = ThetaCis,
                             "Betas" = BetasCis,
                             "Mu" = MuCis,
                             "sigma2" = sigma2Cis,
                             "logPost" = outMCMC$chain$logPostchain,
                             "logPrior" = outMCMC$chain$logPriorchain,
                             "logLik" = outMCMC$chain$logLikchain)
  # Save the output for plotting
  saveRDS(object = resultsMCMC_denovo, file = "~/SigPoisProcess/output/Application/Breast_DeNovo/resultsMCMC_denovo.rds.gzip", compress = "gzip")

  ###############################################################################
  # Save the post burnin chain
  postBurnChain_MCMC_denovo <- list("Signatures" = outMCMC$chain$SIGSchain[-burnin, , ],
                                    "Baselines" = outMCMC$chain$THETAchain[-burnin, , ],
                                    "Betas" = outMCMC$chain$BETASchain[-burnin, , ],
                                    "Mu" = outMCMC$chain$MUchain[-burnin, ],
                                    "sigma2" = outMCMC$chain$SIGMA2chain[-burnin, ],
                                    "logPost" = outMCMC$chain$logPostchain[-burnin],
                                    "logPrior" = outMCMC$chain$logPriorchain[-burnin],
                                    "logLik" = outMCMC$chain$logLikchain[-burnin])
  saveRDS(object = postBurnChain_MCMC_denovo, file = "~/SigPoisProcess/output/Application/Breast_DeNovo/postBurnin_MCMC_denovo.rds.gzip", compress = "gzip")

}


###############################################################################
# Calculate effective sample sizes
mat_allCopy <- t(colSums(data$CopyTrack))[rep(1, sum(mu > 0.01)), ]
ThetaChain_adj <- outMCMC$chain$THETAchain[-c(1:7600), mu > 0.01, ]
array_AllCopy <- array(NA, dim = c(3000, 9, 113))
for(s in 1:3000) ThetaChain_adj[s,,] <- mat_allCopy * ThetaChain_adj[s, , ]

REffects <- get_PosteriorEffectiveSize(outMCMC$chain$SIGSchain[, , mu > 0.01], burnin = 7600)
ThetaEffects <- get_PosteriorEffectiveSize(ThetaChain_adj, burnin = 0)
BetasEffects <- get_PosteriorEffectiveSize(outMCMC$chain$BETASchain[, , mu > 0.01], burnin = 7600)
Sigma2Effects <- get_PosteriorEffectiveSize(outMCMC$chain$SIGMA2chain[, mu > 0.01], burnin = 7600)
MuEffects <- get_PosteriorEffectiveSize(outMCMC$chain$MUchain[, mu > 0.01], burnin = 7600)
lpEffects <- get_PosteriorEffectiveSize(outMCMC$chain$logPostchain, burnin = 7600)

#---- Calculate the posterior effective sample sizes
mat_allCopy <- t(colSums(data$CopyTrack))[rep(1, sum(mu > 0.01)), ]
ThetaChain_adj <- outMCMC$chain$THETAchain[-c(1:7600), mu > 0.01, ]
array_AllCopy <- array(NA, dim = c(3000, 9, 113))
for(s in 1:3000) ThetaChain_adj[s,,] <- mat_allCopy * ThetaChain_adj[s, , ]

REffects <- get_PosteriorEffectiveSize(outMCMC$chain$SIGSchain[, , mu > 0.01], burnin = 7600)
ThetaEffects <- get_PosteriorEffectiveSize(ThetaChain_adj, burnin = 0)
BetasEffects <- get_PosteriorEffectiveSize(outMCMC$chain$BETASchain[, , mu > 0.01], burnin = 7600)
Sigma2Effects <- get_PosteriorEffectiveSize(outMCMC$chain$SIGMA2chain[, mu > 0.01], burnin = 7600)
MuEffects <- get_PosteriorEffectiveSize(outMCMC$chain$MUchain[, mu > 0.01], burnin = 7600)
lpEffects <- get_PosteriorEffectiveSize(outMCMC$chain$logPostchain, burnin = 7600)

