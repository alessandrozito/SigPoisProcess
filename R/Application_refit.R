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
rerun_MCMC_fixedSigs <- FALSE
