# Sys.setenv(OMP_NUM_THREADS = 1)
# Sys.setenv(MKL_NUM_THREADS = 1)
# Sys.setenv(BLAS_NUM_THREADS = 1)
#Sys.setenv(OPENBLAS_NUM_THREADS = 1)
library(RhpcBLASctl)
blas_set_num_threads(12)
omp_set_num_threads(1)
################################################################################
# This file runs the Application in Section 5. It reproduces the following:
# --- Table 2: Perfomance comparison of the models under two scenarios
# --- Figure 2: Track on the signatures
# --- Table 1.2: xxx
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
devtools::document()

create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

#----- Useful functions
add_bin_weights <- function(gr){
  # Import blacklist and regions with gaps
  blacklist <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/hg19-blacklist.v2.bed")
  gap_gr <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/gaps_hg19.bed")

  # Adjust bin sizes to remove blacklist and gaps
  BinWeight <- width(gr)
  score <- gr$score
  # Remove gaps
  overlapGaps <- findOverlaps(gr, gap_gr)
  rangesGaps <- gr[queryHits(overlapGaps)]
  # Find the number of ACGT in the gap regions, and adjust the bin sizes to remove
  # the NNNNN sequences
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rangesGaps)
  countACGT <- rowSums((Biostrings::alphabetFrequency(seqs))[, c("A", "C", "G", "T")])
  BinWeight[queryHits(overlapGaps)] <- countACGT
  # Exlude now the blacklist as well
  overlap_blackList <- findOverlaps(gr, blacklist)
  pairs <- Pairs(gr[queryHits(overlap_blackList)], blacklist[subjectHits(overlap_blackList)])
  grIntersect <- pintersect(pairs)
  BinWeight[queryHits(overlap_blackList)] <- pmax(BinWeight[queryHits(overlap_blackList)] - width(grIntersect), 0)
  score[BinWeight == 0] <- 0 # Replace it with 0, knowing that we will remove it afterwards

  gr$score <- score
  gr$bin_weight <- BinWeight
  return(gr)
}

bin_gr <- function(gr, genome_aggreg, std_chrs){
  gr <- keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  gr <- binnedAverage(bins = keepSeqlevels(genome_aggreg, std_chrs, pruning.mode = "coarse"),
                      numvar = coverage(gr, weight = "score"),
                      varname = "score")
  return(gr)
}

build_CopyTrack <- function(gr_tumor, gr_copy, tilewidth = 2000){
  # Find the desired aggregation level
  genome_aggreg <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                              tilewidth = tilewidth,
                              cut.last.tile.in.chrom = TRUE)
  std_chrs <-  paste0("chr", c(1:22, "X"))
  # Add bin weights to the genome
  gr_CopyTrack <- add_bin_weights(genome_aggreg)
  genome_aggreg_w <- gr_CopyTrack
  std_chrs <-  paste0("chr", c(1:22, "X"))
  samples <- unique(gr_tumor$sample)
  # Filter copies in the tumor
  gr_copy_tumor <- gr_copy[gr_copy$sample %in% unique(gr_tumor$sample)]
  gr_copy_tumor$score[is.na(gr_copy_tumor$score)] <- 2
  # Output to store
  CopyTrack <- matrix(0, nrow = length(genome_aggreg), ncol = length(samples))
  for(s in seq_along(samples)) {
    gr_copy_sample <- bin_gr(gr_copy_tumor[gr_copy_tumor$sample == samples[s]],
                             genome_aggreg, std_chrs)
    # Substitute with 2 if NA. simplifying assumption.
    score <- gr_copy_sample$score
    #score[is.na(score)] <- 2
    score[score < 0.1] <- 0.1
    CopyTrack[, s] <- score/2
  }
  colnames(CopyTrack) <- samples
  mcols(gr_CopyTrack) <- cbind("bin_weight" = genome_aggreg_w$bin_weight, CopyTrack)
  return(gr_CopyTrack[genome_aggreg_w$bin_weight > 0])
}

build_SignalTrack <- function(tilewidth = 2000, type = "avg"){
  # Genome aggregated
  genome_aggreg <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                              tilewidth = tilewidth,
                              cut.last.tile.in.chrom = TRUE)
  std_chrs <-  paste0("chr", c(1:22, "X"))
  gr_SignalTrack <- add_bin_weights(genome_aggreg)

  # ---- GC content
  print("GC")
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$GC <- gr$score

  # ---- Methylation
  print("Methylation")
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Breast-Cancer_Methylation.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$Methyl <- gr$score

  # ---- CTCF + Histone marks programmatically handled
  marks <- c("CTCF", "H3K9me3", "H3K36me3", "H3K27me3", "H3K27ac", "H3K4me1", "H3K4me3")
  for (mark in marks) {
    print(mark)
    if(type == "avg") {
      # Calculate the average between tissue and cell line.
      # Tissue
      file_tissue <- paste0("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Breast-Cancer_tissue_", mark, "_2kb.bigWig")
      gr_tissue <- bin_gr(import(file_tissue), genome_aggreg, std_chrs)
      # Cell line
      file_cell <- paste0("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Breast-Cancer_cell_", mark, "_2kb.bigWig")
      gr_cell <- bin_gr(import(file_cell), genome_aggreg, std_chrs)
      # Average
      gr_SignalTrack$mark <- (gr_tissue$score + gr_cell$score)/2
      colnames(mcols(gr_SignalTrack))[ncol(mcols(gr_SignalTrack))] <- mark
    } else {
      file <- paste0("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Breast-Cancer_", type, "_", mark, "_2kb.bigWig")
      gr <- bin_gr(import(file), genome_aggreg, std_chrs)
      gr_SignalTrack$mark <- gr$score
      colnames(mcols(gr_SignalTrack))[ncol(mcols(gr_SignalTrack))] <- mark
    }
  }

  # ---- Replication Timing
  print("RepliTime")
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$RepliTime <- gr$score

  # ---- Nucleosome Occupancy
  print("NuclOccup")
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/NuclOccupancy/GSM920557_hg19_wgEncodeSydhNsomeK562Sig_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$NuclOccup <- gr$score

  #--------------------
  # Filter out the parts where bin_weight == 0
  mcols(gr_SignalTrack) <- mcols(gr_SignalTrack)[, -1]
  return(gr_SignalTrack[gr_SignalTrack$bin_weight > 0])
}

merge_with_tumor <- function(gr_tumor, gr_SignalTrack) {
  overlaps <- findOverlaps(gr_tumor, gr_SignalTrack)
  Xmat <- matrix(0, nrow = length(gr_tumor), ncol = ncol(mcols(gr_SignalTrack)))
  Xmat[queryHits(overlaps), ] <- as.matrix(mcols(gr_SignalTrack[subjectHits(overlaps)]))
  colnames(Xmat) <- colnames( mcols(gr_SignalTrack))
  mcols(gr_tumor) <- cbind(mcols(gr_tumor), Xmat)
  gr_tumor
}

merge_CopyTrack_tumor <- function(gr_tumor, gr_CopyTrack){
  samples <- unique(gr_tumor$sample)
  copy_number <- rep(NA, length(gr_tumor))
  names(copy_number) <- gr_tumor$sample
  for(s in seq_along(samples)) {
    gr_tumor_sample <- gr_tumor[gr_tumor$sample == samples[s]]
    gr_CopyNum <- gr_CopyTrack[subjectHits(findOverlaps(gr_tumor_sample, gr_CopyTrack))]
    copy_number[names(copy_number) == samples[s]] <-  mcols(gr_CopyNum)[,colnames(mcols(gr_CopyNum)) == samples[s]]
  }
  gr_tumor$copy_number <- copy_number
  return(gr_tumor)
}

standardize_covariates <- function(x, q_min = 0.001, q_max = 0.999){
  y <- pmin(x, quantile(x, q_max))
  y <- pmax(y, quantile(x, q_min))
  scale(y)
}

normalize_covariates <- function(x, q_min = 0.001, q_max = 0.999){
  y <- pmin(x, quantile(x, q_max))
  y <- pmax(y, quantile(x, q_min))
  (y - min(y))/(max(y) - min(y))
}

get_Precoditioning_matrices <- function(out, SignalTrack, bin_weight) {
  p <- nrow(out$Betas)
  K <- ncol(out$Betas)
  PrecondSigma <- array(NA, dim = c(p, p, K))
  theta_sum <- rowSums(out$Thetas)
  W <- exp(SignalTrack %*% out$Betas)
  for(k in 1:K){
    H <- theta_sum[k] * crossprod(SignalTrack * bin_weight * W[, k],
                                  SignalTrack) + 1/out$Mu[k]
    PrecondSigma[, , k] <- solve((H + t(H))/2)
  }
  return(PrecondSigma)
}

get_PosteriorEffectiveSize <- function(chain){
  apply(chain, c(2,3), function(x) coda::effectiveSize(x))
}

generate_SigPoisProcess_init <- function(data,
                                         K = 12,
                                         prior_params = SigPoisProcess_mult.PriorParams(),
                                         nnls = FALSE,
                                         sigStartNNLS = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8",
                                                          "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")){
  # Useful quantities
  J <- ncol(data$CopyTrack)
  I <- length(unique(data$gr_Mutations$channel))
  p <- ncol(data$SignalTrack)
  prior_params <- do.call("SigPoisProcess_mult.PriorParams", prior_params)
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

#----- Load the tumor data
reprocess_data <- FALSE
if(reprocess_data) {
  gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Breast-AdenoCa_snp.rds.gzip")
  # Remove mutations from blacklist
  blacklist <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/hg19-blacklist.v2.bed")
  overlaps <- findOverlaps(gr_tumor, blacklist)
  gr_tumor <- gr_tumor[-queryHits(overlaps)]
  # Remove chrY (we are dealing with breast cancer)
  gr_tumor <- gr_tumor[seqnames(gr_tumor) != "chrY"]

  #----- Build the CopyTrack to record the copy number alteration
  # Load the copy number data
  df_copy <- read_tsv("~/20170119_final_consensus_copynumber_donor")
  gr_copy <- GRanges(
    seqnames = paste0("chr",df_copy$chr),
    ranges = IRanges(start = df_copy$start, end = df_copy$end),
    strand = "*",
    sample = df_copy$sampleID,
    score = df_copy$value)

  # Build the track
  gr_CopyTrack <- build_CopyTrack(gr_tumor, gr_copy, tilewidth = 2000)

  # Merge CopyTrack with Tumor
  #gr_tumor <- merge_CopyTrack_tumor(gr_tumor, gr_CopyTrack)

  #----- Build the SignalTrack and standardize the covariates (after capping at 0.999 quantile)
  gr_SignalTrack_2kb <- build_SignalTrack(tilewidth = 2000, type = "avg")
  gr_SignalTrack_2kb_std <- gr_SignalTrack_2kb
  mcols(gr_SignalTrack_2kb_std)[, -1] <- apply(as.data.frame(mcols(gr_SignalTrack_2kb_std)[, -1]), 2,
                                               function(x) standardize_covariates(x))
  # Merge SignalTrack with tumor
  gr_tumor_2kb_std <- merge_with_tumor(gr_tumor, gr_SignalTrack_2kb_std)

  # Mutations data
  gr_Mutations <- gr_tumor_2kb_std
  gr_Mutations$bin_weight <- NULL

  # Signatures to use
  Cosmic_Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
  sigs_to_use <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8", "SBS20",
                   "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")
  Sigs <- Cosmic_Sigs[, sigs_to_use]

  # Baseline model
  mutMatrix <- getTotalMutations(gr_Mutations)
  mutMatrix <- getTotalMutations(gr_Mutations)
  outBase <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                                K = 20,
                                                alpha = 1.01, a = 1.01)
  CompressiveNMF::plot_SBS_signature(outBase$Signatures)
  # SignalTrack and CopyTrack matrices
  SignalTrack <- as.matrix(mcols(gr_SignalTrack_2kb_std))[, -1]
  CopyTrack <- apply(as.matrix(mcols(gr_CopyTrack))[, -1], 2,
                     function(x) x * gr_CopyTrack$bin_weight)[, colnames(mutMatrix)]

  # Save all data now to avoid re-loading them
  # saveRDS(list("gr_Mutations" = gr_Mutations,
  #              "SignalTrack" = SignalTrack,
  #              "CopyTrack" = CopyTrack),
  #         file = "~/SigPoisProcess/data/ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip")

}


################################################################
# Calculate the MAP from 10 random starting points
################################################################
#--- Load the merged data
file_data <- "~/SigPoisProcess/data/ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip"
data <- readRDS(file_data)

# Control parameters for the MAP
maxiter <- 5000
controls <- SigPoisProcess_mult.controls(maxiter = maxiter, tol = 1e-7)
# Number of signatures
K <- 12

n_start_point <- 10
# Generate the random starting point
regenerate_start_points <- FALSE
if(regenerate_start_points) {
  set.seed(10)
  for(i in 1:n_start_point){
    out_dir <- paste0("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate", sprintf("%02d", i))
    create_directory(out_dir)
    init <- generate_SigPoisProcess_init(data = data, K = K)
    saveRDS(init, paste0(out_dir, "/init.rds"))
  }
  # Generate also the starting point from known signatues and nnls
  out_dir <- paste0("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate", "NNLS")
  create_directory(out_dir)
  init <- generate_SigPoisProcess_init(data = data, nnls = TRUE,
                                       sigStartNNLS = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8",
                                                        "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30"))
  saveRDS(init, paste0(out_dir, "/init.rds"))

}

#---- RUN the MAP now from the chosen starting point
nrepl <- 3 #"NNLS" #3 # <------ Change the value based on the replicate
print(nrepl)
if(is.numeric(nrepl)){
  out_dir <- paste0("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate", sprintf("%02d", nrepl), "/")
} else if (nrepl == "NNLS"){
  out_dir <- paste0("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate", nrepl , "/")
}

rerun_MAP <- FALSE
if(rerun_MAP) {
  # Load initial value
  init <- readRDS(paste0(out_dir, "init.rds"))
  # Find the MAP
  print(paste0("Replicate ", nrepl))
  start_time <- Sys.time()
  outMAP <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
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

#---- RUN the MCMC now from the MAP solution
rerun_MCMC <- FALSE
if(rerun_MCMC){
  # Load the MAP
  outMAP <- readRDS(paste0(out_dir, "MAPSolution.rds.gzip"))
  #------ MCMC
  nsamples_all <- 5000
  nsamples_current <- 0
  nsamples <- 100
  # Run the model storing the output every 100 iterations
  file_name <- paste0(out_dir, "MCMCSolution.rds.gzip")
  set.seed(10)
  while(nsamples_current < nsamples_all) {
    print(paste0("<---- N. samples: ", nsamples_current, " ----->"))

    # Initialize either from the MAP or from the last previous iteration
    if(!file.exists(file_name)) {
      init <- SigPoisProcess_mult.init(R_start = outMAP$Signatures,
                                       Theta_start = outMAP$Thetas,
                                       Betas_start = outMAP$Betas,
                                       Mu_start = outMAP$Mu,
                                       Sigma2_start = outMAP$Sigma2)
      controls <- SigPoisProcess_mult.controls(nsamples = nsamples, burnin = 10)
    } else {
      # Load MCMC chain
      outMCMC <- readRDS(file_name)
      id <- dim(outMCMC$chain$SIGSchain)[1]
      init <- SigPoisProcess_mult.init(R_start = outMCMC$chain$SIGSchain[id, , ],
                                       Theta_start = outMCMC$chain$THETAchain[id, , ],
                                       Betas_start = outMCMC$chain$BETASchain[id, , ],
                                       Mu_start = outMCMC$chain$MUchain[id, ],
                                       Sigma2_start = outMCMC$chain$SIGMA2chain[id, ])
      controls <- SigPoisProcess_mult.controls(nsamples = nsamples, burnin = 10)
    }

    # Run the model for nsamples iterations
    out_temp <-  SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
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
# Now, run the MAP and the MCMC for a fixed set of signatures
################################################################

Sigs_to_use <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8", "SBS20",
                 "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")
CosmicSigs <- SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37[, Sigs_to_use]

# Control parameters for the MAP
maxiter <- 5000
controls <- SigPoisProcess_mult.controls(maxiter = maxiter, tol = 1e-6)
# Number of signatures
K <- length(Sigs_to_use)

out_dir <- "~/SigPoisProcess/output/Application_Breast/SemiSupervised/"
rerun_MAP_fixedSigs <- FALSE
if(rerun_MAP_fixedSigs) {
  set.seed(10)
  start_time <- Sys.time()
  outMAP <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
                                  SignalTrack = data$SignalTrack,
                                  CopyTrack = data$CopyTrack,
                                  K = K,
                                  verbose = TRUE,
                                  update_Signatures = FALSE,
                                  init = SigPoisProcess_mult.init(R_start = CosmicSigs),
                                  controls = controls)
  end_time <- Sys.time()
  outMAP$time <- end_time - start_time

  # Save it
  saveRDS(outMAP,
          file = paste0(out_dir, "MAPSolution_SemiSupervised_v2.rds.gzip"),
          compress = "gzip")
}


rerun_MCMC_fixedSigs <- TRUE
force_zero <- FALSE # <--- force Betas to zero in signatures that are compressed, for the MCMC starting point
start_from_map <- FALSE
if(rerun_MCMC_fixedSigs) {
  # Load the MAP
  outMAP <- readRDS(paste0(out_dir, "MAPSolution_SemiSupervised_v2.rds.gzip"))
  #------ MCMC
  nsamples_all <- 10000
  nsamples_current <- 0
  nsamples <- 100
  # Run the model storing the output every 100 iterations
  if(force_zero) {
    file_name <- paste0(out_dir, "MCMCSolution_SemiSupervised_ForcedZero.rds.gzip")
    print("forced")
  } else if (start_from_map) {
    file_name <- paste0(out_dir, "MCMCSolution_SemiSupervised_v2.rds.gzip")
  } else {
    file_name <- paste0(out_dir, "MCMCSolution_SemiSupervised_noMap.rds.gzip")
    print("Start from random")
  }
  set.seed(10)
  while(nsamples_current < nsamples_all) {
    print(paste0("<---- N. samples: ", nsamples_current, " ----->"))
    # Initialize either from the MAP or from the last previous iteration
    if(!file.exists(file_name)) {
      if(force_zero){
        Betas_forced <- outMAP$Betas
        Betas_forced[, outMAP$Mu < 0.002] <- 0
        init <- SigPoisProcess_mult.init(R_start = outMAP$Signatures,
                                         Theta_start = outMAP$Thetas,
                                         Betas_start = Betas_forced,
                                         Mu_start = outMAP$Mu,
                                         Sigma2_start = outMAP$Sigma2)
      } else if (start_from_map) {
        init <- SigPoisProcess_mult.init(R_start = outMAP$Signatures,
                                         Theta_start = outMAP$Thetas,
                                         Betas_start = outMAP$Betas,
                                         Mu_start = outMAP$Mu,
                                         Sigma2_start = outMAP$Sigma2)
      } else {
        init <- SigPoisProcess_mult.init(R_start = CosmicSigs)
      }

      controls <- SigPoisProcess_mult.controls(nsamples = nsamples, burnin = 10)
    } else {
      # Load MCMC chain
      outMCMC <- readRDS(file_name)
      id <- dim(outMCMC$chain$SIGSchain)[1]
      init <- SigPoisProcess_mult.init(R_start = outMCMC$chain$SIGSchain[id, , ],
                                       Theta_start = outMCMC$chain$THETAchain[id, , ],
                                       Betas_start = outMCMC$chain$BETASchain[id, , ],
                                       Mu_start = outMCMC$chain$MUchain[id, ],
                                       Sigma2_start = outMCMC$chain$SIGMA2chain[id, ])
      controls <- SigPoisProcess_mult.controls(nsamples = nsamples, burnin = 10)
    }

    # Run the model for nsamples iterations
    out_temp <-  SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
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

outMCMCFixed <- readRDS("~/SigPoisProcess/output/Application_Breast/SemiSupervised/MCMCSolution_SemiSupervised.rds.gzip")
outMCMCFixed2 <- readRDS("~/SigPoisProcess/output/Application_Breast/SemiSupervised/MCMCSolution_SemiSupervised_noMap.rds.gzip")
plot(outMCMCFixed$chain$logPostchain, type = "l", ylim = c(-12366500, -12361600))
lines(outMCMCFixed2$chain$logPostchain, type = "l", col = "red")

# outMCMCFixed <- readRDS("~/SigPoisProcess/output/Application_Breast/SemiSupervised/MCMCSolution_SemiSupervised.rds.gzip")
#
# CompressiveNMF::plot_SBS_signature(outMAP$Signatures) + plot_betas2(outMAP$Betas)
for(k in colnames(outMCMCFixed2$Signatures)){
  print(paste0(k, " ", coda::effectiveSize(outMCMCFixed2$chain$MUchain[-c(1:2000), k])))
  print(Sigs_to_use[which(colnames(outMCMCFixed2$Signatures) == k)])
  plot(outMCMCFixed2$chain$MUchain[-c(1:2000), k], type = "l")
  #plot(outMCMC$chain$MUchain[, k], type = "l")
  # abline(v = 6200, col = "red")
  readline(prompt="Press [enter] to continue")
}
#
for(k in colnames(outMCMCFixed$Signatures)){
  for(b in rownames(outMCMCFixed$Betas)){
    print(paste0(k, " ", b, " ", coda::effectiveSize(outMCMCFixed$chain$BETASchain[-c(1:2000), b, k])))
    print(Sigs_to_use[which(colnames(outMCMCFixed$Signatures) == k)])
    plot(outMCMCFixed$chain$BETASchain[-c(1:2000), b, k], type = "l")
    #abline(v = 5000)
    readline(prompt="Press [enter] to continue")
  }
}
plot(outMCMCFixed2$chain$BETASchain[, 9, "SigN06"], type = "l")

#---- make some cool plots
plot_vector_facets <- function(df_Assign) {
  require(ggplot2)
  require(dplyr)

  df_Assign <- df_Assign %>%
    dplyr::mutate(
      best_sig = factor(best_sig, levels = best_sig),
      compressed = ifelse(m == 0 | mu < 0.05, "compressed", "not compressed")
    )

  gg <- ggplot(df_Assign, aes(x = 1, y = 1, size = mu, fill = m, shape = compressed)) +
    geom_point(color = "black", stroke = 0.7) +
    facet_wrap(~ best_sig, ncol = 1, strip.position = "left") +
    scale_size(name = expression(mu[k]), range = c(3, 12)) +
    scale_fill_gradientn(
      name = "N. mutations",
      colours = c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                  "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")
    ) +
    scale_shape_manual(
      name = "Compressed",
      values = c("compressed" = 4, "not compressed" = 21)
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      strip.background = element_rect(fill = 'white', color = NA),
      legend.position = "left",
      panel.spacing = grid::unit(0.03, "lines")
    )

  return(gg)
}

#----- fixed case 1
SigsCisFix <- get_PosteriorCI(outMCMCFixed$chain$SIGSchain[-c(1:2000), , ])
ThetaCisFix <- get_PosteriorCI(outMCMCFixed$chain$THETAchain[-c(1:2000), , ])
BetasCisFix <- get_PosteriorCI(outMCMCFixed$chain$BETASchain[-c(1:2000), , ])
MuFix <- colMeans(outMCMCFixed$chain$MUchain[-c(1:2000), ])
names(MuFix)  <- colnames(CosmicSigs)
colnames(BetasCisFix$mean) <- colnames(BetasCisFix$lowCI) <- colnames(BetasCisFix$highCI) <- colnames(CosmicSigs)
p1 <- plot_betas_capped(BetasCisFix$mean, BetasCisFix$lowCI, BetasCisFix$highCI) + theme(legend.position = "none")

# Fraction of mutations attributed to the signature
bigProd <- SigsCisFix$mean[data$gr_Mutations$channel, ] * t(ThetaCisFix$mean[, data$gr_Mutations$sample]) *
  exp(as.matrix(mcols(data$gr_Mutations)[, -c(1:3)]) %*% BetasCisFix$mean)
mm <- apply(bigProd, 1, which.max)
df_AssignFix <- data.frame(data$gr_Mutations) %>%
  mutate(best_sig = colnames(CosmicSigs)[mm],
         #best_sig = paste0("SigN", sprintf("%02d", best_sig)),
         best_sig = factor(best_sig, sort(colnames(CosmicSigs)))) %>%
  group_by(sample, best_sig) %>%
  summarize(m = n(), .groups = "drop") %>%
  tidyr::complete(sample, best_sig, fill = list(m = 0)) %>%
  ungroup() %>%
  group_by(best_sig) %>%
  summarize(m = sum(m)) %>%
  left_join(data.frame(best_sig = names(MuFix),
                       mu = MuFix), by = "best_sig")

p_mu1 <- plot_vector_facets(df_AssignFix)
p_mu1 + p1 + patchwork::plot_layout(widths = c(0.5,2))

#----- fixed case 2
SigsCisFix2 <- get_PosteriorCI(outMCMCFixed2$chain$SIGSchain[-c(1:2000), , ])
ThetaCisFix2 <- get_PosteriorCI(outMCMCFixed2$chain$THETAchain[-c(1:2000), , ])
BetasCisFix2 <- get_PosteriorCI(outMCMCFixed2$chain$BETASchain[-c(1:2000), , ])
MuFix2 <- colMeans(outMCMCFixed2$chain$MUchain[-c(1:2000), ])
names(MuFix2)  <- colnames(CosmicSigs)
colnames(BetasCisFix2$mean) <- colnames(BetasCisFix2$lowCI) <- colnames(BetasCisFix2$highCI) <- colnames(CosmicSigs)
p2 <- plot_betas_capped(BetasCisFix2$mean, BetasCisFix2$lowCI, BetasCisFix2$highCI) + theme(legend.position = "none")

# Fraction of mutations attributed to the signature
bigProd2 <- SigsCisFix2$mean[data$gr_Mutations$channel, ] * t(ThetaCisFix2$mean[, data$gr_Mutations$sample]) *
  exp(as.matrix(mcols(data$gr_Mutations)[, -c(1:3)]) %*% BetasCisFix2$mean)
mm2 <- apply(bigProd2, 1, which.max)
df_AssignFix2 <- data.frame(data$gr_Mutations) %>%
  mutate(best_sig = colnames(CosmicSigs)[mm2],
         #best_sig = paste0("SigN", sprintf("%02d", best_sig)),
         best_sig = factor(best_sig, sort(colnames(CosmicSigs)))) %>%
  group_by(sample, best_sig) %>%
  summarize(m = n(), .groups = "drop") %>%
  tidyr::complete(sample, best_sig, fill = list(m = 0)) %>%
  ungroup() %>%
  group_by(best_sig) %>%
  summarize(m = sum(m)) %>%
  left_join(data.frame(best_sig = names(MuFix2),
                       mu = MuFix2), by = "best_sig")
p_mu2 <- plot_vector_facets(df_AssignFix2)
p_mu1 + p1 + p_mu2 + p2 + patchwork::plot_layout(widths = c(0.3, 2, 0.3,2))



p1 + p2

mu_tile <- colMeans(outMCMCFixed$chain$MUchain[-c(1:1000), ])


CompressiveNMF::plot_SBS_signature(SigsCisFix$mean)

plot_betas_capped <- function(Betas_sol, lowCI, highCI) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)

  df_betas <- as.data.frame(Betas_sol) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "Beta", -Covariate)

  df_low <- as.data.frame(lowCI) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "lowCI", -Covariate)

  df_high <- as.data.frame(highCI) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "highCI", -Covariate)

  df_betas <- df_betas %>%
    left_join(df_low, by = c("Covariate", "Signature")) %>%
    left_join(df_high, by = c("Covariate", "Signature")) %>%
    mutate(
      Covariate = factor(Covariate, levels = unique(Covariate)),
      Signature = as.factor(Signature),
      Beta = ifelse(lowCI < 0 & highCI > 0, NA, Beta)
    ) %>%
    # Add capped value for colors only
    mutate(
      Beta_capped = case_when(
        is.na(Beta) ~ NA_real_,
        Beta < -1 ~ -1,
        Beta >  1 ~  1,
        TRUE      ~ Beta
      )
    )

  oldopt <- options()
  options(ggplot2.continuous.colour="ignore", ggplot2.continuous.fill="ignore")

  plot_out <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta_capped)) +
    geom_tile(color = "black", width = 1, height = 1, na.rm = TRUE) +
    scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"),
      values = scales::rescale(c(-1, 0, 1)),
      limits = c(-1, 1),
      na.value = "gray90"
    ) +
    geom_text(aes(label = ifelse(is.na(Beta), "", round(Beta, 2))),
              size = 3, na.rm = TRUE) +
    scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1, size = 10, colour = "black"),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks = element_blank(),
      panel.grid = element_blank())

  options(oldopt)
  return(plot_out)
}
plot_betas_capped_size <- function(Betas_sol, lowCI, highCI, sig_vec) {
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)

  # Prep data
  df_betas <- as.data.frame(Betas_sol) %>%
    rownames_to_column("Covariate") %>%
    tidyr::gather("Signature", "Beta", -Covariate)

  df_low <- as.data.frame(lowCI) %>%
    rownames_to_column("Covariate") %>%
    tidyr::gather("Signature", "lowCI", -Covariate)

  df_high <- as.data.frame(highCI) %>%
    rownames_to_column("Covariate") %>%
    tidyr::gather("Signature", "highCI", -Covariate)

  df_betas <- df_betas %>%
    left_join(df_low, by = c("Covariate", "Signature")) %>%
    left_join(df_high, by = c("Covariate", "Signature")) %>%
    mutate(
      Covariate = factor(Covariate, levels = unique(Covariate)),
      Signature = as.factor(Signature),
      Beta = ifelse(lowCI < 0 & highCI > 0, NA, Beta)
    ) %>%
    mutate(
      Beta_capped = case_when(
        is.na(Beta) ~ NA_real_,
        Beta < -1 ~ -1,
        Beta > 1 ~ 1,
        TRUE ~ Beta
      )
    ) %>%
    # Join tile size (match sig_vec by Signature)
    left_join(
      tibble(Signature = names(sig_vec), tile_size_raw = as.numeric(sig_vec)),
      by = "Signature"
    ) %>%
    # Rescale tile_size to [0.2, 1] for width/height
    mutate(
      tile_size = rescale(tile_size_raw, to = c(0.2, 1))
    )

  options(ggplot2.continuous.colour="ignore", ggplot2.continuous.fill="ignore")

  p <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta_capped)) +
    geom_tile(aes(width = tile_size, height = tile_size),
              color = "black", na.rm = TRUE) +
    scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"),
      values = scales::rescale(c(-1, 0, 1)),
      limits = c(-1, 1),
      na.value = "gray90"
    ) +
    geom_text(aes(label = ifelse(is.na(Beta), "", round(Beta, 2))),
              size = 3, na.rm = TRUE) +
    scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1, size = 10, colour = "black"),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks = element_blank(),
      panel.grid = element_blank())

  return(p)
}
plot_betas_capped_size(BetasCisFix$mean, BetasCisFix$lowCI, BetasCisFix$highCI, mu_tile)


if(FALSE){

controls <- SigPoisProcess_mult.controls(maxiter = 200, tol = 1e-7)
grRepli <- data$gr_Mutations
mcols(grRepli) <- mcols(grRepli)[, c("tumor", "sample", "channel", "RepliTime")]
sigRepli <- matrix(data$SignalTrack[, "RepliTime"])
colnames(sigRepli) <- "RepliTime"
outMAP_repli <- SigPoisProcess_multCN(gr_Mutations = grRepli,
                                SignalTrack = sigRepli,
                                CopyTrack = data$CopyTrack,
                                K = ncol(CosmicSigs),
                                verbose = TRUE,
                                update_Signatures = FALSE,
                                init = SigPoisProcess_mult.init(R_start = CosmicSigs),
                                controls = controls)


cbind(outMAP_repli$Mu, c(outMAP_repli$Betas))
#

# sort(c("nnls" = tail(outMAP_nnls$sol$trace, 1),
#   "map1" = tail(outMAP1$sol$trace, 1),
#   "map2" = tail(outMAP2$sol$trace, 1),
#   "map3" = tail(outMAP3$sol$trace, 1)))
#
# outMAP_nnls <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/ReplicateNNLS/MAPSolution.rds.gzip")
# outMAP1 <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate01/MAPSolution.rds.gzip")
# outMAP2 <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate02/MAPSolution.rds.gzip")
# outMAP3 <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate03/MAPSolution.rds.gzip")
#
# pnnls <- CompressiveNMF::plot_SBS_signature(outMAP_nnls$Signatures) + plot_betas2(outMAP_nnls$Betas)
#CompressiveNMF::plot_SBS_signature(outMAP1$Signatures) + plot_betas2(outMAP1$Betas)
# p1 <- CompressiveNMF::plot_SBS_signature(outMAP2$Signatures) + plot_betas2(outMAP2$Betas)
# p2 <- CompressiveNMF::plot_SBS_signature(outMAP3$Signatures) + plot_betas2(outMAP3$Betas)
# round(outMAP$Mu, 4)
# p1
# p2
# match_to_RefSigs(outMAP1$Signatures, CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
# SigPoisProcess::match_to_RefSigs(outMAP3$Signatures, CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
#
# round(outMAP1$Mu,  4)
# outMCMC <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/ReplicateNNLS/MCMCSolution.rds.gzip")
# library(patchwork)
# outMAP <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate03/MAPSolution.rds.gzip")
# outMCMC <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate03/MCMCSolution.rds.gzip")
R <- apply(outMCMC$chain$SIGSchain[-c(1:7000), ,], c(2,3), mean)
Theta <- apply(outMCMC$chain$THETAchain[-c(1:7000), ,], c(2,3), mean)
Betas <- apply(outMCMC$chain$BETASchain[-c(1:7000), ,], c(2,3), mean)

CIs <- get_PosteriorCI(outMCMC$chain$SIGSchain[-c(1:7000), ,])
BetasCIs <- get_PosteriorCI(outMCMC$chain$BETASchain[-c(1:7000), ,])
CompressiveNMF::plot_SBS_signature(CIs$mean, lowCI = CIs$lowCI, highCI = CIs$highCI) + plot_betas2(Betas)
CompressiveNMF::plot_SBS_signature(R) + CompressiveNMF::plot_SBS_signature(outMAP_nnls$Signatures)
SigPoisProcess::match_to_RefSigs(R[, mu > 0.01], CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
#plot(outMCMC$chain$MUchain[, 8], type = "l")
eff <- get_PosteriorEffectiveSize(outMCMC$chain$BETASchain[-c(1:8000), ,])
problems <- which(eff<10, arr.ind = TRUE)
for(i in 1:nrow(problems)){
  print(paste0(colnames(Betas)[problems[i, 2]], " ",
               rownames(Betas)[problems[i, 1]], " ",
               coda::effectiveSize(outMCMC$chain$BETASchain[-c(1:8000), problems[i, 1], problems[i, 2]])))
  plot(outMCMC$chain$BETASchain[-c(1:8000), problems[i, 1], problems[i, 2]], type = "l")
  readline(prompt="Press [enter] to continue")
}

mu <- colMeans(outMCMC$chain$MUchain[-c(1:7000), ])
for(k in colnames(R[, mu>0.01])){
  for(b in rownames(Betas[, mu>0.01])){
    print(paste0(k, " ", b, " ", coda::effectiveSize(outMCMC$chain$BETASchain[-c(1:7000), b, k])))
    plot(outMCMC$chain$BETASchain[-c(1:7000), b, k], type = "l")
    #abline(v = 5000)
    readline(prompt="Press [enter] to continue")
  }
}

eff <- get_PosteriorEffectiveSize(outMCMC$chain$BETASchain[-c(1:7000), ,mu>0.01])
plot(outMCMC$chain$BETASchain[, "RepliTime", "SigN10"], type = "l")
for(k in colnames(R)){
    print(paste0(k, " ", coda::effectiveSize(outMCMC$chain$MUchain[-c(1:7000), k])))
    plot(outMCMC$chain$MUchain[-c(1:7000), k], type = "l")
    #plot(outMCMC$chain$MUchain[, k], type = "l")
   # abline(v = 6200, col = "red")
    readline(prompt="Press [enter] to continue")

}

outSim <- readRDS("../output/Simulation/Scenario_A_indep/Simulation_10/output_mcmc_FullModel.rds.gzip")
plot(outSim$chain$logPostchain)
plot(outMCMC$chain$logPostchain[-c(1:8000)], type = "l")
CompressiveNMF::plot_SBS_signature(outSim$init$R_start) + CompressiveNMF::plot_SBS_signature(outSim$Signatures)
# plot(outMCMC$chain$logPostchain, type = "l")

l1 <- outMCMC$chain$logPostchain
l2 <- outMCMC$chain$logPostchain
plot(l1[-c(1:6000)], type = "l")
lines(l2, col = "red")
# CompressiveNMF::plot_SBS_signature(R) + plot_betas2(Betas)
# CompressiveNMF::plot_SBS_signature(outMAP$Signatures) + plot_betas2(outMAP$Betas)
# CompressiveNMF::plot_SBS_signature(R) + CompressiveNMF::plot_SBS_signature(outMAP$Signatures)
#
# bigProd <- R[data$gr_Mutations$channel, ] * t(Theta[, data$gr_Mutations$sample]) * exp(as.matrix(mcols(data$gr_Mutations)[, -c(1:3)]) %*% Betas)
#
# bigProbs <- t(apply(bigProd, 1, function(x) x/sum(x)))
#
# mm <- apply(bigProd, 1, which.max)
# mutAssign <- rep(0, 12)
# mutAssign[as.numeric(names(table(mm)))] <- c(table(mm))
# cbind(mutAssign, round(colMeans(outMCMC$chain$MUchain[-c(1:500), ]), 3))
mu <- colMeans(outMCMC$chain$MUchain[-c(1:8000), ])
CompressiveNMF::plot_SBS_signature(R[, mu > 0.01]) + plot_betas2(Betas[, mu > 0.01])
SigPoisProcess::match_to_RefSigs(R[, mu > 0.01], CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
devtools::document()
plot_96SBSsig(R[, mu > 0.01]) + plot_betas2(Betas[, mu > 0.01])






# PCAtrack <- princomp(data$SignalTrack)
# cumsum(PCAtrack$sdev^2/sum(PCAtrack$sdev^2))
# PCAtrack$loadings
# plot(PCAtrack$scores[1:1000, 1], type = "l")
# lines(PCAtrack$scores[1:1000, 2], type = "l", col = "red")


#--------- Reproduce the figure for the Research statement
#--- Load the merged data
file_data <- "~/SigPoisProcess/data/ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip"
data <- readRDS(file_data)

# Control parameters for the MAP
maxiter <- 5000
controls <- SigPoisProcess_mult.controls(maxiter = maxiter, tol = 1e-7)
# Number of signatures
K <- 12

Cosmic_Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
sigs_to_use <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8",
                 "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")
Sigs <- Cosmic_Sigs[, sigs_to_use]
mutMatrix <- getTotalMutations(data$gr_Mutations)
BaseNMF <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                              K = 0,
                                              S = 1e6 * Sigs + 1,
                                              alpha = 1.01, a = 1.01)



set.seed(10)
outFixed <- SigPoisProcess_multCN(gr_Mutations = data$gr_Mutations,
                                  SignalTrack = data$SignalTrack,
                                  CopyTrack = data$CopyTrack,
                                  K = ncol(Sigs),
                                  method = "map",
                                  verbose = TRUE,
                                  init = SigPoisProcess_mult.init(R_start = Sigs),
                                  controls = SigPoisProcess_mult.controls(maxiter = 200, tol = 1e-6),
                                  update_Signatures = FALSE)
saveRDS(outFixed, "~/SigPoisProcess/output/Application_Breast/Plot_for_Statement.rds")
CompressiveNMF::plot_SBS_signature(outFixed$Signatures) + plot_betas(outFixed$Betas)

gr_SignalTrack <-gr_SignalTrack_2kb_std
MultNMF <- outMCMC

chr <- "chrX"
region_start <- 1e6
region_end <- 5e6

# Make the plot now.
calculate_ChannelProbs_region <- function(chr,
                                          region_start, region_end,
                                          j = 1,
                                          MultNMF,
                                          gr_SignalTrack,
                                          CopyTrack) {
  # Extract solution
  Theta <- MultNMF$Thetas
  Betas <- MultNMF$Betas
  R <- MultNMF$Signatures

  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, gr_SignalTrack)
  gr_subset <- gr_SignalTrack[subjectHits(overlaps)]
  # Estimate the loadings in each subset
  SignalTrack_subset <- as.matrix(mcols(gr_subset)[, -1])
  CopyTrack_subset <- CopyTrack[subjectHits(overlaps), j]
  ExpSolBetas <- exp(SignalTrack_subset %*% Betas)
  # Big Storage
  ChannelProbs <- array(NA, dim = c(nrow(ExpSolBetas), ncol(ExpSolBetas), 96),
                        dimnames = list(NULL, colnames(ExpSolBetas), rownames(R)))

  ExpSolBetas_copy <- apply(ExpSolBetas, 2, function(x) x * CopyTrack_subset)
  for(i in 1:96){
    SigTheta <- t(R[i, ] * Theta[, j])
    bigProd <- SigTheta[rep(1, nrow(ExpSolBetas_copy)), ] * ExpSolBetas_copy
    sig_probs <-  bigProd#/rowSums(bigProd)
    ChannelProbs[, , i] <- sig_probs
  }
  ChannelProbs
}

get_signal_regions <- function(chr, region_start, region_end, gr_SignalTrack) {
  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, gr_SignalTrack)
  gr_subset <- gr_SignalTrack[subjectHits(overlaps)]
  return(gr_subset)
}

j <- "DO220823"#"DO218709" #"DO1006"
chr <- "chr1"
region_start <- 10e6#15e6
region_end   <- 17e6#22e6
ChannelProbs <- calculate_ChannelProbs_region(chr,
                                              region_start,
                                              region_end,
                                              j = j,
                                              outFixed,
                                              gr_SignalTrack_2kb_std, data$CopyTrack)

i <- "A[C>A]C"

dimnames(BaseNMF$mapOutput$R) <- dimnames(Sigs)
colnames(BaseNMF$mapOutput$Theta) <- colnames(mutMatrix)
rownames(BaseNMF$mapOutput$Theta) <- colnames(Sigs)
pr_const <- 1800 * BaseNMF$mapOutput$R[i, ] * BaseNMF$mapOutput$Theta[, j]/(sum(gr_SignalTrack$bin_weight))
df_const <- data.frame(t(pr_const)[rep(1,nrow(ChannelProbs)), ]) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start,
         method = "Standard NMF")

df_mod <- data.frame(ChannelProbs[, , i])
colnames(df_mod) <- colnames(BaseNMF$mapOutput$R)


window <- 10
p_intensity <- df_mod %>%
  mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
                                                   align = "center"))) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start,
         method = "PoissonProcess") %>%
  bind_rows(df_const) %>%
  gather(key = "Sig", value = "Prob", - region, -method) %>%
  filter(Sig %in% c("SBS1", "SBS3", "SBS5", "SBS8")) %>%#head(names(sort(pr_const, decreasing = TRUE)), 4))%>%
  ggplot() +
  geom_line(aes(x = region, y = Prob, color = Sig, linetype = method), linewidth = 0.7) +
  scale_color_manual(name = "Signature",
                     values = c("gray50", "red", "darkorange", "blue")) +
  ylab("Mutation intensity") +
  #xlab(paste0("Genomic regions in ",  chr, "-",
  #            region_start,":", region_end, " (2kb)")) +
  xlab(paste0("Genomic regions in ",  chr, " (2kb)")) +
  facet_wrap(~paste0("Mutation type ", i, " - Patient ", j)) +
  theme_minimal(base_size = 10) +
  scale_x_continuous(expand = c(0,0))

# Now, make the plot of the various histone marks below it.
subset_tracks <- get_signal_regions(chr, region_start, region_end, gr_SignalTrack_2kb_std)
subset_tracks$bin_weight <- NULL
df_tracks <- as.data.frame(mcols(subset_tracks)) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start) %>%
  gather(key = "covariate", value = "score", -region)

p_track <- ggplot(df_tracks) +
  geom_raster(aes(x = region + 1000, y = covariate, fill = score)) +
  #facet_wrap(.~ pos2, nrow = 1) +
  #scale_fill_gradientn(name = "Signal track\n(standardized)",
  #                     colours =c("antiquewhite", "#F6C866", "#F2AB67", "#EF8F6B",
  #                                "#ED7470", "#BF6E97", "#926AC2",
  #                                "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
  scale_fill_gradient2(name = "Signal track\n(standardized)",
                       mid = "white", low = "blue", high = "red", midpoint = 0,
                       limits = c(-3.5, 3.5), oob = scales::squish)+
  theme_minimal(base_size = 10) +
  scale_x_continuous(expand = c(0,0)) +
  #theme(axis.title = element_blank())
  ylab("Genomic covariate")+
  xlab(paste0("Genomic regions in ",  chr, " (2kb)"))

p_intensity/p_track + plot_layout(heights = c(2, 1.5))

# Now plot just the mutational signatures
p_sigs <- CompressiveNMF::plot_SBS_signature(Sigs[, c("SBS1", "SBS3", "SBS5", "SBS8")])+
  theme(axis.text.x = ggplot2::element_blank())
colnames(outFixed$Betas) <- colnames(Sigs)
p_betas <- plot_betas2(outFixed$Betas[, c("SBS1", "SBS3", "SBS5", "SBS8")])+
  theme(axis.text.y = ggplot2::element_blank())
p_sigs + p_betas+ plot_layout(widths = c(2, 2))


}







# plot_betas3 <- function(Betas_sol, lowCI, highCI) {
#
#   df_betas <- as.data.frame(Betas_sol) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "Beta", -Covariate)
#
#   df_low <- as.data.frame(lowCI) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "lowCI", -Covariate)
#
#   df_high <- as.data.frame(highCI) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "highCI", -Covariate)
#
#   # Merge together, ensure same ordering
#   df_betas <- df_betas %>%
#     left_join(df_low, by = c("Covariate", "Signature")) %>%
#     left_join(df_high, by = c("Covariate", "Signature")) %>%
#     mutate(
#       Covariate = factor(Covariate, levels = unique(Covariate)),
#       Signature = as.factor(Signature),
#       markX = ifelse(lowCI < 0 & highCI > 0, TRUE, FALSE)
#     )
#
#   max_abs_beta <- max(abs(df_betas$Beta), na.rm = TRUE)
#
#   plot_out <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta)) +
#     geom_tile(color = "black", width = 0.95, height = 0.95) + # black borders
#     scale_fill_gradientn(
#       colours = c("#2166AC", "white", "#B2182B"),
#       limits = c(-max_abs_beta, max_abs_beta)
#     ) +
#     # Print Beta values
#     geom_text(aes(label = round(Beta, 2)), size = 3) +
#     # Print "X" where CI includes 0
#     geom_text(
#       data = subset(df_betas, markX),
#       aes(label = "X"),
#       color = "black", fontface = "bold", size = 7, vjust=0.5
#     ) +
#     scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
#     scale_x_discrete(position = "top") +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1, size = 10),
#       axis.title = element_blank(),
#       plot.margin = margin(0, 0, 0, 0),
#       axis.ticks.x = element_blank()
#     )
#
#   return(plot_out)
# }
#
# plot_betas3(BetasCIs$mean, BetasCIs$lowCI, BetasCIs$highCI)
# plot_betas4 <- function(Betas_sol, lowCI, highCI) {
#   require(dplyr)
#   require(tidyr)
#   require(ggplot2)
#
#   df_betas <- as.data.frame(Betas_sol) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "Beta", -Covariate)
#
#   df_low <- as.data.frame(lowCI) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "lowCI", -Covariate)
#
#   df_high <- as.data.frame(highCI) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "highCI", -Covariate)
#
#   df_betas <- df_betas %>%
#     left_join(df_low, by = c("Covariate", "Signature")) %>%
#     left_join(df_high, by = c("Covariate", "Signature")) %>%
#     mutate(
#       Covariate = factor(Covariate, levels = unique(Covariate)),
#       Signature = as.factor(Signature),
#       Beta = ifelse(lowCI < 0 & highCI > 0, NA, Beta)
#     )
#
#   max_abs_beta <- max(abs(df_betas$Beta), na.rm = TRUE)
#
#   oldopt <- options()
#   options(ggplot2.continuous.colour="ignore", ggplot2.continuous.fill="ignore")
#
#   plot_out <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta)) +
#     geom_tile(color = "black", width = 1, height = 1, na.rm = TRUE) +
#     scale_fill_gradientn(
#       colours = c("#2166AC", "white", "#B2182B"),
#       limits = c(-max_abs_beta, max_abs_beta),
#       na.value = "white"
#     ) +
#     geom_text(aes(label = ifelse(is.na(Beta), "", round(Beta, 2))), size = 3, na.rm = TRUE) +
#     scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
#     scale_x_discrete(position = "top") +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1, size = 10),
#       axis.title = element_blank(),
#       plot.margin = margin(0, 0, 0, 0),
#       axis.ticks.x = element_blank()
#     )
#
#   options(oldopt)
#
#   return(plot_out)
# }
# plot_betas5 <- function(Betas_sol, lowCI, highCI) {
#
#   df_betas <- as.data.frame(Betas_sol) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "Beta", -Covariate)
#
#   df_low <- as.data.frame(lowCI) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "lowCI", -Covariate)
#
#   df_high <- as.data.frame(highCI) %>%
#     tibble::rownames_to_column("Covariate") %>%
#     tidyr::gather(key = "Signature", value = "highCI", -Covariate)
#
#   df_betas <- df_betas %>%
#     left_join(df_low, by = c("Covariate", "Signature")) %>%
#     left_join(df_high, by = c("Covariate", "Signature")) %>%
#     mutate(
#       Covariate = factor(Covariate, levels = unique(Covariate)),
#       Signature = as.factor(Signature),
#       Beta = ifelse(lowCI < 0 & highCI > 0, NA, Beta)
#     )
#
#   max_abs_beta <- max(abs(df_betas$Beta), na.rm = TRUE)
#
#   oldopt <- options()
#   options(ggplot2.continuous.colour="ignore", ggplot2.continuous.fill="ignore")
#
#   plot_out <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta)) +
#     geom_tile(color = "black", width = 1, height = 1, na.rm = TRUE) +
#     scale_fill_gradientn(
#       colours = c("#2166AC", "white", "#B2182B"),
#       limits = c(-max_abs_beta, max_abs_beta),
#       na.value = "gray90"
#     ) +
#     geom_text(aes(label = ifelse(is.na(Beta), "", round(Beta, 2))), size = 3, na.rm = TRUE) +
#     scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
#     scale_x_discrete(position = "top") +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1, size = 10,
#                                  colour = "black"),
#       axis.title = element_blank(),
#       plot.margin = margin(0, 0, 0, 0),
#       axis.ticks.x = element_blank(), axis.text.y = element_blank()
#     )
#
#   options(oldopt)
#
#   return(plot_out)
# }
#
# plot_sigs <- function(signatures,
#          lowCI = NULL,
#          highCI = NULL,
#          palette = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) {
#
#   signatures <- as.matrix(signatures)
#   names_sig <- rownames(signatures)
#   df_plot <- data.frame(signatures) %>%
#     dplyr::mutate(
#       Channel = names_sig,
#       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
#                       function(x) paste0(x[c(1,3,7)], collapse = "")),
#       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
#                        function(x) paste0(x[c(3,4,5)], collapse = "")),
#       Mutation = as.factor(Mutation)
#     ) %>%
#     tidyr::gather(key = "Sig", value = "Prob", -Channel, -Triplet, -Mutation)
#
#   if(!is.null(lowCI) & !is.null(highCI)){
#     df_plot <- df_plot %>%
#       dplyr::left_join(data.frame(lowCI) %>%
#                          dplyr::mutate(Channel = names_sig,
#                                        Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
#                                                        function(x) paste0(x[c(1,3,7)], collapse = "")),
#                                        Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
#                                                         function(x) paste0(x[c(3,4,5)], collapse = "")),
#                                        Mutation = as.factor(Mutation)) %>%
#                          tidyr::gather(key = "Sig", value = "lowCI", -Channel, -Triplet, -Mutation),
#                        by = c("Channel", "Triplet", "Mutation", "Sig")) %>%
#       dplyr::left_join(data.frame(highCI) %>%
#                          dplyr::mutate(Channel = names_sig,
#                                        Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
#                                                        function(x) paste0(x[c(1,3,7)], collapse = "")),
#                                        Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
#                                                         function(x) paste0(x[c(3,4,5)], collapse = "")),
#                                        Mutation = as.factor(Mutation)) %>%
#                          tidyr::gather(key = "Sig", value = "highCI", -Channel, -Triplet, -Mutation),
#                        by = c("Channel", "Triplet", "Mutation", "Sig"))
#   }
#
#   p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Triplet, y = Prob, fill = Mutation))+
#     ggplot2::geom_bar(stat = "identity", width = 0.7) +
#     ggplot2::facet_grid(Sig~Mutation, scales = "free", switch = "y") + # Horizontal facet labels on left
#     ggplot2::theme_minimal()+
#     ggplot2::scale_fill_manual(values = palette)+
#     ggplot2::theme(
#       legend.position = "none",
#       axis.title = ggplot2::element_blank(),
#       axis.text.y = ggplot2::element_blank(),
#       axis.text.x = ggplot2::element_text(angle = 90, color = "gray35",
#                                           vjust = .5, size = 6.5, margin = ggplot2::margin(t = -4)),
#       panel.grid = ggplot2::element_blank(),
#       panel.spacing.x = ggplot2::unit(0, "lines"),
#       panel.spacing.y = ggplot2::unit(0,"lines")
#     )
#
#   if(!is.null(lowCI) & !is.null(highCI)){
#     p <- p +
#       ggplot2::geom_linerange(ggplot2::aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey65")
#   }
#
#   return(p)
#
# }
#
# pBeta <- plot_betas6(BetasCIs$mean[, mu>0.01], BetasCIs$lowCI[, mu>0.01], BetasCIs$highCI[, mu>0.01])
# pSig <- plot_sigs(CIs$mean[, mu>0.01], CIs$lowCI[, mu>0.01], CIs$highCI[, mu>0.01])
# pSig + pBeta





















