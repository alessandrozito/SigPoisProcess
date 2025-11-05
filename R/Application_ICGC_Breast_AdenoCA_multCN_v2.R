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


build_SignalTrack <- function(tilewidth = 2000, type = "tissue"){
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
    file <- paste0("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Breast-Cancer_", type, "_", mark, "_2kb.bigWig")
    gr <- bin_gr(import(file), genome_aggreg, std_chrs)
    gr_SignalTrack$mark <- gr$score
    colnames(mcols(gr_SignalTrack))[ncol(mcols(gr_SignalTrack))] <- mark
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

standardize_covariates <- function(x, q_min = 0.001, q_max = 0.99){
  y <- pmin(x, quantile(x, q_max))
  y <- pmax(y, quantile(x, q_min))
  scale(y)
}

normalize_covariates <- function(x, q_min = 0.001, q_max = 0.99){
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

plot_betas2 <- function(Betas_sol) {

  df_betas <- as.data.frame(Betas_sol) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "Beta", -Covariate) %>%
    dplyr::mutate(
      Covariate = factor(Covariate, levels = unique(Covariate)),
      Signature = as.factor(Signature)
    )

  max_abs_beta <- max(abs(df_betas$Beta), na.rm = TRUE)

  plot_out <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta)) +
    geom_tile(color = "white", width = 0.95, height = 0.95) +
    scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"),
      limits = c(-max_abs_beta, max_abs_beta)
    ) +
    geom_text(aes(label = round(Beta, 2)), size = 3) +
    scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks.x = element_blank()
    )

  return(plot_out)
}


#------

#----- Load the tumor data
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
gr_SignalTrack_2kb <- build_SignalTrack(tilewidth = 2000, type = "tissue")
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
outBase <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                              K = 0,
                                              S = 1e6 * Sigs + 1,
                                              alpha = 1.01, a = 1.01)

mutMatrix <- getTotalMutations(gr_Mutations)
outBase <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                              K = 20,
                                              alpha = 1.01, a = 1.01)
CompressiveNMF::plot_SBS_signature(outBase$Signatures)
# SignalTrack and CopyTrack matrices
SignalTrack <- as.matrix(mcols(gr_SignalTrack_2kb_std))[, -1]
CopyTrack <- apply(as.matrix(mcols(gr_CopyTrack))[, -1], 2,
                   function(x) x * gr_CopyTrack$bin_weight)[, colnames(mutMatrix)]


# Number of starting points for our iterations
n_start_points <- 10
#registerDoParallel(n_start_points)
maxiter <- 2000
K <- length(sigs_to_use) # 13 in this case

# Void quantities for copy-number only
gr_Mutations_void <- gr_Mutations
mcols(gr_Mutations_void) <- mcols(gr_Mutations_void)[, c(1, 2, 3, 4)]
Betas_start_void <- matrix(0, nrow = 1, ncol = K)
SignalTrack_void <- as.matrix(SignalTrack[, c(1)])

#----------------------------------------------------------------
# Model 1 - signature fixed as in Monganella et al (2016)
#----------------------------------------------------------------
run_model_fixed <- FALSE
run_model_free <- TRUE

if (run_model_fixed) {

  #---------------------------------------------------- M1.1 Copy-Number only

  set.seed(10)
  out_m1.1 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_void,
                                    SignalTrack = SignalTrack_void,
                                    CopyTrack = CopyTrack,
                                    K = K,
                                    init = SigPoisProcess_mult.init(R_start = Sigs,
                                                                    Betas_start = Betas_start_void),
                                    controls = SigPoisProcess_mult.controls(maxiter = maxiter,
                                                                            tol = 1e-6),
                                    update_Betas = FALSE,
                                    update_Signatures = FALSE)
  saveRDS(out_m1.1,
          file = "~/SigPoisProcess/output/Application_Breast/Result_MAP_Breast_Fixed_CopyOnly.rds.gzip",
          compress = "gzip")

  #---------------------------------------------------- M1.2 All covariates
  set.seed(10)
  out_m1.2 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                                    SignalTrack = SignalTrack,
                                    CopyTrack = CopyTrack,
                                    K = K,
                                    init = SigPoisProcess_mult.init(R_start = Sigs),
                                    controls = SigPoisProcess_mult.controls(maxiter = maxiter,
                                                                            tol = 1e-6),
                                    update_Signatures = FALSE)

  saveRDS(out_m1.2,
          file = "~/SigPoisProcess/output/Application_Breast/Result_MAP_Breast_Fixed_FullCovariates.rds.gzip",
          compress = "gzip")

}

out_m1.2 <- readRDS("~/SigPoisProcess/output/Application_Breast/Result_MAP_Breast_Fixed_FullCovariates.rds.gzip")
CompressiveNMF::plot_SBS_signature(out_m1.2$Signatures) + plot_betas2(out_m1.2$Betas)
if(run_model_free) {

  #---------------------------------------------------- M2.1 Copy-Number only
  set.seed(10)
  out_m2.1 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_void,
                                    SignalTrack = SignalTrack_void,
                                    CopyTrack = CopyTrack,
                                    K = K,
                                    init = SigPoisProcess_mult.init(Betas_start = Betas_start_void),
                                    controls = SigPoisProcess_mult.controls(maxiter = maxiter,
                                                                            tol = 1e-6),
                                    update_Betas = FALSE)
  saveRDS(out_m2.1,
          file = "~/SigPoisProcess/output/Application_Breast/Result_MAP_Breast_Free_CopyOnly.rds.gzip",
          compress = "gzip")


  #---------------------------------------------------- MAP for all covariates
  set.seed(10)
  outMAP <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                                    SignalTrack = SignalTrack,
                                    CopyTrack = CopyTrack,
                                    K = K,
                                    controls = SigPoisProcess_mult.controls(maxiter = maxiter, tol = 1e-6))
  saveRDS(outMAP,
          file = "~/SigPoisProcess/output/Application_Breast/Result_MAP_Breast_Free_FullCovariates.rds.gzip",
          compress = "gzip")

  #------ MCMC
  nsamples_all <- 10000
  nsamples_current <- 0
  nsamples <- 100
  # Run the model storing the output every 100 iterations
  file_name <- "~/SigPoisProcess/output/Application_Breast/Result_MCMC_Breast_Free_FullCovariates.rds.gzip"
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
    out_temp <-  SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                                       SignalTrack = SignalTrack,
                                       CopyTrack = CopyTrack,
                                       K = K,
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
  }
  #----
}

outMAP <- readRDS("~/SigPoisProcess/output/Application_Breast/Result_MAP_Breast_Free_FullCovariates.rds.gzip")
outMCMC <- readRDS(file_name)
CompressiveNMF::plot_SBS_signature(outMAP$Signatures) + plot_betas2(outMAP$Betas)
CompressiveNMF::plot_SBS_signature(outMCMC$Signatures) + plot_betas2(outMCMC$Betas)
plot(outMCMC$chain$BETASchain[4700:5900, "H3K27me3", "SigN04"], type = "l")
plot(outMCMC$chain$BETASchain[, 1, "SigN12"], type = "l")
plot(outMCMC$chain$BETASchain[, "RepliTime", "SigN06"], type = "l")
plot(outMCMC$chain$MUchain[4900:6000, 10], type = "l")
effectiveSize(outMCMC$chain$MUchain[, 10])
plot(outMCMC$chain$THETAchain[, 7, 10], type = "l")

plot(outMCMC$chain$logPostchain, type = "l")
plot(outMCMC$chain$logLikchain, type = "l")
plot(outMCMC$chain$logPriorchain, type = "l")
acf(outMCMC$chain$logPostchain[1600:2700])
outMCMC$chain$logPostchain[1600:2700]
ss <- apply(outMCMC$chain$SIGSchain[4900:5900, ,], c(2,3), mean)
bb <- apply(outMCMC$chain$BETASchain[4900:5900, ,], c(2,3), mean)
match_to_RefSigs(ss, Cosmic_Sigs)
CompressiveNMF::plot_SBS_signature(ss) + plot_betas2(bb)


for(k in 1:13){
  for(l in 1:11){
    plot(outMCMC$chain$BETASchain[4700:5500, l, k], type = "l")
    print(effectiveSize(outMCMC$chain$BETASchain[4700:5500, l, k]))
    readline(prompt="Press [enter] to continue")
  }
}

CompressiveNMF::plot_SBS_signature(ss) + CompressiveNMF::plot_SBS_signature(outMAP$Signatures)
CompressiveNMF::plot_SBS_signature(ss) + CompressiveNMF::plot_SBS_signature(outMCMC$Signatures)
plot_betas2(outMAP$Betas) + plot_betas2(bb)

LambdaPred <- rowSums(Reconstruct_Lambda(SignalTrack = SignalTrack,
                                         CopyTrack = CopyTrack,
                                         R = outMCMC$Signatures,
                                         Theta = outMCMC$Thetas,
                                         Betas = outMCMC$Betas))
LambdaTrack <- gr_SignalTrack_2kb_std
LambdaTrack$LambdaPred <- LambdaPred

# # plot it
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e4, cut.last.tile.in.chrom = TRUE)
over <- findOverlaps(LambdaTrack, hg19_Mb)
df_track <- as.data.frame(LambdaTrack) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred))

df_tumor_all <- as.data.frame(gr_Mutations) %>%
  mutate(region = subjectHits(findOverlaps(gr_Mutations, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  left_join(df_track, by = "region")

plot(df_tumor_all$region, df_tumor_all$n, type = "l", col = "gray")#cex = 0.5)
lines(df_tumor_all$region, df_tumor_all$Lambda, col = "red")
plot(df_tumor_all$Lambda, df_tumor_all$n, cex = 0.5)
abline(a = 0, b = 1)

df_tumor_all[which.max(df_tumor_all$Lambda), ]



