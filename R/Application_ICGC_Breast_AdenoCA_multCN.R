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
  library(dplyr)
  library(tidyr)
  library(ggplot2)

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
rerun <- FALSE
#registerDoParallel(n_start_points)
maxiter <- 2000
K <- length(sigs_to_use) # 13 in this case

if (rerun) {
  #--------- Model 1 - Keep Signatures fixed to the values in Monganella et al
  #-- M1.1 Copy-Number only
  gr_Mutations_void <- gr_Mutations
  mcols(gr_Mutations_void) <- mcols(gr_Mutations_void)[, c(1, 2, 3, 4)]
  Betas_start_void <- matrix(0, nrow = 1, ncol = K)
  SignalTrack_void <- as.matrix(SignalTrack[, c(1)])

  set.seed(10)
  out_m1.1 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_void,
                                    SignalTrack = SignalTrack_void,
                                    CopyTrack = CopyTrack,
                                    K = K,
                                    init = SigPoisProcess_mult.init(R_start = Sigs,
                                                                    Betas_start = Betas_start_void),
                                    controls = SigPoisProcess_mult.controls(maxiter = maxiter),
                                    update_Betas = FALSE,
                                    update_Signatures = FALSE)
  saveRDS(out_m1.1,
          file = "~/SigPoisProcess/output/Application_Breast/out_FixSig_noCovariates.rds.gzip",
          compress = "gzip")

  #-- M1.2 All covariates
  set.seed(10)
  out_m1.2 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                               SignalTrack = SignalTrack,
                               CopyTrack = CopyTrack,
                               K = K,
                               init = SigPoisProcess_mult.init(R_start = Sigs),
                               controls = SigPoisProcess_mult.controls(maxiter = maxiter, tol = 1e-6),
                               update_Signatures = FALSE)

  saveRDS(out_m1.2,
          file = "~/SigPoisProcess/output/Application_Breast/out_FixSig_Covariates.rds.gzip",
          compress = "gzip")

  #--------- Model 2 - Let the model learn the signatures.
  # We try 10 different starting points

  #-- M2.1 Copy-Number only
  set.seed(10)
  out_m2.1 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations_void,
                          SignalTrack = SignalTrack_void,
                          CopyTrack = CopyTrack,
                          K = K,
                          init = SigPoisProcess_mult.init(Betas_start = Betas_start_void),
                          controls = SigPoisProcess_mult.controls(maxiter = maxiter),
                          update_Betas = FALSE)
  saveRDS(out_m2.1,
          file = "~/SigPoisProcess/output/Application_Breast/out_FreeSig_noCovariates.rds.gzip",
          compress = "gzip")


  #-- M2.2 All covariates
  set.seed(10)
  out_m2.2 <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                                    SignalTrack = SignalTrack,
                                    CopyTrack = CopyTrack,
                                    K = K,
                                    controls = SigPoisProcess_mult.controls(maxiter = maxiter))
  saveRDS(out_m2.2,
          file = "~/SigPoisProcess/output/Application_Breast/out_FreeSig_Covariates.rds.gzip",
          compress = "gzip")


}


#------ Find the Maximum-a-posteriori
set.seed(42)
maxiter <- 2000
start_MAP <- Sys.time()
outMAP <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                                SignalTrack = SignalTrack,
                                CopyTrack = CopyTrack,
                                K = 13,
                                controls = SigPoisProcess_mult.controls(maxiter = maxiter,
                                                                        tol = 1e-6))
end_MAP <- Sys.time()
outMAP$time <- end_MAP - start_MAP

saveRDS(outMAP,
        file = "~/SigPoisProcess/output/Application_Breast/Result_MAP2_BreastCopyNum.rds.gzip",
        compress = "gzip")


#----- Run now the mcmc using map as as starting point
nsamples <- 1500
burnin <- 100
# Set inital quantities
init <- SigPoisProcess_mult.init(R_start = outMAP$Signatures,
                                 Theta_start = outMAP$Thetas,
                                 Betas_start = outMAP$Betas,
                                 Mu_start = outMAP$Mu,
                                 Sigma2_start = outMAP$Sigma2)
controls <- SigPoisProcess_mult.controls(nsamples = nsamples, burnin = burnin)

set.seed(42)
start_MCMC <- Sys.time()
outMCMC <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                                 SignalTrack = SignalTrack,
                                 CopyTrack = CopyTrack,
                                 K = K,
                                 method = "mcmc",
                                 init = init,
                                 controls = controls)
end_MCMC <- Sys.time()
outMCMC$time <- end_MCMC - start_MCMC

saveRDS(outMCMC,
        file = "~/SigPoisProcess/output/Application_Breast/Result_MCMC2_BreastCopyNum.rds.gzip",
        compress = "gzip")

if(FALSE){
outMAP <- readRDS("~/SigPoisProcess/output/Application_Breast/Result_MAP2_BreastCopyNum.rds.gzip")
outMCMC <- readRDS("~/SigPoisProcess/output/Application_Breast/Result_MCMC_BreastCopyNum.rds.gzip")

CompressiveNMF::plot_SBS_signature(outMAP$Signatures) + plot_betas(outMAP$Betas)
CompressiveNMF::plot_SBS_signature(outMAP$Signatures[, outMAP$Mu > 0.01]) + plot_betas2(outMAP$Betas[, outMAP$Mu > 0.01])

bigProd <- outMAP$Signatures[gr_Mutations$channel, ] * t(outMAP$Thetas[, gr_Mutations$sample]) * exp(as.matrix(mcols(gr_Mutations)[, -c(1:3)]) %*% outMAP$Betas)

bigProbs <- t(apply(bigProd, 1, function(x) x/sum(x)))

mm <- apply(bigProbs, 1, which.max)
cbind(table(mm), round(outMAP$Mu, 3))







tt <- data$gr_Mutations[data$gr_Mutations$sample == "Patient_01" & data$gr_Mutations$channel =="A[C>A]A"]
start(tt)

match_to_RefSigs(outMAP$Signatures, Cosmic_Sigs)
round(outMAP$Mu, 4)

p <- sample(1:11, 1); k <- sample(1:13, 1)
plot(outMCMC$chain$BETASchain[, "H3K4me3", 2], type = "l")
plot(outMCMC$chain$BETASchain[, "RepliTime", 2], type = "l")

plot(outMCMC$chain$BETASchain[-c(1:200), "H3K4me3", 2], outMCMC$chain$logPostchain[-c(1:200)])
plot(outMCMC$chain$BETASchain[-c(1:200), "RepliTime", 2], outMCMC$chain$logPostchain)


plot(outMCMC$chain$BETASchain[, "H3K4me3", 2], outMCMC$chain$logPostchain)
plot(outMCMC$chain$BETASchain[, "RepliTime", 2], outMCMC$chain$logPostchain)
plot(outMCMC$chain$BETASchain[, "H3K4me3", 2], outMCMC$chain$BETASchain[, "RepliTime", 2])


pairs(outMCMC$chain$BETASchain[, c("RepliTime", "H3K4me3"), c(2,3)])
plot(outMCMC$chain$BETASchain[, "H3K4me3", 2], type = "l")
plot(outMCMC$chain$BETASchain[, "H3K4me3", 3], type = "l")
plot(outMCMC$chain$BETASchain[, "RepliTime", 3], type = "l")

plot(outMCMC$chain$BETASchain[, "H3K4me3", 3], outMCMC$chain$BETASchain[, "RepliTime", 3])

# Effective Sample Sizes
plot(outMCMC$chain$logPostchain, type = "l")

plot(outMCMC$chain$MUchain[, 1], type = "l", ylim =c(0 , max(outMCMC$chain$MUchain) + 1))
for(k in 2:13){
  lines(outMCMC$chain$MUchain[, k], type = "l", col = k)
}

library(coda)
burnin <- 200
effectiveSize(outMCMC$chain$MUchain[-c(1:burnin), ])
plot(outMCMC$chain$MUchain[, 10], type = "l")

get_PosteriorEffectiveSize(outMCMC$chain$BETASchain[-c(1:burnin), ,])
ci_beta <- get_PosteriorCI(outMCMC$chain$BETASchain[-c(1:burnin), ,])
ci_sigs <- get_PosteriorCI(outMCMC$chain$SIGSchain[-c(1:burnin), ,])
ci_theta <- get_PosteriorCI(outMCMC$chain$THETAchain[-c(1:burnin), ,])

p_betas <- plot_betasCI(ci_beta$mean, ci_beta$lowCI, ci_beta$highCI)
p_sigs <- CompressiveNMF:::plot_SBS_signature(ci_sigs$mean,
                                              lowCI = ci_sigs$lowCI,
                                              highCI = ci_sigs$highCI)

p_sigs + p_betas
match_to_RefSigs(ci_sigs$mean, Cosmic_Sigs)


LambdaPred <- rowSums(Reconstruct_Lambda(SignalTrack = SignalTrack,
                                         CopyTrack = CopyTrack,
                                         R = outMAP$Signatures,
                                         Theta = outMAP$Thetas,
                                         Betas = outMAP$Betas))
LambdaTrack <- gr_SignalTrack_2kb_std
LambdaTrack$LambdaPred <- LambdaPred

# # plot it
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)
over <- findOverlaps(LambdaTrack, hg19_Mb)
df_track <- as.data.frame(LambdaTrack) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred))

df_tumor_all <- as.data.frame(gr_Mutations) %>%
  mutate(region = subjectHits(findOverlaps(gr_Mutations, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n())
plot(df_tumor_all$region, slider::slide_dbl(df_tumor_all$n, mean, .before = 1, .after = 1))
lines(df_track$region, slider::slide_dbl(df_track$Lambda, mean, .before = 1, .after = 1),
      col = "red")

plot(df_tumor_all$region, df_tumor_all$n)
lines(df_track$region, df_track$Lambda, col = "red")

plot(df_tumor_all$region[1:500], df_tumor_all$n[1:500])
lines(df_track$region[1:500], df_track$Lambda[1:500], col = "red")

which.max(df_track$Lambda)

tt <- df_tumor_all %>%
  left_join(df_track, by = "region")
plot(tt$n, tt$Lambda)
abline(a = 0, b = 1)

plot(tt$region, tt$n)
lines(tt$region, tt$Lambda, col = "red")

plot(tt$n, tt$Lambda, ylim = c(0, 100), xlim = c(0,100))
abline(a=0, b = 1)


#------------- Make now the plot for Breast cancer
set.seed(10)
outFixed <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
                                  SignalTrack = SignalTrack,
                                  CopyTrack = CopyTrack,
                                  K = K,
                                  init = SigPoisProcess_mult.init(R_start = Sigs),
                                  controls = SigPoisProcess_mult.controls(maxiter = 200, tol = 1e-6),
                                  update_Signatures = FALSE)

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
region_start <- 15e6
region_end   <- 22e6
ChannelProbs <- calculate_ChannelProbs_region(chr,
                                              region_start,
                                              region_end,
                                              j = j,
                                              outFixed,
                                              gr_SignalTrack_2kb_std, CopyTrack)

i <- "A[C>A]C"

dimnames(BaseNMF$mapOutput$R) <- dimnames(Sigs)
colnames(BaseNMF$mapOutput$Theta) <- colnames(mutMatrix)
rownames(BaseNMF$mapOutput$Theta) <- colnames(Sigs)
pr_const <- 2000 *BaseNMF$mapOutput$R[i, ] * BaseNMF$mapOutput$Theta[, j]/(sum(gr_SignalTrack$bin_weight))
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
  ylab("Topographic covariate")+
  xlab(paste0("Genomic regions in ",  chr, " (2kb)"))

p_intensity/p_track + plot_layout(heights = c(2, 1.5))

# Now plot just the mutational signatures
p_sigs <- CompressiveNMF::plot_SBS_signature(Sigs[, c("SBS1", "SBS3", "SBS5", "SBS8")])+
  theme(axis.text.x = ggplot2::element_blank())
colnames(outFixed$Betas) <- colnames(Sigs)
p_betas <- plot_betas(outFixed$Betas[, c("SBS1", "SBS3", "SBS5", "SBS8")])+
  theme(axis.text.y = ggplot2::element_blank())
p_sigs + p_betas+ plot_layout(widths = c(3.5, 2))



}
# plot(outMCMC$chain$BETASchain[,1, 1], type = "l",
#      ylim =c(-0.01 + min(outMCMC$chain$BETASchain) , max(outMCMC$chain$BETASchain) +0.01))
# for(k in 1:13){
#   print(k)
#   for(p in 1:11){
#     lines(outMCMC$chain$BETASchain[, p, k], type = "l", col = k)
#     readline(prompt="Press [enter] to continue")
#   }
# }
#


# Make now the plot from scatch

# Build table for supplement
df_histones <- read_tsv("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/Files_to_use_breast.tsv") %>%
  as.data.frame()


colnames(df_histones)
df_histones %>%
  mutate(Repository = "ENCODE") %>%
  dplyr::select(Repository, target, output_type, classification, accession) %>%
  mutate(accession = paste0(accession, ".bigWig")) %>%
  arrange(classification, target)





# out_m1.2 <- readRDS(file = "~/SigPoisProcess/output/Application_Breast/out_FixSig_Covariates.rds.gzip")
# CompressiveNMF::plot_SBS_signature(out_m1.2$Signatures) + plot_betas(out_1.2$Betas)
# round(out_m2.2$Mu, 3)
# diff(out_m2.2$sol$trace)
# plot(out_m2.2$sol$iter)
# plot(out_m1.2$sol$trace)
#
# match_to_RefSigs(out_m2.2$Signatures, Cosmic_Sigs)



#
# gr_Mutations <- gr_tumor_2kb_std[gr_tumor_2kb_std$sample %in% patients_to_use]
# gr_Mutations$bin_weight <- NULL
#
#
# colnames(MutMatrix) == colnames(CopyTrack)
#
# table(Mutations$sample == names(copy_number))
# levels(Mutations$sample) == colnames(CopyTrack)
# levels(Mutations$sample) == colnames(MutMatrix)
#
# head(CopyTrack)
#
# set.seed(10)
# out <- SigPoisProcess_multCN(gr_Mutations = gr_Mutations,
#                              SignalTrack = SignalTrack,
#                              CopyTrack = CopyTrack,
#                              K = 10,
#                              controls = SigPoisProcess_mult.controls(maxiter = 100))
#
# out$Signatures
# plot_betas(out$Betas)
# CompressiveNMF::plot_SBS_signature(out$Signatures) + plot_betas(out$Betas)
# match_to_RefSigs(out$Signatures, CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)

#--- Now, plot the signalTrack at Mb Scale

# # Let's use only 30 patients
# samples <- unique(gr_tumor$sample)
# set.seed(10)
# patients_to_use <- sample(samples, size = 30)
# CopyTrack <- CopyTrack[, patients_to_use]



#
# # #--- predict the track in model M2.2
# out <- out_m2.2
# LambdaPred <- rowSums(Reconstruct_Lambda(SignalTrack, CopyTrack[, colnames(out$Thetas)], out$Signatures,
#                                  out$Thetas, out$Betas))
# LambdaTrack <- gr_SignalTrack_2kb_std
# LambdaTrack$LambdaPred <- LambdaPred
#
# # # plot it
# genome <- BSgenome.Hsapiens.UCSC.hg19
# chrom_lengths <- seqlengths(genome)[1:23]
# hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
# over <- findOverlaps(LambdaTrack, hg19_Mb)
# df_track <- as.data.frame(LambdaTrack) %>%
#   mutate(region = subjectHits(over)) %>%
#   group_by(region) %>%
#   summarize(Lambda = sum(LambdaPred))
#
# df_tumor_all <- as.data.frame(gr_Mutations) %>%
#   mutate(region = subjectHits(findOverlaps(gr_Mutations, hg19_Mb))) %>%
#   group_by(region) %>%
#   summarise(n = n())
# plot(df_tumor_all$region, slider::slide_dbl(df_tumor_all$n, mean, .before = 1, .after = 1))
# lines(df_track$region, slider::slide_dbl(df_track$Lambda, mean, .before = 1, .after = 1),
#       col = "red")
#
# plot(df_tumor_all$region, df_tumor_all$n)
# lines(df_track$region, df_track$Lambda, col = "red")
#
# tt <- intersect(df_tumor_all$region, df_track$region)
# plot(df_track$Lambda[df_track$region %in% tt],
#      df_tumor_all$n[df_tumor_all$region %in% tt], cex = .2)
# abline(a = 0, b = 1)
#
# # #
# # tt <- intersect(df_tumor_all$region, df_track$region)
# # plot(df_track$Lambda[df_track$region %in% tt],
# #      df_tumor_all$n[df_tumor_all$region %in% tt], cex = .2)
# # abline(a = 0, b = 1)
#
# out_m2.1 <- readRDS(file = "~/SigPoisProcess/output/Application_Breast/out_FreeSig_noCovariates.rds.gzip")
# LambdaPred2 <- rowSums(Reconstruct_Lambda(SignalTrack_void,
#                                           CopyTrack[, colnames(out_m2.1$Thetas)],
#                                           out_m2.1$Signatures,
#                                           out_m2.1$Thetas, out_m2.1$Betas))
# LambdaTrack2 <- gr_SignalTrack_2kb_std
# LambdaTrack2$LambdaPred <- LambdaPred2
# df_track2 <- as.data.frame(LambdaTrack2) %>%
#   mutate(region = subjectHits(over)) %>%
#   group_by(region) %>%
#   summarize(Lambda = sum(LambdaPred))
#
# plot(df_tumor_all$region, df_tumor_all$n)
# lines(df_track$region, df_track$Lambda, col = "red")
# lines(df_track2$region, df_track2$Lambda, col = "blue")
#
# tt <- intersect(df_tumor_all$region, df_track$region)
# plot(df_track2$Lambda[df_track2$region %in% tt],
#      df_tumor_all$n[df_tumor_all$region %in% tt], cex = .2)
# abline(a = 0, b = 1)
#
# tt <- intersect(df_tumor_all$region, df_track2$region)
# plot(df_track2$Lambda[df_track2$region %in% tt],
#      df_tumor_all$n[df_tumor_all$region %in% tt], cex = .2)
# abline(a = 0, b = 1)
#
# sqrt(mean((df_track2$Lambda[df_track2$region %in% tt] - df_tumor_all$n[df_tumor_all$region %in% tt])^2))
# sqrt(mean((df_track$Lambda[df_track$region %in% tt] - df_tumor_all$n[df_tumor_all$region %in% tt])^2))
#
#
#
# sum(outBase$Signatures %*%outBase$Theta)
#
#
# hist(exp(SignalTrack %*% out_m2.2$Betas))
# mutMatrix <- getTotalMutations(gr_Mutations)
# MutEst <- Reconstruct_CountMatrix(SignalTrack, CopyTrack[, colnames(out_m2.2$Thetas)],
#                         R = out_m1.2$Signatures,
#                         Theta = out_m1.2$Thetas,
#                         Betas = out_m1.2$Betas)
#
# plot(MutEst, mutMatrix)
# abline(a = 0, b = 1)
#
#
# CopyTrack2 <- apply(as.matrix(mcols(gr_CopyTrack))[, -1], 2, function(x)
#   x * gr_CopyTrack$bin_weight)
#
# df_track3 <- as.data.frame(LambdaTrack2) %>%
#   mutate(region = subjectHits(over),
#          Lambda = rowSums(CopyTrack),
#          Lambda2 = rowSums(CopyTrack2)) %>%
#   group_by(region) %>%
#   summarize(Lambda = sum(Lambda),
#             Lambda2 = sum(Lambda2))
#
# df_track3 <- as.data.frame(LambdaTrack2) %>%
#   mutate(region = subjectHits(over),
#          Lambda = rowSums(CopyTrack)) %>%
#   group_by(region) %>%
#   summarize(Lambda = sum(Lambda))
#
#
# plot(df_tumor_all$region, scale(df_tumor_all$n), cex = 0.5)
# lines(df_track3$region, scale(df_track3$Lambda2), col = "red")
# lines(df_track3$region, scale(df_track3$Lambda), col = "blue")
#
# tmp <- df_tumor_all%>%
#   left_join(df_track3, by = "region")
#
# plot(tmp$Lambda2, tmp$n); cor(tmp$Lambda2, tmp$n)
# plot(tmp$Lambda, tmp$n); cor(tmp$Lambda, tmp$n)
#
#
# gr_tumor <- gr_Mutations
# samples <- unique(gr_tumor$sample)
# gr_tumor$copy_number <- NA_real_
# for (s in samples) {
#   idx_tumor <- which(gr_tumor$sample == s)
#   idx_copy  <- which(gr_copy$sample == s)
#
#   hits <- findOverlaps(gr_tumor[idx_tumor], gr_copy[idx_copy])
#
#   # Map 'query' (positions in tumor, per sample) to 'subject' in copy
#   # (the indices in hits are relative to the *local* indices above)
#   gr_tumor$copy_number[idx_tumor[queryHits(hits)]] <- gr_copy$score[idx_copy[subjectHits(hits)]]
# }
# gr_tumor$copy_number[is.na(gr_tumor$copy_number)] <- 2
# table(gr_tumor$copy_number, useNA = "ifany")
#
# # Let's do it by patient now.
# df_tumor_all <- as.data.frame(gr_tumor) %>%
#   mutate(region = subjectHits(findOverlaps(gr_tumor, hg19_Mb))) %>%
#   group_by(region) %>%
#   summarise(n = n(),
#             copy = pmin(mean(copy_number), 8))
# plot(df_tumor_all$copy, df_tumor_all$n)
# par(mfrow = c(1, 1))
# plot(df_tumor_all$region, df_tumor_all$n)
# lines(df_tumor_all$region, df_tumor_all$copy * 83.74, type = "l", col = "red")
# lines(df_track3$region, df_track3$Lambda * 1.835e-06 , col = "blue")
#
# lm(df_tumor_all$n ~ df_tumor_all$copy - 1)
# lm(n ~ Lambda - 1, df_tumor_all %>% left_join(df_track3, by = "region"))
# plot(df_tumor_all$copy, df_tumor_all$n)
#
# plot(df_tumor_all$region, scale(df_tumor_all$n))
# lines(df_tumor_all$region, scale(df_tumor_all$copy), type = "l", col = "red")
#
# plot(df_tumor_all$region[1:400], scale(df_tumor_all$copy[1:400]))
# lines(df_tumor_all$region[1:400][-c(1,2,399,400)], zoo::rollmean(scale(df_tumor_all$n[1:400]), k =5))
#
# plot(df_tumor_all$copy, df_tumor_all$n)
#
# plot(df_track3$Lambda[tt], df_tumor_all$copy)
#
#
