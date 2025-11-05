library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

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

build_SignalTrack <- function(tilewidth = 10000){
  # Genome aggregated
  genome_aggreg <- genome_5kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                                            tilewidth = tilewidth,
                                            cut.last.tile.in.chrom = TRUE)
  std_chrs <-  paste0("chr", c(1:22, "X"))
  gr_SignalTrack <- add_bin_weights(genome_aggreg)

  # Load all covariates and aggregate them at the desired resolution.
  #---- GC content
  print("GC")
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$GC <- gr$score

  #---- Methylation
  print("Methylation")
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Stomach-AdenoCA_Methylation.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$Methyl <- gr$score

  #---- CTCF
  print("CTCF")
  file <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_CTCF_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$CTCF <- gr$score

  #---- H3K9me3
  print("H3K9me3")
  file <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K9me3_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$H3K9me3 <- gr$score

  #---- H3K36me3
  print("H3K36me3")
  file <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K36me3_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$H3K36me3 <- gr$score

  #---- H3K27me3
  print("H3K27me3")
  file <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K27me3_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$H3K27me3 <- gr$score

  #---- H3K27ac
  print("H3K27ac")
  file <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K27ac_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$H3K27ac <- gr$score

  #---- H3K4me1
  print("H3K4me1")
  file <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K4me1_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$H3K4me1 <- gr$score

  #---- H3K4me3
  print("H3K4me3")
  file <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K4me3_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$H3K4me3 <- gr$score

  #---- Replication Timing
  print("RepliTime")
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$RepliTime <- gr$score

  #---- Nucleosome Occupancy
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
#------

# Load gr_tumor and remove chrY and the mutations in the blacklist
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")
# Remove mutations from blacklist
blacklist <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/hg19-blacklist.v2.bed")
overlaps <- findOverlaps(gr_tumor, blacklist)
gr_tumor <- gr_tumor[-queryHits(overlaps)]
gr_tumor <- gr_tumor[seqnames(gr_tumor) != "chrY"]

# ------ 2 kb aggregation
gr_SignalTrack_2kb <- build_SignalTrack(tilewidth = 2000)

#################### 5 kb aggregation

#-- Get the SignalTrack
gr_SignalTrack_5kb <- build_SignalTrack(tilewidth = 5000)

# Version 1 - Standardize covariates (after capping at 0.999 quantile)
gr_SignalTrack_5kb_std <- gr_SignalTrack_5kb
mcols(gr_SignalTrack_5kb_std)[, -1] <- apply(as.data.frame(mcols(gr_SignalTrack_5kb_std)[, -1]), 2,
                                             function(x) standardize_covariates(x))

# Version 2 - normalize between 0 and 1 (after capping at 0.999 quantile)
gr_SignalTrack_5kb_norm <- gr_SignalTrack_5kb
mcols(gr_SignalTrack_5kb_norm)[, -1] <- apply(as.data.frame(mcols(gr_SignalTrack_5kb_norm)[, -1]), 2,
                                             function(x) normalize_covariates(x))

#-- Match with gr_tumor
gr_tumor_5kb_norm <- merge_with_tumor(gr_tumor, gr_SignalTrack_5kb_norm)
gr_tumor_5kb_std <- merge_with_tumor(gr_tumor, gr_SignalTrack_5kb_std)

# Signatures to be kept fixed
Sigs_to_include <- c("SBS1", "SBS2", "SBS3", "SBS5","SBS13", "SBS17a",
                     "SBS17b", "SBS18", "SBS21", "SBS22a", "SBS26", "SBS44")
Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37[, Sigs_to_include]

#------ Try the normalized solution
bin_weight <- gr_SignalTrack_5kb$bin_weight
mcols(gr_tumor_5kb_norm)[, "bin_weight"] <- NULL
SignalTrack_5kb_norm <- as.matrix(mcols(gr_SignalTrack_5kb_norm))[, -1]

set.seed(10)
out_5kb_norm <- SigPoisProcess_mult(gr_Mutations = gr_tumor_5kb_norm,
                           SignalTrack = SignalTrack_5kb_norm,
                           bin_weight = bin_weight,
                           K = 12,
                           controls = SigPoisProcess_mult.controls(maxiter = 100),
                           update_Signatures = FALSE,
                           init = SigPoisProcess_mult.init(R_start = Sigs))

t(out_5kb_norm$Betas)

#------ Try the standardized solution
bin_weight <- gr_SignalTrack_5kb$bin_weight
mcols(gr_tumor_5kb_std)[, "bin_weight"] <- NULL
SignalTrack_5kb_std <- as.matrix(mcols(gr_SignalTrack_5kb_std))[, -1]

# Start with 3 covariates
cols <- c("tumor", "sample", "channel", "GC",
  "Methyl", "CTCF", "H3K9me3", "H3K36me3",
  "H3K27me3", "H3K27ac", "H3K4me1", "H3K4me3",
  "RepliTime", "NuclOccup")

cols_tmp <- c("tumor", "sample", "channel", "Methyl", "H3K9me3", "RepliTime")

gr_tmp <- gr_tumor_5kb_std
mcols(gr_tmp) <- mcols(gr_tmp)[, cols_tmp]
track_tmp <- SignalTrack_5kb_std[, cols_tmp[-c(1:3)]]

set.seed(10)
out_5kb_std <- SigPoisProcess_mult(gr_Mutations = gr_tmp,
                                    SignalTrack = track_tmp,
                                    bin_weight = bin_weight,
                                    K = 12,
                                    controls = SigPoisProcess_mult.controls(maxiter = 100),
                                    update_Signatures = FALSE,
                                    init = SigPoisProcess_mult.init(R_start = Sigs))

# Try the fully-unsupervised solution
bin_weight <- gr_SignalTrack_5kb$bin_weight
mcols(gr_tumor_5kb_std)[, "bin_weight"] <- NULL
SignalTrack_5kb_std <- as.matrix(mcols(gr_SignalTrack_5kb_std))[, -1]
set.seed(10)
out_5kb_std <- SigPoisProcess_mult(gr_Mutations = gr_tumor_5kb_std,
                                   SignalTrack = SignalTrack_5kb_std,
                                   bin_weight = bin_weight,
                                   K = 12,
                                   controls = SigPoisProcess_mult.controls(maxiter = 300))
plot_96SBSsig(out_5kb_std$Signatures)

plot(out_5kb_std$sol$trace[-c(1:2)])
diff(out_5kb_std$sol$trace)
out_5kb_std$sol$maxdiff

match_to_RefSigs(out_5kb_std$Signatures, CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
t(out_5kb_std$Betas)

# Check now if thw matrix is well reconstructed
MutMatrix <- getTotalMutations(gr_tumor_5kb_std)
ExpW <- bin_weight * exp(SignalTrack_5kb_std %*% out_5kb_std$Betas)

MutMatrix_est <- out_5kb_std$Signatures %*% (out_5kb_std$Thetas * colSums(ExpW))
plot(MutMatrix_est, MutMatrix)
abline(a = 0, b = 1)

#################### 2 kb aggregation
#-- Get the SignalTrack
gr_SignalTrack_2kb <- build_SignalTrack(tilewidth = 2000)

# Version 1 - Standardize covariates (after capping at 0.999 quantile)
gr_SignalTrack_2kb_std <- gr_SignalTrack_2kb
mcols(gr_SignalTrack_2kb_std)[, -1] <- apply(as.data.frame(mcols(gr_SignalTrack_2kb_std)[, -1]), 2,
                                             function(x) standardize_covariates(x))

# Version 2 - normalize between 0 and 1 (after capping at 0.999 quantile)
gr_SignalTrack_2kb_norm <- gr_SignalTrack_2kb
mcols(gr_SignalTrack_2kb_norm)[, -1] <- apply(as.data.frame(mcols(gr_SignalTrack_2kb_norm)[, -1]), 2,
                                              function(x) normalize_covariates(x))

#-- Match with gr_tumor
gr_tumor_2kb_norm <- merge_with_tumor(gr_tumor, gr_SignalTrack_2kb_norm)
gr_tumor_2kb_std <- merge_with_tumor(gr_tumor, gr_SignalTrack_2kb_std)

#------- Standardized
bin_weight <- gr_SignalTrack_2kb$bin_weight
mcols(gr_tumor_2kb_std)[, "bin_weight"] <- NULL
SignalTrack_2kb_std <- as.matrix(mcols(gr_SignalTrack_2kb_std))[, -1]


# Start with 3 covariates
cols <- c("tumor", "sample", "channel", "GC",
          "Methyl", "CTCF", "H3K9me3", "H3K36me3",
          "H3K27me3", "H3K27ac", "H3K4me1", "H3K4me3",
          "RepliTime", "NuclOccup")

cols_tmp <- c("tumor", "sample", "channel", "Methyl", "H3K9me3",
              "RepliTime", "H3K36me3")

gr_tmp <- gr_tumor_2kb_std
mcols(gr_tmp) <- mcols(gr_tmp)[, cols_tmp]
track_tmp <- SignalTrack_2kb_std[, cols_tmp[-c(1:3)]]

# Start from CompressiveNMF
set.seed(10)
MutMatrix <- SigPoisProcess::getTotalMutations(gr_tmp)
outStart <- CompressiveNMF::CompressiveNMF_map(MutMatrix)
CompressiveNMF::CompressiveNMF_map(MutMatrix)

R_start <- outStart$Signatures
Theta_start <- outStart$Theta/sum(bin_weight)

set.seed(10)
out_2kb_std <- SigPoisProcess_mult(gr_Mutations = gr_tmp,
                                   SignalTrack = track_tmp,
                                   bin_weight = bin_weight,
                                   K = ncol(R_start),
                                   controls = SigPoisProcess_mult.controls(maxiter = 30),
                                   init = SigPoisProcess_mult.init(R_start = R_start,
                                                                   Theta_start = Theta_start))

plot(out_2kb_std$sol$trace[-1])
t(out_2kb_std$Betas)
plot_96SBSsig(out_2kb_std$Signatures)


PrecondSigma <- get_Precoditioning_matrices(out_2kb_std, track_tmp, bin_weight)

X <- as.matrix(mcols(gr_tmp)[, -c(1:3)])
channel_id <- as.numeric(gr_tmp$channel)
gr_tmp$sample <- factor(gr_tmp$sample, levels = colnames(out_2kb_std$Thetas))
sample_id <- as.numeric(gr_tmp$sample)

set.seed(10)
mala_res_all <- .PoissonProcess_Bayes_MALA(R_start = out_2kb_std$Signatures,
                                           Theta_start = out_2kb_std$Thetas,
                                           Betas_start = out_2kb_std$Betas,
                                           Mu_start = c(out_2kb_std$Mu),
                                           X = X,
                                           SignalTrack = track_tmp,
                                           bin_weight = bin_weight,
                                           nsamples = 200,
                                           burn_in = 100,
                                           channel_id = channel_id - 1,
                                           sample_id = sample_id - 1,
                                           SigPrior = out_2kb_std$PriorParams$SigPrior,
                                           PrecondSigma = PrecondSigma,
                                           eps_step = 1.3,
                                           update_R = TRUE,
                                           update_Theta = TRUE,
                                           update_Betas = TRUE,
                                           a = 1.1,
                                           a0 = 1.1 * 37 + 1,
                                           b0 = 0.01 * 1.1 * 37)

save(gr_tmp, X, track_tmp, bin_weight, out_2kb_std, mala_res_all,
     file = "../output/Application_Stomach/Application_for_jeff.rdata")

plot(mala_res_all$BETASchain[, 1, 2], type = "l")

library(coda)
burnin <- c(1:50)
1-rejectionRate(as.mcmc(mala_res_all$BETASchain[-burnin, 2, 1]))


plot(mala_res_all$BETASchain[-burnin, 2, 10], type = "l")

apply(mala_res_all$BETASchain[-burnin, ,], c(2,3), function(x) 1-rejectionRate(as.mcmc(x)))

plot(mala_res_all$MUchain[-burnin, 2], type = "l")

plot(mala_res_all$MUchain[-burnin, 1], type = "l",
     ylim = c(0, max(mala_res_all$MUchain[-burnin, ])))
for(j in 2:10) {
  lines(mala_res_all$MUchain[-burnin, j], type = "l", col = j)
}
effectiveSize(mala_res_all$MUchain[-burnin, ])

p <- 1
k <- 1
plot(mala_res_all$BETASchain[-burnin, p, 1], type = "l",
     ylim = c(min(mala_res_all$BETASchain[-burnin, p, ]),
              max(mala_res_all$BETASchain[-burnin, p, ])))
for(k in 2:10) {
  #k <- k+1
  lines(mala_res_all$BETASchain[-burnin, p, k], type = "l", col = k)
}

Betas_start

apply(mala_res_all$BETASchain[-burnin, ,], c(2,3), function(x) effectiveSize(x))

effectiveSize(mala_res_all$BETASchain[-burnin, p, ])

plot(mala_res_all$THETAchain[-burnin, 4, 8] * sum(bin_weight), type = "l")

effectiveSize(mala_res_all$THETAchain[-burnin, 1, 3])
apply(mala_res_all$THETAchain[-burnin, , ] * sum(bin_weight),
      c(2,3), function(x) effectiveSize(x))
apply(mala_res_all$SIGSchain[-burnin, , ], c(2,3), function(x) effectiveSize(x))

apply(mala_res_all$SIGSchain[-burnin, , ], c(2,3), function(x) effectiveSize(x))
effectiveSize(mala_res_all$SIGSchain[-burnin, 1, 5])

SigMean <- apply(mala_res_all$SIGSchain[-burnin, ,], c(2, 3), mean)
rownames(SigMean) <- rownames(CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
colnames(SigMean) <- paste0("SigN", sprintf("%02d", 1:10))
CompressiveNMF::plot_SBS_signature(SigMean)



#################### at 2 kb aggregation, test the following.

# Fix the signatures to those in cosmic (Stomach-AdenoCA)
Sigs_to_include <- c("SBS1", "SBS2", "SBS3", "SBS5","SBS13", "SBS17a",
                     "SBS17b", "SBS18", "SBS21", "SBS22a", "SBS26", "SBS44")
Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37[, Sigs_to_include]

# Version 1 - Standardize covariates (after capping at 0.999 quantile)#-- Get the SignalTrack
gr_SignalTrack_2kb <- build_SignalTrack(tilewidth = 2000)
gr_SignalTrack_2kb_std <- gr_SignalTrack_2kb
mcols(gr_SignalTrack_2kb_std)[, -1] <- apply(as.data.frame(mcols(gr_SignalTrack_2kb_std)[, -1]), 2,
                                             function(x) standardize_covariates(x))

# Merge with gr_tumor
gr_tumor_2kb_std <- merge_with_tumor(gr_tumor, gr_SignalTrack_2kb_std)
bin_weight <- gr_SignalTrack_2kb$bin_weight
mcols(gr_tumor_2kb_std)[, "bin_weight"] <- NULL
SignalTrack_2kb_std <- as.matrix(mcols(gr_SignalTrack_2kb_std))[, -1]

# 1) Model with no covariates (standard NMF)
# Constant model
MutMatrix <- getTotalMutations(gr_tumor_2kb_std)
set.seed(10)
SigPrior <- 1e6 * Sigs[, Sigs_to_include] + 1
BaseNMF <- CompressiveNMF::CompressiveNMF_map(MutMatrix,
                                              S = SigPrior,
                                              K = 0, a = 1.01, alpha = 1.01)

# 2) Model with multiplicative effects
bin_weight <- gr_SignalTrack_2kb$bin_weight
mcols(gr_tumor_2kb_std)[, "bin_weight"] <- NULL
SignalTrack_2kb_std <- as.matrix(mcols(gr_SignalTrack_2kb_std))[, -1]

set.seed(10)
MultNMF_std <- SigPoisProcess_mult(gr_Mutations = gr_tumor_2kb_std,
                                   SignalTrack = SignalTrack_2kb_std,
                                   bin_weight = bin_weight,
                                   K = 12,
                                   controls = SigPoisProcess_mult.controls(maxiter = 100),
                                   update_Signatures = FALSE,
                                   init = SigPoisProcess_mult.init(R_start = Sigs))

diff(MultNMF_std$sol$trace[-1])
t(MultNMF_std$Betas)
Betas <- MultNMF_std$Betas
colnames(Betas) <- Sigs_to_include
round(t(Betas), 3)
plot_96SBSsig(Sigs)

# Make the plot now.
calculate_ChannelProbs_region <- function(chr, region_start, region_end,
                                      j = 1,
                                      MultNMF, gr_SignalTrack) {
  # Extract solution
  Theta <- MultNMF$Thetas
  Betas <- MultNMF$Betas
  R <- MultNMF$Signatures

  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, gr_SignalTrack)
  gr_subset <- gr_SignalTrack[subjectHits(overlaps)]
  # Estimate the probabilities in each subset
  SignalTrack_subset <- as.matrix(mcols(gr_subset)[, -1])
  ExpSolBetas <- exp(SignalTrack_subset %*% Betas)
  # Big Storage
  ChannelProbs <- array(NA, dim = c(nrow(ExpSolBetas), ncol(ExpSolBetas), 96),
               dimnames = list(NULL, colnames(ExpSolBetas), rownames(R)))
  for(i in 1:96){
    SigTheta <- t(R[i, ] * Theta[, j])
    bigProd <- SigTheta[rep(1, nrow(ExpSolBetas)), ] * ExpSolBetas
    sig_probs <-  bigProd/rowSums(bigProd)
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

j <- "DO217814"
chr <- "chr1"
region_start <- 15e6
region_end   <- 22e6
ChannelProbs <- calculate_ChannelProbs_region(chr, region_start,
                                              region_end, j = j, MultNMF_std,
                                              gr_SignalTrack_2kb_std)

i <- "A[C>A]C"
df_mod <- data.frame(ChannelProbs[, , i])
colnames(df_mod) <- colnames(BaseNMF$Signatures)

pr_const <- BaseNMF$Signatures[i, ] * BaseNMF$Theta[, j]
df_const <- data.frame(t(pr_const/sum(pr_const))[rep(1,nrow(ChannelProbs)), ]) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start,
         method = "Standard NMF")

window <- 10
p_plots <- df_mod %>%
  mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
                                                   align = "center"))) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start,
         method = "PoissonProcess") %>%
  bind_rows(df_const) %>%
  gather(key = "Sig", value = "Prob", - region, -method) %>%
  filter(Sig %in% head(names(sort(pr_const, decreasing = TRUE)), 4))%>%
  ggplot() +
  geom_line(aes(x = region, y = Prob, color = Sig, linetype = method), linewidth = 0.7) +
  scale_color_manual(name = "Signature",
                     values = c("gray25", "darkorange", "red", "blue")) +
  ylab("Signature attribution\nprobability") +
  #xlab(paste0("Genomic regions in ",  chr, "-",
  #            region_start,":", region_end, " (2kb)")) +
  xlab(paste0("Genomic regions in ",  chr, " (2kb)")) +
  facet_wrap(~paste0("Mutation type ", i, ", Patient ", j)) +
  theme_minimal(base_size = 10) +
  scale_x_continuous(expand = c(0,0))

# Now, make the plot of the various histone marks below it.
subset_tracks <- get_signal_regions(chr, region_start, region_end, gr_SignalTrack_2kb_std)
subset_tracks$bin_weight <- NULL
df_tracks <- as.data.frame(mcols(subset_tracks)) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start) %>%
  gather(key = "covariate", value = "score", -region)

p_tracks <- ggplot(df_tracks) +
  geom_raster(aes(x = region + 1000, y = covariate, fill = score)) +
  #facet_wrap(.~ pos2, nrow = 1) +
  scale_fill_gradientn(name = "Signal track\n(standardized)",
                       colours =c("antiquewhite", "#F6C866", "#F2AB67", "#EF8F6B",
                                  "#ED7470", "#BF6E97", "#926AC2",
                                  "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
  theme_minimal(base_size = 10) +
  scale_x_continuous(expand = c(0,0)) +
  #theme(axis.title = element_blank())
  ylab("Epigenetic mark")+
  xlab(paste0("Genomic regions in ",  chr, " (2kb)"))



# Now plot just the mutational signatures
p_sigs <- CompressiveNMF::plot_SBS_signature(Sigs[, c("SBS1", "SBS18", "SBS3", "SBS5")])

# Finally, plot the coefficients
df_betas <- as.data.frame(Betas[, c("SBS1", "SBS18", "SBS3", "SBS5")]) %>%
  rownames_to_column("Covariate") %>%
  gather(key = "Signature", value = "Beta", -Covariate) %>%
  dplyr::mutate(
    Covariate = factor(Covariate, levels=unique(Covariate)),
    Signature = as.factor(Signature),
    AltCol = ifelse((as.numeric(Covariate) %% 2) == 0, "gray93", "gray97"))

p_coefs <- ggplot2::ggplot(df_betas, aes(x = Covariate, y = Signature)) +
  ggplot2::geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
  ggplot2::scale_fill_manual(
    values = c("gray93" = "gray93", "gray97" = "gray97"),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_point(aes(fill = Beta, size = abs(Beta)),
                      shape = 21, color = "grey30", stroke = 0.2,
                      na.rm = TRUE) +
  ggplot2::scale_fill_gradientn(name = "Beta",
                                colours =c("blue","lightblue", "white", "red")) +
  #scale_fill_gradient2(mid = "white", low = "blue", high = "red")+
  ggplot2::scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
    axis.title = ggplot2::element_blank(),
    #axis.text.y = element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0),
    axis.ticks.x = ggplot2::element_blank())+
    guides(size = "none")

right_col <- p_plots / p_tracks + plot_layout(heights = c(2, 1))
p_sigs + p_coefs + right_col + plot_layout(widths = c(1.5, 0.8, 2))
p_sigs + p_coefs + plot_layout(widths = c(1.5, 0.5))
#--------------------------------------------------------
j <- "DO217814"
gr_SignalTrack_2kb
gr_tumor_2kb <- merge_with_tumor(gr_tumor, gr_SignalTrack_2kb)

df_tmp <- data.frame(gr_tumor_2kb) %>%
  filter(sample == "DO217814")

df_do <- data.frame(gr_tumor_2kb) %>%
  filter(sample == "DO217814") %>%
  group_by(channel) %>%
  summarise(n = n(),
            H3K9me3 = sum(H3K9me3),
            CTCF = sum(GC))

plot(df_do$H3K9me3, df_do$CTCF)

plot(df_do$H3K9me3, df_do$n)

MM <- df_do %>%
  column_to_rownames("channel") %>%
  as.matrix()
plot_96SBSsig(MM)

#--------------------------------------------------------
# Now, experience with Bayesian sampler.
gr_tmp <- gr_tumor_2kb_std
mcols(gr_tmp) <- mcols(gr_tmp)[, cols_tmp]
track_tmp <- SignalTrack_2kb_std[, cols_tmp[-c(1:3)]]

set.seed(10)
out_2kb_std <- SigPoisProcess_mult(gr_Mutations = gr_tmp,
                                   SignalTrack = track_tmp,
                                   bin_weight = bin_weight,
                                   K = 12,
                                   controls = SigPoisProcess_mult.controls(maxiter = 150),
                                   update_Signatures = FALSE,
                                   init = SigPoisProcess_mult.init(R_start = Sigs))

gr_tmp$sample <- factor(gr_tmp$sample, levels = colnames(out_2kb_std$Thetas))
X <- as.matrix(mcols(gr_tmp)[, -c(1:3)])
channel_id <- as.numeric(gr_tmp$channel)
sample_id <- as.numeric(gr_tmp$sample)

SignalTrack <- track_tmp
SigPrior <- matrix(1.1, 96, 12)

R_start <- out_2kb_std$Signatures
Theta_start <- out_2kb_std$Thetas
Betas_start <- out_2kb_std$Betas

nsamples = 100
burn_in = 100

set.seed(10)
out_2kb_std <- SigPoisProcess_mult(gr_Mutations = gr_tmp,
                                   SignalTrack = track_tmp,
                                   bin_weight = bin_weight,
                                   K = 12,
                                   controls = SigPoisProcess_mult.controls(maxiter = 2))

R_start <- out_2kb_std$Signatures
Theta_start <- out_2kb_std$init$Theta_start
Betas_start <- out_2kb_std$init$Betas_start
set.seed(10)
result <- .PoissonProcess_Bayes(R_start = R_start,
                                Theta_start = Theta_start,
                                Betas_start = Betas_start,
                                X = X,
                                SignalTrack = SignalTrack,
                                bin_weight = bin_weight,
                                nsamples = 3000,
                                burn_in = burn_in,
                                channel_id = channel_id - 1,
                                sample_id = sample_id - 1,
                                SigPrior = matrix(0.5, 96, K),
                                update_R = TRUE,
                                update_Theta = TRUE,
                                update_Betas = TRUE,
                                a = 1,
                                a0 = 1 * 37 + 1,
                                b0 = 0.01 * 1 * 37)

burnin <- 1:2000
resMean <- apply(result$SIGSchain[-burnin,,], c(2,3), mean)
rownames(resMean) <- rownames(CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
colnames(resMean) <- paste0("Sig", 1:K) #dimnames(R_start)
CompressiveNMF::plot_SBS_signature(resMean)

plot(result$BETASchain[ ,3,8], type = "l")
plot(result$THETAchain[ ,1,1], type = "l")
plot(result$MUchain[-burnin, 10], type = "l")

colMeans(result$MUchain[-burnin, ])
mm <- max(result$MUchain[-burnin, ])
plot(result$MUchain[-burnin, 1], type = "l", ylim = c(0, mm))
for(j in 2:K) {
  lines(result$MUchain[-burnin, j], type = "l", col = j)
}
coda::effectiveSize(result$MUchain[-burnin, ])

plot(result$BETASchain[ ,3,8], type = "l")
plot(result$THETAchain[ ,1,1], type = "l")
plot(result$MUchain[-burnin, 10], type = "l")

pp <- 3
mm1 <- max(result$BETASchain[-burnin, pp,])
mm2 <- min(result$BETASchain[-burnin, pp,])
plot(result$BETASchain[-burnin, pp, 1], type = "l", ylim = c(mm2, mm1))
for(j in 2:K) {
  lines(result$BETASchain[-burnin, pp, j], type = "l", col = j)
}

resBetas <- apply(result$BETASchain[-burnin, ,], c(2,3), mean)
colnames(resBetas) <- colnames(resMean)
rownames(resBetas) <- colnames(X)
t(resBetas)


R <- out_2kb_std$Signatures
Theta <- out_2kb_std$Thetas
Betas <- out_2kb_std$Betas

Mu <- out_2kb_std$Mu
ExpXbeta <- exp(X %*% Betas)
Probs_un <- R[channel_id, ] * t(Theta[, sample_id]) * ExpXbeta
Probs <- Probs_un/rowSums(Probs_un)

W <- sample_Categorical(Probs)


Rnew <- sample_R(R, W, SigPrior, channel_id - 1)

(Theta_new <- sample_Theta(t(Theta), Betas, SignalTrack, bin_weight,
                          W, Mu, 1.1, sample_id - 1))


#--------------------------------------------------------
# Let's try the MALA algorithm then, since adaptive metropolis failed.

set.seed(10)
out_2kb_std <- SigPoisProcess_mult(gr_Mutations = gr_tmp,
                                   SignalTrack = track_tmp,
                                   bin_weight = bin_weight,
                                   K = 10,
                                   controls = SigPoisProcess_mult.controls(maxiter = 1000))

plot_96SBSsig(out_2kb_std$Signatures)

Betas <- out_2kb_std$Betas; Betas_start <- Betas
R <- out_2kb_std$Signatures
Theta <- out_2kb_std$Thetas
Mu <- out_2kb_std$Mu
ExpXbeta <- exp(X %*% Betas)

out <- out_2kb_std
Probs_un <- R[channel_id, ] * t(Theta[, sample_id]) * ExpXbeta
Probs <- Probs_un/rowSums(Probs_un)
set.seed(10)
W <- sample_Categorical(Probs)

out <- out_2kb_std
PrecondSigma <- get_Precoditioning_matrices(out, SignalTrack, bin_weight)
set.seed(10)
mala_res <- sample_Betas_MALA(Betas_start = Betas_start,
                              W = W,
                              X = X,
                              SignalTrack = SignalTrack,
                              Theta = t(Theta),
                              bin_weight = bin_weight,
                              PrecondSigma = PrecondSigma,
                              Mu = Mu,
                              nsamples = 100, eps_step = 1.2)

plot(mala_res[, 1, 1], type = "l")

Betas <- out_2kb_std$Betas; Betas_start <- Betas
R_start <- out_2kb_std$Signatures
Theta_start <- out_2kb_std$Thetas
Mu_start <- out_2kb_std$Mu
SigPrior <- matrix(1.1, 96, 10)

X <- as.matrix(mcols(gr_tmp)[, -c(1:3)])
channel_id <- as.numeric(gr_tmp$channel)
gr_tmp$sample <- factor(gr_tmp$sample, levels = colnames(out_2kb_std$Thetas))
sample_id <- as.numeric(gr_tmp$sample)

set.seed(10)
mala_res_all <- .PoissonProcess_Bayes_MALA(R_start = out_2kb_std$Signatures,
                                           Theta_start = out_2kb_std$Thetas,
                                           Betas_start = out_2kb_std$Betas,
                                           Mu_start = out_2kb_std$Mu,
                                           X = X,
                                           SignalTrack = SignalTrack,
                                           bin_weight = bin_weight,
                                           nsamples = 500,
                                           burn_in = 100,
                                           channel_id = channel_id - 1,
                                           sample_id = sample_id - 1,
                                           SigPrior = SigPrior,
                                           PrecondSigma = PrecondSigma,
                                           eps_step = 1.3,
                                           update_R = TRUE,
                                           update_Theta = TRUE,
                                           update_Betas = TRUE,
                                           a = 1.1,
                                           a0 = 1.1 * 37 + 1,
                                           b0 = 0.01 * 1.1 * 37)

plot(mala_res_all$BETASchain[, 1, 4], type = "l")

library(coda)
burnin <- c(1:50)
1-rejectionRate(as.mcmc(mala_res_all$BETASchain[-burnin,2,1]))

mean(apply(mala_res_all$BETASchain[-burnin, ,], c(2,3), function(x)
  1-rejectionRate(as.mcmc(x))))

plot(mala_res_all$MUchain[, 6], type = "l")


plot(mala_res_all$MUchain[-burnin, 1], type = "l",
     ylim = c(0, max(mala_res_all$MUchain[-burnin, ])))
for(j in 2:10) {
  lines(mala_res_all$MUchain[-burnin, j], type = "l", col = j)
}
effectiveSize(mala_res_all$MUchain[-burnin, ])

p <- 1
k <- 1
plot(mala_res_all$BETASchain[-burnin, p, 1], type = "l",
     ylim = c(min(mala_res_all$BETASchain[-burnin, p, ]),
              max(mala_res_all$BETASchain[-burnin, p, ])))
for(k in 2:10) {
  #k <- k+1
  lines(mala_res_all$BETASchain[-burnin, p, k], type = "l", col = k)
}

Betas_start

apply(mala_res_all$BETASchain[-burnin, ,], c(2,3), function(x) effectiveSize(x))

effectiveSize(mala_res_all$BETASchain[-burnin, p, ])

plot(mala_res_all$THETAchain[-burnin, 4, 35], type = "l")

effectiveSize(mala_res_all$THETAchain[-burnin, , ])
apply(mala_res_all$THETAchain[-burnin, , ], c(2,3), function(x) effectiveSize(x))
apply(mala_res_all$SIGSchain[-burnin, , ], c(2,3), function(x) effectiveSize(x))

apply(mala_res_all$SIGSchain[-burnin, , ], c(2,3), function(x) effectiveSize(x))
effectiveSize(mala_res_all$SIGSchain[-burnin, 1, 5])

SigMean <- apply(mala_res_all$SIGSchain[-burnin, ,], c(2, 3), mean)
rownames(SigMean) <- rownames(CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
colnames(SigMean) <- paste0("SigN", sprintf("%02d", 1:10))
CompressiveNMF::plot_SBS_signature(SigMean)

effectiveSize(result$BETASchain[, , 1])

0.000015391857 * crossprod(exp(SignalTrack %*% beta_new) * bin_weight, SignalTrack)
0.000015391857 * crossprod(exp(SignalTrack %*% Betas_start[, 1]) * bin_weight, SignalTrack)

# Adaptive metropolis. One proposal matrix for each
K <-12; p <- 3
Sigmas_prop <- array(NA, dim = c(K, p, p))
for(k in 1:K){
  Sigmas_prop[k, ,] <- diag(p)
}



set.seed(10)
Betas_start = Betas;
resBetas <- sample_Betas(Betas_start = Betas_start,
                         W = W,
                         X = X,
                         SignalTrack = SignalTrack,
                         Theta = t(Theta),
                         bin_weight = bin_weight,
                         Mu = Mu,
                         nsamples = 1000)

plot(resBetas[, 2, 5], type = "l")

plot(density(resBetas[, 1, 1]))
acf(resBetas[, 2, 5])
1 - coda::rejectionRate(coda::as.mcmc(resBetas[, 2, 5]))
library(coda)
effectiveSize(resBetas[, 2, 5])

M <- matrix(
  c(0.0144, -0.0173, 0.1014,
    -0.0173, 0.0209, -0.1221,
    0.1014, -0.1221, 0.7141),
  nrow = 3, byrow = TRUE
)

v1 <- c(-0.1244, 0.1483, -0.8693)
v2 <- c(-0.1225, 0.1491, -0.8702)


bin_weight %*%SignalTrack


#---- EPIGENETIC MARKS
patient <- "DO217814"




Signal
















