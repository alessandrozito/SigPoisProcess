library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
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


build_CopyTrack <- function(gr_tumor, gr_copy, genome_aggreg){
  # Add bin weights to the genome
  gr_CopyTrack <- add_bin_weights(genome_aggreg)
  std_chrs <-  paste0("chr", c(1:22, "X"))
  samples <- unique(gr_tumor$sample)
  # Filter copies in the tumor
  gr_copy_tumor <- gr_copy[gr_copy$sample %in% unique(gr_tumor$sample)]
  # Output to store
  CopyTrack <- matrix(0, nrow = length(genome_aggreg), ncol = length(samples))
  for(s in seq_along(samples)) {
    gr_copy_sample <- bin_gr(gr_copy_tumor[gr_copy_tumor$sample == samples[s]],
                             genome_aggreg, std_chrs)
    # Substitute with 2 if NA. simplifying assumption.
    score <- gr_copy_sample$score
    score[is.na(score)] <- 2
    score[score < 0.1] <- 0.1
    CopyTrack[, s] <- score/2
  }
  colnames(CopyTrack) <- samples
  mcols(gr_CopyTrack) <- cbind("bin_weight" = genome_aggreg_w$bin_weight, CopyTrack)
  return(gr_CopyTrack[genome_aggreg_w$bin_weight > 0])
}

gr_CopyTrack <- build_CopyTrack(gr_tumor, gr_copy, genome_aggreg)

table(is.na(CopyTrack))

table(CopyTrack[genome_aggreg_w$bin_weight == 0, 22], useNA = "ifany")
table(genome_aggreg_w$bin_weight[is.na(CopyTrack[, 11])], useNA = "ifany")

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
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Stomach-AdenoCA_Methylation.bigWig"
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
  file <- "~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig"
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

gr_tumor <- merge_CopyTrack_tumor(gr_tumor, gr_CopyTrack)

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

plot_betas <- function(Betas_sol){
  df_betas <- as.data.frame(Betas_sol) %>%
    rownames_to_column("Covariate") %>%
    gather(key = "Signature", value = "Beta", -Covariate) %>%
    dplyr::mutate(
      Covariate = factor(Covariate, levels=unique(Covariate)),
      Signature = as.factor(Signature),
      AltCol = ifelse((as.numeric(Covariate) %% 2) == 0, "gray93", "gray97"))

  plot_out <- ggplot2::ggplot(df_betas, aes(x = Covariate, y = Signature)) +
    ggplot2::geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
    ggplot2::scale_fill_manual(
      values = c("gray93" = "gray93", "gray97" = "gray97"),
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(aes(fill = Beta, size = abs(Beta)),
                        shape = 21, color = "grey30", stroke = 0.2,
                        na.rm = TRUE) +
    #ggplot2::scale_fill_gradientn(name = "Beta",
    #                              colours =c("blue","lightblue", "white", "red")) +
    scale_fill_gradient2(mid = "white", low = "blue", high = "red", midpoint = 0)+
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
  return(plot_out)
}
#------

# Load gr_tumor and remove chrY and the mutations in the blacklist
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Breast-AdenoCa_snp.rds.gzip")
# Remove mutations from blacklist
blacklist <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/hg19-blacklist.v2.bed")
overlaps <- findOverlaps(gr_tumor, blacklist)
gr_tumor <- gr_tumor[-queryHits(overlaps)]
gr_tumor <- gr_tumor[seqnames(gr_tumor) != "chrY"]

# ------ 2 kb aggregation
gr_SignalTrack_2kb <- build_SignalTrack(tilewidth = 2000, type = "tissue")

# Standardize covariates (after capping at 0.999 quantile)
gr_SignalTrack_2kb_std <- gr_SignalTrack_2kb
mcols(gr_SignalTrack_2kb_std)[, -1] <- apply(as.data.frame(mcols(gr_SignalTrack_2kb_std)[, -1]), 2,
                                             function(x) standardize_covariates(x))

# Match with tumor
gr_tumor_2kb_std <- merge_with_tumor(gr_tumor, gr_SignalTrack_2kb_std)

# Run the model
bin_weight <- gr_SignalTrack_2kb$bin_weight
mcols(gr_tumor_2kb_std)[, "bin_weight"] <- NULL
SignalTrack_2kb_std <- as.matrix(mcols(gr_SignalTrack_2kb_std))[, -1]


#----- Model 1 - Keep signatures fixed at those predicted by the paper,
# Start from the value at compressiveNMF.
Cosmic_Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
sigs_to_use <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8", "SBS20",
                 "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")
# Run the model
Sigs <- Cosmic_Sigs[, sigs_to_use]
set.seed(10)
out_fixedSigs <- SigPoisProcess_mult(gr_Mutations = gr_tumor_2kb_std,
                                     SignalTrack = SignalTrack_2kb_std,
                                     bin_weight = bin_weight,
                                     K = ncol(Sigs),
                                     update_Signatures = FALSE,
                                     controls = SigPoisProcess_mult.controls(maxiter = 1000,
                                                                             tol = 1e-6),
                                     prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                                    c0 =  100, d0 = 1),
                                     init = SigPoisProcess_mult.init(R_start = Sigs))

saveRDS(out_fixedSigs, file = "~/SigPoisProcess/output/Sigfixed_Breast.rds")
out_fixedSigs <- readRDS("~/SigPoisProcess/output/Sigfixed_Breast.rds")

#------ What is the estimate of the mutation rate along the genome?
colnames(out_fixedSigs$Signatures) <- colnames(out_fixedSigs$Betas) <- colnames(Sigs)
p1 <- CompressiveNMF::plot_SBS_signature(out_fixedSigs$Signatures)
p2 <- plot_betas(out_fixedSigs$Betas)
p1 + p2
colnames(out_fixedSigs$Mu) <- colnames(Sigs)
out_fixedSigs$Mu

round(t(out_fixedSigs$Betas), 3)
#---- Plot now the predicted track
ExpTrackPred <- apply(exp(SignalTrack_2kb_std %*% out_fixedSigs$Betas), 2, function(x) bin_weight * x)
LambdaPred <- rep(0, length(bin_weight))
for(j in 1:ncol(out_fixedSigs$Theta)){
  print(j)
  tmp <- out_fixedSigs$Theta[, j] * t(ExpTrackPred)
  for(i in 1:96){
    LambdaPred <- LambdaPred +  c(crossprod(out_fixedSigs$Signatures[i, ], tmp))
  }
}
LambdaTrack <- gr_SignalTrack_2kb_std
LambdaTrack$LambdaPred <- LambdaPred

# Take regions of 1Mb.
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
over <- findOverlaps(LambdaTrack, hg19_Mb)
df_track <- as.data.frame(LambdaTrack) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred))

df_tumor_all <- as.data.frame(gr_tumor) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n())
plot(df_tumor_all$region, df_tumor_all$n)
lines(df_track$region, df_track$Lambda, col = "red")

tt <- intersect(df_tumor_all$region, df_track$region)
plot(df_track$Lambda[df_track$region %in% tt],
     df_tumor_all$n[df_tumor_all$region %in% tt])
abline(a = 0, b = 1)



set.seed(10)
out_FreeSigs <- SigPoisProcess_mult(gr_Mutations = gr_tumor_2kb_std,
                                     SignalTrack = SignalTrack_2kb_std,
                                     bin_weight = bin_weight,
                                     K = 12,
                                     controls = SigPoisProcess_mult.controls(maxiter = 2000,
                                                                             tol = 1e-6),
                                     prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma",
                                                                                    c0 =  100, d0 = 1))
saveRDS(out_FreeSigs, file = "~/SigPoisProcess/output/SigFree_Breast.rds")


diff(out_FreeSigs$sol$trace)
# This just went up. Very strange
Betas_sol <- out_fixedSigs$Betas
colnames(Betas_sol) <- colnames(Sigs)


p1 <- CompressiveNMF::plot_SBS_signature(out_FreeSigs$Signatures)
p2 <- plot_betas(out_FreeSigs$Betas)
p1 + p2
colnames(out_FreeSigs$Mu) <- colnames(out_FreeSigs$Signatures)

match_to_RefSigs(out_FreeSigs$Signatures, Cosmic_Sigs)


#---- Plot now the predicted track
ExpTrackPred_free <- apply(exp(SignalTrack_2kb_std %*% out_FreeSigs$Betas), 2, function(x) bin_weight * x)
LambdaPred_free <- rep(0, length(bin_weight))
for(j in 1:ncol(out_FreeSigs$Theta)){
  print(j)
  tmp <- out_FreeSigs$Theta[, j] * t(ExpTrackPred_free)
  for(i in 1:96){
    LambdaPred_free <- LambdaPred_free +  c(crossprod(out_FreeSigs$Signatures[i, ], tmp))
  }
}
LambdaTrack_free <- gr_SignalTrack_2kb_std
LambdaTrack_free$LambdaPred <- LambdaPred_free

# Take regions of 1Mb.
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
over <- findOverlaps(LambdaTrack_free, hg19_Mb)
df_track_free <- as.data.frame(LambdaTrack_free) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred))

plot(df_tumor_all$region, df_tumor_all$n)
lines(df_track_free$region, df_track_free$Lambda, col = "red")
lines(df_track$region, df_track$Lambda, col = "blue")
plot(df_track$Lambda, df_track_free$Lambda)


#----- Check copy number alteration now
df_copy <- read_tsv("~/20170119_final_consensus_copynumber_donor") %>%
  filter(sampleID %in% samples)

length(unique(df_copy$sampleID))

head(df_copy)

gr_copy <- GRanges(
  seqnames = paste0("chr",df_copy$chr),
  ranges = IRanges(start = df_copy$start, end = df_copy$end),
  strand = "*",
  sample = df_copy$sampleID,
  score = df_copy$value)

# Append it for all subjects.
gr_copy_subject <- gr_copy[gr_copy$sample == samples[j]]
hg19_Mb
overlaps <- findOverlaps(hg19_Mb, gr_copy[gr_copy$sample == samples[j]])
table(table(queryHits(overlaps)))
gr_copy_subject$value[is.na(gr_copy_subject$value)] <- 2
numvar <- coverage(gr_copy_subject, weight = "value")
rm(gr)
# 1kb
print("1kb")
binned_score_1Mb <- binnedAverage(bins = hg19_Mb,
                                  numvar = numvar,
                                  varname = "score")


samples <- unique(gr_tumor$sample)
gr_tumor$copy_number <- NA_real_
for (s in samples) {
  idx_tumor <- which(gr_tumor$sample == s)
  idx_copy  <- which(gr_copy$sample == s)

  hits <- findOverlaps(gr_tumor[idx_tumor], gr_copy[idx_copy])

  # Map 'query' (positions in tumor, per sample) to 'subject' in copy
  # (the indices in hits are relative to the *local* indices above)
  gr_tumor$copy_number[idx_tumor[queryHits(hits)]] <- gr_copy$score[idx_copy[subjectHits(hits)]]
}
gr_tumor$copy_number[is.na(gr_tumor$copy_number)] <- 2
table(gr_tumor$copy_number, useNA = "ifany")

# Let's do it by patient now.
df_tumor_all <- as.data.frame(gr_tumor) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n(),
            copy = mean(copy_number))

plot(df_tumor_all$copy, df_tumor_all$n)
par(mfrow = c(1, 1))
plot(df_tumor_all$region, scale(df_tumor_all$n))
lines(df_tumor_all$region, scale(df_tumor_all$copy), type = "l", col = "blue")

samples <- unique(gr_tumor$sample)
j <- sample(1:113, 1)#j <- 26
df_tumor_all <- as.data.frame(gr_tumor) %>%
  filter(sample == samples[j]) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor[gr_tumor$sample == samples[j]],
                                           hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n(),
            copy = mean(copy_number))
par(mfrow = c(1, 2))
plot(df_tumor_all$region, scale(df_tumor_all$n))
lines(df_tumor_all$region, scale(df_tumor_all$copy), type = "l", col = "red")
boxplot(n ~ round(copy), data = df_tumor_all)
plot(log2(df_tumor_all$copy), df_tumor_all$n)
plot(df_tumor_all$copy, df_tumor_all$n)
fit <- glm(n ~ copy - 1, df_tumor_all %>% mutate(copy = copy/2 + 0.1),
           family = poisson(link = "identity"), start = 1)
fit <- glm(n ~ copy, df_tumor_all %>% mutate(copy = copy/2 + 0.1),
           family = poisson(link = "log"))
plot(predict(fit, type = "response"), df_tumor_all$n)
abline(a = 0, b = 1)
summary(fit)

plot(df_tumor_all$copy, df_tumor_all$n)
# Merge with hg19 genome
samples <- unique(gr_tumor$sample)
gr_tumor$copy_number <- NA_real_
for (s in samples) {
  idx_tumor <- which(gr_tumor$sample == s)
  idx_copy  <- which(gr_copy$sample == s)
  hits <- findOverlaps(gr_tumor[idx_tumor], gr_copy[idx_copy])
  gr_tumor$copy_number[idx_tumor[queryHits(hits)]] <- gr_copy$value[idx_copy[subjectHits(hits)]]
}
gr_tumor$copy_number[is.na(gr_tumor$copy_number)] <- 2
table(gr_tumor$copy_number, useNA = "ifany")

gr_copy_subject <- gr_copy[gr_copy$sample == samples[j]]
seqnames(gr_copy_subject)
hg19_Mb
overlaps <- findOverlaps(hg19_Mb, gr_copy[gr_copy$sample == samples[j]])
table(table(queryHits(overlaps)))
gr_copy_subject$value[is.na(gr_copy_subject$score)] <- 2
numvar <- coverage(gr_copy_subject, weight = "value")
rm(gr)
# 1kb
print("1kb")
binned_score_1Mb <- binnedAverage(bins = hg19_Mb,
                                  numvar = numvar,
                                  varname = "score")


plot(df_tumor_all$n)
plot(binned_score_1Mb$score, type = "l")
lines(scale(binned_score_1Mb[df_tumor_all$region]$score), col = "red")
plot(binned_score_1Mb[df_tumor_all$region]$score, df_tumor_all$n)

#------------ Try now with the copy-number matrix version

# Two matrices: the whole track of copy numbers by patient, and




#------------------------------------------------------------
# Start from CompressiveNMF
SigPrior <- 1e6 * Cosmic_Sigs[, sigs_to_use] + 1
set.seed(42)
MutMatrix <- SigPoisProcess::getTotalMutations(gr_tumor_2kb_std)
outStart <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = SigPrior, K = 0)
CompressiveNMF::plot_SBS_signature(outStart$Signatures)
match_to_RefSigs(outStart$Signatures, Cosmic_Sigs)

R_start <- outStart$Signatures
Theta_start <- outStart$Theta/sum(bin_weight)




# ---- Model 2 - Let signatures vary, starting from CompressiveNMF solution



# ----- Model 3 - fully unsupervised model
set.seed(10)
out_unSup <- SigPoisProcess_mult(gr_Mutations = gr_tumor_2kb_std,
                                     SignalTrack = SignalTrack_2kb_std,
                                     bin_weight = bin_weight,
                                     K = ncol(Sigs),
                                     controls = SigPoisProcess_mult.controls(maxiter = 100),
                                     prior_params = SigPoisProcess_mult.PriorParams(shrinkage = "mu_sigma"))
saveRDS(out_unSup, file = "~/SigPoisProcess/output/SigUnSup_Breast.rds")

CompressiveNMF::plot_SBS_signature(out_unSup$Signatures)

match_to_RefSigs(out_unSup$Signatures, Cosmic_Sigs)
out_unSup$sol$trace
diff(out_unSup$sol$trace)

plot_betas(out_unSup$Betas)













