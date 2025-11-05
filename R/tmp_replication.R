# Tmp replication

data_path <- "~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/"
files <- list.files(path = data_path)
files <- files[grepl("Nhek", files)]
files_full <- file.path(data_path, files)

import_bw <- function(x) rtracklayer::import(x, format="bigWig")
# Import all files and store as a named list
track_list <- setNames(lapply(files_full, import_bw), files)


sum_gr <- track_list$wgEncodeUwRepliSeqNhekSumSignalRep1.bigWig
wave_gr <- track_list$wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig
hits <- findOverlaps(sum_gr, wave_gr)

# Indices of overlapping intervals
query_idx <- queryHits(hits)   # indices in sum_gr
subject_idx <- subjectHits(hits) # indices in wave_gr

# Corresponding GRanges from each
overlap_sum <- sum_gr[query_idx]
overlap_wave <- wave_gr[subject_idx]

plot(overlap_wave[1:2000]$score, overlap_sum[1:2000]$score)

hist(track_list$wgEncodeUwRepliSeqNhekSumSignalRep1.bigWig$score)
hist(track_list$wgEncodeUwRepliSeqNhekS1PctSignalRep1.bigWig$score)

migrate_to_1kb <- function(gr, std_chr = paste0("chr", c(1:22, "X"))){
  gr <- keepSeqlevels(gr, std_chr, pruning.mode = "coarse")
  numvar <- coverage(gr, weight = "score")
  binned_score <- binnedAverage(bins = keepSeqlevels(genome_1kb, std_chr, pruning.mode = "coarse"),
                                numvar = numvar,
                                varname = "score")
  binned_score
}

# Assign to short names for clarity
G1b <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekG1bPctSignalRep1.bigWig)
S1  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS1PctSignalRep1.bigWig)
S2  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS2PctSignalRep1.bigWig)
S3  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS3PctSignalRep1.bigWig)
S4  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS4PctSignalRep1.bigWig)
G2  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekG2PctSignalRep1.bigWig)

plot(G1b[5000:1000]$score, S4[5000:1000]$score)

g1b_vec <- mcols(G1b)$score
s1_vec  <- mcols(S1)$score
s2_vec  <- mcols(S2)$score
s3_vec  <- mcols(S3)$score
s4_vec  <- mcols(S4)$score
g2_vec  <- mcols(G2)$score

old_signal <- ((0.917*g1b_vec) + (0.750*s1_vec) + (0.583*s2_vec) +
  (0.417*s3_vec) + (0.250*s4_vec) + (0*g2_vec))/(0.917 + 0.750 + 0.583 + 0.417 + 0.250 + 0)

new_signal <- ((0*g1b_vec) + (0.250*s1_vec) + (0.417*s2_vec) +
                 (0.583*s3_vec) + (0.750*s4_vec) + (0.917*g2_vec))#/(0.917 + 0.750 + 0.583 + 0.417 + 0.250 + 0)

new_sign2 <- ((0.250*s1_vec) + (0.417*s2_vec))/new_signal * 100
hist(new_signal)
sc1 <- old_signal[1:1000]
sc2 <- new_signal[1:1000]
plot(sc1, sc2)

grX <- keepSeqlevels(genome_1kb, paste0("chr", c(1:22, "X")), pruning.mode = "coarse")
grX$score <- new_sign2

old_signal <- import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig")
plot(old_signal[1:3000]$score, grX$score[1:3000])

grX2 <- grX
grX2$score <- 1 - grX$score/100
hist(grX$score)


blacklist <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/hg19-blacklist.v2.bed")
overlaps <- findOverlaps(grX, blacklist)
grX[queryHits(overlaps)]$score <- 0

hist(log2(grX[-queryHits(overlaps)]$score), breaks = 100)

hist(g1b_vec + s1_vec + s2_vec + s3_vec + s4_vec + g2_vec, breaks = 100)

grX$score <- (s2_vec + s3_vec + s4_vec)/100
hist(grX[-queryHits(overlaps)]$score, breaks = 100)
hist(s1_vec + s2_vec)


hist(g1b_vec)
hist(s1_vec)
hist(s2_vec)
hist(s3_vec)

hist(scale(grX[-queryHits(overlaps)]$score), breaks = 100)
hist(scale(log2(grX[-queryHits(overlaps)]$score/mean(grX[-queryHits(overlaps)]$score) + 0.1)), breaks = 100)

grGC <- import("~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig")
overlapsGC <- findOverlaps(grGC, blacklist)
grGC[queryHits(overlapsGC)]$score <- 0
sc <- grGC[-queryHits(overlapsGC)]$score[grGC[-queryHits(overlapsGC)]$score > 0]
hist(log2(sc/mean(sc)), breaks = 100)
hist(log2(grGC$score))
hist(scale(sc))
ns <- new_signal[new_signal > 0]
hist(scale(ns))

nt <- old_signal$score
nt <- nt[nt > 0]
plot(density(scale(sc)))
lines(density(scale(log2(dnase/mean(dnase)))))
lines(density(scale(pl)))
lines(density(scale(ns)), col = "red")
lines(density(scale(nt)), col = "blue")

pl <- plogis(pmax(g1b_vec/100, 0.98) + 0.01)
hist(pl)
# Let's check at mutations
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")
hist(pl)
table(is.na(pl))

gr_DNase <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Stomach-AdenoCA_male_DNAse-seq.bigWig")
hist(gr_DNase$score, breaks = 100)
which.max(gr_DNase$score)
gr_DNase[121486]
dnase <- pmin(gr_DNase[-queryHits(overlapsGC)]$score, 7) + 0.01
hist(scale(log2(dnase/mean(dnase))), breaks = 100)

transform_gr_scores <- function(gr,
                                trans_log2 = TRUE,
                                reference_cutoff = NULL,
                                scale_score = TRUE,
                                cap_quantile99.99 = TRUE, ...){
  blacklist <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/hg19-blacklist.v2.bed")
  overlaps <- findOverlaps(gr, blacklist)
  gr[queryHits(overlaps)]$score <- 0
  scores <- gr$score
  if(cap_quantile99.99) {
    scores <- pmin(scores, quantile(scores, 0.9999))
  }
  if(trans_log2){
    if(is.null(reference_cutoff)){
      scores[scores>1e-18] <- log2(scores[scores>1e-18]/mean(scores))
    } else {
      scores[scores>1e-18] <- log2(scores[scores>1e-18]/reference_cutoff)
    }
  }
  if(scale_score) {
    scores <- scale(scores)
  }
  gr$score <- scores
  return(gr)
}

merge_scaled_covariate <- function(gr_tumor, gr_cov, score_name = "score"){
  scores <- rep(0, length(gr_tumor))
  bin_weight <- gr_cov$bin_weight
  overlaps <- findOverlaps(gr_tumor, gr_cov)
  scores[queryHits(overlaps)] <- gr_cov[subjectHits(overlaps)]$score_scaled
  gr_tumor$score <- scores
  colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- score_name
  gr_tumor$score_plus <- pmax(0, scores)
  colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- paste0(score_name, "_plus")
  gr_tumor$score_minus <- pmax(0, -scores)
  colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- paste0(score_name, "_minus")
  return(gr_tumor)
}

transform_gr_scores_1kb <- function(gr,
                                    trans_log2 = TRUE,
                                    nugget_log2 = 0.1,
                                    reference_cutoff = NULL,
                                    cap_quantile99.99 = TRUE, ...){
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

  ids <- BinWeight > 0
  # Transform the scores now
  scores_tr <- score[ids]# * BinWeight[BinWeight > 0]/1000 # 1000 because we are dealing with 1kb intervals
  if(cap_quantile99.99) {
    scores_tr <- pmin(scores_tr, quantile(scores_tr, 0.9999))
  }
  if(trans_log2){
    if(is.null(reference_cutoff)){
      scores_tr <- log2((scores_tr)/mean(scores_tr * BinWeight[ids]/1000) + nugget_log2)
    } else {
      scores_tr <- log2(scores_tr/reference_cutoff + nugget_log2)
    }
  }
  # Scale now the score
  scores_scaled <- scale(scores_tr)#(scores_tr - median(scores_tr))/sd(scores_tr)
  # Add positive and negative part of the score
  enrich <- pmax(0, scores_scaled)
  depl <- pmax(0, -scores_scaled)
  # Merge with gr object
  gr$bin_weight <- BinWeight
  gr$score_tr <- 0
  gr$score_scaled <- 0
  gr$score_plus <- 0
  gr$score_minus <- 0

  gr[ids]$score_tr <- scores_tr
  gr[ids]$score_scaled <- scores_scaled
  gr[ids]$score_plus <- pmax(0, scores_scaled)
  gr[ids]$score_minus <- pmax(0, -scores_scaled)
  return(gr)
}

gr_DNase <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Stomach-AdenoCA_male_DNAse-seq.bigWig")
gr_DNase2 <- transform_gr_scores_1kb(gr_DNase, trans_log2 = TRUE)
grGC <- import("~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig")
grGC2 <- transform_gr_scores_1kb(grGC, trans_log2 = FALSE)
grRepli <- import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig")
grRepli2 <- transform_gr_scores_1kb(grRepli, trans_log2 = FALSE)
grH3K9me3 <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K9me3_1kb.bigWig")
grH3K9me3_2 <- transform_gr_scores_1kb(grH3K9me3, trans_log2 = TRUE, reference_cutoff = 1)
grH3K27ac <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K27ac_1kb.bigWig")
grH3K27ac_2 <- transform_gr_scores_1kb(grH3K27ac, trans_log2 = TRUE, reference_cutoff = 1)
grCTCF <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_CTCF_1kb.bigWig")
grCTCF_2 <- transform_gr_scores_1kb(grCTCF, trans_log2 = TRUE, reference_cutoff = 1)

gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")
# Exclude mutations
blacklist <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/hg19-blacklist.v2.bed")
gr_tumor <- gr_tumor[-queryHits(findOverlaps(gr_tumor, blacklist))]
# Merge with covariates
gr_tumor <- merge_scaled_covariate(gr_tumor, gr_DNase2, "DNase")
gr_tumor <- merge_scaled_covariate(gr_tumor, grGC2, "GC")
gr_tumor <- merge_scaled_covariate(gr_tumor, grRepli2, "RepliTime")
gr_tumor <- merge_scaled_covariate(gr_tumor, grH3K9me3_2, "H3K9me3")
gr_tumor <- merge_scaled_covariate(gr_tumor, grCTCF_2, "CTCF")
gr_tumor <- merge_scaled_covariate(gr_tumor, grH3K27ac_2, "H3K27ac")

gr_cov <- grCTCF_2

scores <- rep(0, length(gr_tumor))
overlaps <- findOverlaps(gr_tumor, gr_cov)
scores[queryHits(overlaps)] <- gr_cov[subjectHits(overlaps)]$score_scaled
gr_tumor$score <- scores
colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- score_name


plot(density(gr_tumor$DNase), ylim = c(0, 1))
lines(density(gr_tumor$GC), col = "red")
lines(density(gr_tumor$RepliTime), col = "blue")
lines(density(gr_tumor$H3K9me3), col = "forestgreen")
lines(density(gr_tumor$CTCF), col = "green")
lines(density(gr_tumor$H3K27ac), col = "gray")

plot(density(gr_tumor$DNase_plus[gr_tumor$DNase_plus >0]), ylim = c(0, 1))
lines(density(gr_tumor$DNase_minus[gr_tumor$DNase_minus >0]), col = "red")

plot(density(gr_tumor$GC_plus[gr_tumor$GC_plus >0]), ylim = c(0, 1))
lines(density(gr_tumor$GC_minus[gr_tumor$GC_minus >0]), col = "red")

plot(density(gr_tumor$RepliTime_plus), ylim = c(0, 1))
lines(density(gr_tumor$RepliTime_minus), col = "red")

plot(density(gr_tumor$RepliTime_plus), ylim = c(0, 1))
lines(density(gr_tumor$RepliTime_minus), col = "red")

plot(density(gr_tumor$H3K9me3_plus), ylim = c(0, 1))
lines(density(gr_tumor$H3K9me3_minus), col = "red")

lines(density(gr_tumor$RepliTime), col = "blue")
lines(density(gr_tumor$H3K9me3), col = "forestgreen")
lines(density(gr_tumor$CTCF), col = "green")


hist(pmin(gr_DNase$score, quantile(gr_DNase$score, 0.9999)), breaks = 100)
hist(grGC2$score, breaks = 30)
mean(pmin(gr_DNase$score, 3))

gr <- gr_DNase
sc <- grH3K9me3_2$score
mean(sc)
hist(sc)
mean(sc[sc>1e-18])

gc0 <- grGC[grGC$score == 0]
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, blacklist)
print(seqs)
width(blacklist)
freqs <- Biostrings::alphabetFrequency(seqs)
blacklist[834]

# Remove blacklist, remove regions with NN only


gap_url <- "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz"
gap_df <- read_tsv(gap_url, col_names = FALSE)
library(GenomicRanges)
gap_gr <- GRanges(seqnames = gap_df$X2,
                  ranges   = IRanges(start=gap_df$X3 + 1, end=gap_df$X4))

export(gap_gr, con = "~/SigPoisProcess/data/regions_to_exclude/gaps_hg19.bed")
head(gap_gr)
overlaps <- findOverlaps(gr_tumor, blacklist)

overs <- findOverlaps(grGC, gap_gr)

grH3K4me3 <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K4me3_1kb.bigWig")

gr <- grCTCF
# Remove blacklist
BinSizes <- width(gr)
score <- gr$score
# Remove gaps
overlapGaps <- findOverlaps(gr, gap_gr)
rangesGaps <- gr[queryHits(overlapGaps)]
# Find the number of ACGT in the gap regions, and adjust the bin sizes to remove
# the NNNNN sequences
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rangesGaps)
countACGT <- rowSums((Biostrings::alphabetFrequency(seqs))[, c("A", "C", "G", "T")])
BinSizes[queryHits(overlapGaps)] <- countACGT
# Exlude now the blacklist as well
overlap_blackList <- findOverlaps(gr, blacklist)
pairs <- Pairs(gr[queryHits(overlap_blackList)], blacklist[subjectHits(overlap_blackList)])
grIntersect <- pintersect(pairs)
BinSizes[queryHits(overlap_blackList)] <- pmax(BinSizes[queryHits(overlap_blackList)] - width(grIntersect), 0)
score[BinSizes == 0] <- 0 # Replace it with 0, knowing that we will remove it afterwards




hist(score[BinSizes > 0] * BinSizes[BinSizes > 0]/1000)

sc2 <- score[BinSizes > 0] * BinSizes[BinSizes > 0]/1000
hist(sc2, breaks = 100)
mn <- sum(score[BinSizes > 0] * BinSizes[BinSizes > 0])/sum(BinSizes[BinSizes > 0])
mean()
quantile(sc2)
sc2 <- pmin(sc2, quantile(sc2, 0.9999))
mean(sc2)
sc3 <- scale(log2(sc2 + 0.1))
hist(sc3)
enrich <- pmax(0, sc3)
depl <- pmax(0, -sc3)

enrich2 <- enrich[1:3000]
depl2 <- depl[1:3000]
stats::filter(enrich2, rep(1/3, 3), sides = 2)
plot(stats::filter(enrich2, rep(1/10, 10), sides = 2), type = "l")
lines(stats::filter(depl2, rep(1/10, 10), sides = 2), type = "l", col = "red")

grIntersect



which.max(gr[queryHits(overlapGaps)]$score)
grs <- gr[queryHits(overlapGaps)][156792]
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, grs)







