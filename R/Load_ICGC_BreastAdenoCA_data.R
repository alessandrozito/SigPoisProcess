# This file prepares the dataset for for the application in Section 5

# Load packages
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

create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

############################################################
# Part 1 - functions to build the dataset
############################################################

#----- Function to add a weight to each bin, at 2kb resolution
add_bin_weights <- function(gr){
  # Import blacklist and regions with gaps
  blacklist <- rtracklayer::import("~/SigPoisProcess/data/data_for_application/hg19-blacklist.v2.bed")
  gap_gr <- rtracklayer::import("~/SigPoisProcess/data/data_for_application/gaps_hg19.bed")

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

#----- Function to aggregare a genomic ranges object into a smaller one
bin_gr <- function(gr, genome_aggreg, std_chrs){
  gr <- keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  gr <- binnedAverage(bins = keepSeqlevels(genome_aggreg, std_chrs, pruning.mode = "coarse"),
                      numvar = coverage(gr, weight = "score"),
                      varname = "score")
  return(gr)
}

#----- Function to build the matrix tracking the number of copies in the genome
# of each patient
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

#----- Function to build the matrix of covariates for each region of the genome
build_SignalTrack <- function(tilewidth = 2000, type = "avg"){
  # Genome aggregated
  genome_aggreg <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                              tilewidth = tilewidth,
                              cut.last.tile.in.chrom = TRUE)
  std_chrs <-  paste0("chr", c(1:22, "X"))
  gr_SignalTrack <- add_bin_weights(genome_aggreg)

  # ---- GC content
  print("GC")
  file <- "~/SigPoisProcess/data/data_for_application/gc_content_1kb.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$GC <- gr$score

  # ---- Methylation
  print("Methylation")
  file <- "~/SigPoisProcess/data/data_for_application/Breast-Cancer_Methylation.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$Methyl <- gr$score

  # ---- CTCF + Histone marks programmatically handled
  marks <- c("CTCF", "H3K9me3", "H3K36me3", "H3K27me3", "H3K27ac", "H3K4me1", "H3K4me3")
  for (mark in marks) {
    print(mark)
    if(type == "avg") {
      # Calculate the average between tissue and cell line.
      # Tissue
      file_tissue <- paste0("~/SigPoisProcess/data/data_for_application/Breast-Cancer_tissue_", mark, "_2kb.bigWig")
      gr_tissue <- bin_gr(import(file_tissue), genome_aggreg, std_chrs)
      # Cell line
      file_cell <- paste0("~/SigPoisProcess/data/data_for_application/Breast-Cancer_cell_", mark, "_2kb.bigWig")
      gr_cell <- bin_gr(import(file_cell), genome_aggreg, std_chrs)
      # Average
      gr_SignalTrack$mark <- (gr_tissue$score + gr_cell$score)/2
      colnames(mcols(gr_SignalTrack))[ncol(mcols(gr_SignalTrack))] <- mark
    } else {
      file <- paste0("~/SigPoisProcess/data/data_for_application/Breast-Cancer_", type, "_", mark, "_2kb.bigWig")
      gr <- bin_gr(import(file), genome_aggreg, std_chrs)
      gr_SignalTrack$mark <- gr$score
      colnames(mcols(gr_SignalTrack))[ncol(mcols(gr_SignalTrack))] <- mark
    }
  }

  # ---- Replication Timing
  print("RepliTime")
  file <- "~/SigPoisProcess/data/data_for_application/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig"
  gr <- bin_gr(import(file), genome_aggreg, std_chrs)
  gr_SignalTrack$RepliTime <- gr$score

  # ---- Nucleosome Occupancy
  print("NuclOccup")
  file <- "~/SigPoisProcess/data/data_for_application/GSM920557_hg19_wgEncodeSydhNsomeK562Sig_1kb.bigWig"
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


############################################################
# Part 2 - build the dataset
############################################################

#----- Load the tumor data
gr_tumor <- readRDS("~/SigPoisProcess/data/data_for_application/Breast-AdenoCa_snp.rds.gzip")
# Remove mutations from blacklist
blacklist <- rtracklayer::import("~/SigPoisProcess/data/data_for_application/hg19-blacklist.v2.bed")
overlaps <- findOverlaps(gr_tumor, blacklist)
gr_tumor <- gr_tumor[-queryHits(overlaps)]
# Remove chrY (we are dealing with breast cancer)
gr_tumor <- gr_tumor[seqnames(gr_tumor) != "chrY"]

#----- Build the CopyTrack to record the copy number alteration
# Load the copy number data from ICGC
df_copy <- read_tsv("~/SigPoisProcess/data/data_for_application/20170119_final_consensus_copynumber_donor")
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

# Matrix version of the signal and copytrack
mutMatrix <- getTotalMutations(gr_Mutations)
SignalTrack <- as.matrix(mcols(gr_SignalTrack_2kb_std))[, -1]
CopyTrack <- apply(as.matrix(mcols(gr_CopyTrack))[, -1], 2,
                   function(x) x * gr_CopyTrack$bin_weight)[, colnames(mutMatrix)]

data <- list("gr_Mutations" = gr_Mutations,
             "SignalTrack" = SignalTrack,
             "CopyTrack" = CopyTrack,
             "gr_CopyTrack" = gr_CopyTrack,
             "gr_SignalTrack" = gr_SignalTrack_2kb_std)
#saveRDS(data, file = "~/SigPoisProcess/data/ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip")
rm(list = setdiff(ls(), c("data", names(Filter(is.function, mget(ls(), .GlobalEnv))))), envir = .GlobalEnv)


