# This file dowloads the data from Replication timing and aggregates them under the usual
# 1kb window.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

df_repli <- read_tsv("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/Replication_time_files.tsv")

df_repli <- df_repli %>%
  filter(View == "WaveSignal",
         `Cell line` %in% c("MCF-7", "HepG2", "NHEK")) %>%
  as.data.frame()

genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Load the files
names_bigWig <- list.files("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/")
names_bigWig <- names_bigWig[grepl("[.]bigWig", names_bigWig)]

# Replication timing
for(i in 1:length(names_bigWig)){
  print(i)
  gr <- rtracklayer::import(paste0("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/", names_bigWig[i]))
  std_chrs_tmp <- std_chrs[std_chrs %in% levels(seqnames(gr))]
  gr <- keepSeqlevels(gr, std_chrs_tmp, pruning.mode = "coarse")
  numvar <- coverage(gr, weight = "score")
  binned_score <- binnedAverage(bins = keepSeqlevels(genome_1kb, std_chrs_tmp, pruning.mode = "coarse"),
                                numvar = numvar,
                                varname = "score")
  #rm(numvar)
  rm(gr)
  export.bw(binned_score,
            con = paste0("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/", names_bigWig[i]))
}



