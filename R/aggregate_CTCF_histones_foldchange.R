# This file aggregates the bigWig file by splitting the genome in
# regions of 1kb/10kb/100kb and calculating the average in each window
# to ensure an equal splitting between regions
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

#--------------------------------------------------------------

# Load the genome and split it
genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)

genome_10kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 10000,
                         cut.last.tile.in.chrom = TRUE)

genome_100kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 100000,
                         cut.last.tile.in.chrom = TRUE)

std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Load the files
names_bigWig <- list.files("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/")
names_bigWig <- names_bigWig[grepl("[.]bigWig", names_bigWig)]

# Files specific to cancers
for(i in 1:length(names_bigWig)){
  print(i)
  gr <- rtracklayer::import(paste0("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/", names_bigWig[i]))
  gr <- keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  # Calculate coverage
  numvar <- coverage(gr, weight = "score")
  rm(gr)

  # 1kb
  print("1kb")
  binned_score_1kb <- binnedAverage(bins = genome_1kb,
                                    numvar = numvar,
                                    varname = "score")
  export.bw(binned_score_1kb,
            con = paste0("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_1kb/", names_bigWig[i]))

  # 10kb
  print("10kb")
  binned_score_10kb <- binnedAverage(bins = genome_10kb,
                                    numvar = numvar,
                                    varname = "score")
  export.bw(binned_score_10kb,
            con = paste0("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_10kb/", names_bigWig[i]))

  # 100kb
  print("100kb")
  binned_score_100kb <- binnedAverage(bins = genome_100kb,
                                    numvar = numvar,
                                    varname = "score")
  export.bw(binned_score_100kb,
            con = paste0("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_100kb/", names_bigWig[i]))

}

# hist(log2(gr$score), breaks = 40)
#
# hist(log2(binned_score_1kb$score + 0.05), breaks = 40)
#
# plot(log2(binned_score_10kb$score + 0.05), type = "l")
# hist(pmax(log2(binned_score_10kb$score), -5), breaks = 30)
#
# plot(log2(binned_score_100kb$score + 0.05), type = "l")
# hist(pmax(log2(binned_score_100kb$score + 0.05), -5), breaks = 30)
#
# (ress4 <- sum(exp(0.2 * log2(gr$score + 0.05)) * width(gr)))
# (ress3 <- sum(exp(0.2 * log2(binned_score_500kb$score + 0.05)) * width(binned_score_500kb)))
# (ress <- sum(exp(0.2 * log2(binned_score_100kb$score + 0.05)) * width(binned_score_100kb)))
# (ress2 <- sum(exp(0.2 * log2(binned_score_10kb$score + 0.05)) * width(binned_score_10kb)))
#
#
# plot(log2(binned_score_100kb$score + 0.05), type = "l")
# hist(pmax(log2(binned_score_100kb$score + 0.05), -5), breaks = 30)
