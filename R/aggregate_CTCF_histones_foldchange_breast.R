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
genome_2kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                         tilewidth = 2000,
                         cut.last.tile.in.chrom = TRUE)

genome_5kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                          tilewidth = 5000,
                          cut.last.tile.in.chrom = TRUE)


std_chrs <- paste0("chr", c(1:22, "X"))

# Load the files
names_bigWig <- list.files("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/")
names_bigWig <- names_bigWig[grepl("[.]bigWig", names_bigWig)]

# Data for breast
df_names <- read_tsv("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/Files_to_use_breast.tsv")
names_bigWig <- names_bigWig[gsub(".bigWig", "", names_bigWig) %in% df_names$accession]

# Files specific to cancers
for(i in 1:length(names_bigWig)){
  print(i)
  gr <- rtracklayer::import(paste0("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/", names_bigWig[i]))
  gr <- keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  # Calculate coverage
  numvar <- coverage(gr, weight = "score")
  rm(gr)

  # 2kb
  print("2kb")
  binned_score_2kb <- binnedAverage(bins = genome_2kb,
                                    numvar = numvar,
                                    varname = "score")
  export.bw(binned_score_2kb,
            con = paste0("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_2kb/",
                         names_bigWig[i]))

  # 5kb
  print("5kb")
  binned_score_5kb <- binnedAverage(bins = genome_5kb,
                                     numvar = numvar,
                                     varname = "score")
  export.bw(binned_score_5kb,
            con = paste0("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_5kb/",
                         names_bigWig[i]))

}
