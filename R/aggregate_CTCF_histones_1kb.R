# This file aggregates the bigWig file by splitting the genome in
# regions of 1kb and calculating the average in each window
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
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Load the files
names_bigWig <- list.files("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bigWig/")
names_bigWig <- names_bigWig[grepl("[.]bigWig", names_bigWig)]


# Files specific to cancers
for(i in 1:length(names_bigWig)){
  print(i)
  gr <- rtracklayer::import(paste0("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bigWig/", names_bigWig[i]))
  gr <- keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  #numvar <- coverage(gr, weight = "score")
  binned_score <- binnedAverage(bins = genome_1kb,
                                numvar = coverage(gr, weight = "score"),
                                varname = "score")
  #rm(numvar)
  rm(gr)
  export.bw(binned_score,
            con = paste0("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/CTCF_histones_1kb/", names_bigWig[i]))
}




