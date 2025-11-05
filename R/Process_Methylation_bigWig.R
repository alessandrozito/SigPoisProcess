# This file downloads all Methylation covariate as presented by the paper
# Lofer et al. 2023 nature
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

df_methyl <- read_tsv(file = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/suppl/filelist.txt")

tissues <- c("Gastric", "Breast", "Liver", "Prostate", "Keratinocyte", "Epidermal")

df_methyl <- df_methyl %>%
  filter(grepl("bigwig", Name),
         !grepl("hg38", Name),
         grepl(paste(tissues, collapse="|"), Name)) %>%
  filter(!grepl("Prostate-Smooth-Muscle", Name),
         !grepl("Liver-Macrophages", Name)) %>%
  as.data.frame()

# Call now the sh file to download everything
genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

names_bigWig <- list.files("~/SigPoisProcess/data/ENCODE_pvalues/Methylation_data/Geo_methylation_bigWig/")
names_bigWig <- names_bigWig[grepl("[.]bigwig", names_bigWig)]

for(i in 1:length(names_bigWig)){
  print(i)
  gr <- rtracklayer::import(paste0("~/SigPoisProcess/data/ENCODE_pvalues/Methylation_data/Geo_methylation_bigWig/", names_bigWig[i]))
  std_chrs_tmp <- std_chrs[std_chrs %in% levels(seqnames(gr))]
  gr <- keepSeqlevels(gr, std_chrs_tmp, pruning.mode = "coarse")
  binned_score <- binnedAverage(bins = keepSeqlevels(genome_1kb, std_chrs_tmp, pruning.mode = "coarse"),
                                numvar = coverage(gr, weight = "score"),
                                varname = "score")
  #rm(numvar)
  rm(gr)
  export.bw(binned_score,
            con = paste0("~/SigPoisProcess/data/ENCODE_pvalues/Methylation_data/Geo_methylation_bigWig/", names_bigWig[i]))
}





