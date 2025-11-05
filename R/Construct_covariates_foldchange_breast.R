# This file constructs all covariates along the genome by averaging all
# marks for the same sex and same tissue/cell line.
# All files should have been processed already, and made in the same
# format.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

############################################# Histone marks and CTCF
df_histones <- read_tsv("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/Files_to_use_breast.tsv") %>%
  mutate(donor_sex = case_when(grepl("female", biosample_summary) ~ "female",
                               TRUE ~ "male"),
         adult = grepl("adult", biosample_summary),
         class_all =  case_when(grepl("tissue", classification) ~ "tissue",
                            TRUE ~ "cell")) %>%
  as.data.frame()

mark <- unique(df_histones$target)
tumors <- unique(df_histones$`Cancer Type`)
class <- unique(df_histones$class_all)

dir_out <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/"
dir_2kb <- "~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_2kb/"
for(m in 1:length(mark)) {
  for(s in 1:length(class)) {
    for(t in 1:length(tumors)){
      print(c(mark[m], class[s], tumors[t]))
      break
      # Find all accession files
      accession_files <- df_histones %>%
        dplyr::filter(target == mark[m],
                      class_all == class[s],
               `Cancer Type`== tumors[t]) %>%
        pull(accession)
      print(accession_files)
      if(length(accession_files) > 0){
        # Load first file
        gr <- rtracklayer::import(paste0(dir_2kb, accession_files[1], ".bigWig"))
        # Iterate across all others
        if(length(accession_files) > 1){
          for(i in 2:length(accession_files)){
            gr_tmp <- rtracklayer::import(paste0(dir_2kb, accession_files[i], ".bigWig"))
            if(length(gr_tmp)!=length(gr)){
              break
            }
            # running mean
            gr$score <- ((i - 1) * gr$score + gr_tmp$score)/i
          }
        }
        nm <- paste0(dir_out, tumors[t], "_", class[s], "_", mark[m], "_2kb.bigWig")
        export.bw(gr, con = nm)
      }
    }
  }
}

dir_5kb <- "~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_5kb/"
for(m in 1:length(mark)) {
  for(s in 1:length(class)) {
    for(t in 1:length(tumors)){
      print(c(mark[m], class[s], tumors[t]))
      # Find all accession files
      accession_files <- df_histones %>%
        filter(target == mark[m],
               class_all == class[s],
               `Cancer Type`== tumors[t]) %>%
        pull(accession)
      if(length(accession_files) > 0){
        # Load first file
        gr <- rtracklayer::import(paste0(dir_5kb, accession_files[1], ".bigWig"))
        # Iterate across all others
        if(length(accession_files) > 1){
          for(i in 2:length(accession_files)){
            gr_tmp <- rtracklayer::import(paste0(dir_5kb, accession_files[i], ".bigWig"))
            if(length(gr_tmp)!=length(gr)){
              break
            }
            # running mean
            gr$score <- ((i - 1) * gr$score + gr_tmp$score)/i
          }
        }
        nm <- paste0(dir_out, tumors[t], "_", class[s], "_", mark[m], "_5kb.bigWig")
        export.bw(gr, con = nm)
      }
    }
  }
}


gr1 <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Breast-Cancer_cell_H3K4me3_2kb.bigWig")
gr2 <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Breast-Cancer_tissue_H3K4me3_2kb.bigWig")

gr1 <- import("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_2kb/ENCFF216DKX.bigWig")
gr2 <- import("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_2kb/ENCFF094MQS.bigWig")
plot(gr1[1:5000]$score, gr2[1:5000]$score)
cor(gr1$score, gr2$score)
plot(0.5 *(gr1[1:10000]$score + gr2[1:10000]$score), type = "l")
hist(pmin(0.5 *(gr1$score + gr2$score), 4))


quantile(gr1$score)
quantile(gr2$score)
cor(gr1$score, gr2$score)
hist(gr1$score)
hist(gr2$score)


