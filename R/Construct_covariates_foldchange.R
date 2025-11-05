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
df_histones <- read_tsv("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/Files_to_use.tsv") %>%
  mutate(donor_sex = case_when(grepl("female", biosample_summary) ~ "female",
                               TRUE ~ "male"),
         adult = grepl("adult", biosample_summary)) %>%
  as.data.frame()

mark <- unique(df_histones$target)
tumors <- unique(df_histones$`Cancer Type`)
sex <- unique(df_histones$donor_sex)

dir_out <- "~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/"
dir_1kb <- "~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_1kb/"
for(m in 1:length(mark)) {
  for(s in 1:length(sex)){
    for(t in 1:length(tumors)){
      print(c(mark[m], sex[s], tumors[t]))
      # Find all accession files
      accession_files <- df_histones %>%
        filter(target == mark[m],
               donor_sex == sex[s],
               `Cancer Type`== tumors[t]) %>%
        pull(accession)
      if(length(accession_files) > 0){
        # Load first file
        gr <- rtracklayer::import(paste0(dir_1kb, accession_files[1], ".bigWig"))
        # Iterate across all others
        if(length(accession_files) > 1){
          for(i in 2:length(accession_files)){
            gr_tmp <- rtracklayer::import(paste0(dir_1kb, accession_files[i], ".bigWig"))
            if(length(gr_tmp)!=length(gr)){
              break
            }
            # running mean
            gr$score <- ((i - 1) * gr$score + gr_tmp$score)/i
          }
        }
        nm <- paste0(dir_out, tumors[t], "_", sex[s], "_", mark[m], "_1kb.bigWig")
        export.bw(gr, con = nm)
      }
    }
  }
}

