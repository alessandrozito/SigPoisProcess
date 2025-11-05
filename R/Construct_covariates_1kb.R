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
df_histones <- read_tsv("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bigWigfiles_to_use.tsv")[-29, ] %>%
  mutate(donor_sex = case_when(grepl("female", biosample_summary) ~ "female",
                               TRUE ~ "male"),
         adult = grepl("adult", biosample_summary)) %>%
  as.data.frame()

mark <- unique(df_histones$target)
tumors <- unique(df_histones$`Cancer Type`)
sex <- unique(df_histones$donor_sex)

dir_out <- "~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/"
dir_1kb <- "~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/CTCF_histones_1kb/"
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
        nm <- paste0(dir_out, tumors[t], "_", sex[s], "_", mark[m], ".bigWig")
        export.bw(gr, con = nm)
      }
    }
  }
}


############################################# POLR2A
df_polr2a <- read_tsv("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/POLR2A_bigWig/POLR2A_accession_files.tsv", skip = 1) %>%
  filter(`Preferred Default`,
         !is.na(`Simple biosample summary`)) %>%
  dplyr::select(Organ, Accession, `Simple biosample summary`)%>%
  # Separate the bio sample summary into parts
  mutate(
    `Cancer Type` = case_when(grepl("stomach", Organ)  ~ "Stomach-AdenoCA",
                           grepl("prostate", Organ)  ~ "Prost-AdenoCA",
                           grepl("breast", Organ)  ~ "Breast-Cancer",
                           grepl("skin", Organ)  ~ "Skin-Melanoma",
                           grepl("liver", Organ)  ~ "Liver-HCC"),
    donor_sex = if_else(str_detect(`Simple biosample summary`, "female"), "female", "male")) %>%
  as.data.frame()

mark <- "POLR2A"
tumors <- unique(df_polr2a$`Cancer Type`)
sex <- unique(df_polr2a$donor_sex)

dir_out <- "~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/"
dir_1kb <- "~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/POLR2A_1kb/"

for(m in 1:length(mark)) {
  for(s in 1:length(sex)){
    for(t in 1:length(tumors)){
      print(c(mark[m], sex[s], tumors[t]))
      # Find all accession files
      accession_files <- df_polr2a %>%
        filter(donor_sex == sex[s],
               `Cancer Type`== tumors[t]) %>%
        pull(Accession)
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
        nm <- paste0(dir_out, tumors[t], "_", sex[s], "_", mark[m], ".bigWig")
        export.bw(gr, con = nm)
      }
    }
  }
}


#-------------------------------------------- DNAse-seq
df_DNAse <- read_tsv("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/DNAse-seq_BigWig.tsv") %>%
  # Separate the bio sample summary into parts
  mutate(
    `Cancer Type` = case_when(grepl("stomach", organ)  ~ "Stomach-AdenoCA",
                              grepl("prostate", organ)  ~ "Prost-AdenoCA",
                              grepl("breast", organ)  ~ "Breast-Cancer",
                              grepl("skin", organ)  ~ "Skin-Melanoma",
                              grepl("liver", organ)  ~ "Liver-HCC"),
    donor_sex = if_else(str_detect(`biosample_summary`, "female"), "female", "male")) %>%
  as.data.frame()

mark <- "DNAse-seq"
tumors <- unique(df_DNAse$`Cancer Type`)
sex <- unique(df_DNAse$donor_sex)

dir_out <- "~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/"
dir_1kb <- "~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_1kb/"

for(m in 1:length(mark)) {
  for(s in 1:length(sex)){
    for(t in 1:length(tumors)){
      print(c(mark[m], sex[s], tumors[t]))
      # Find all accession files
      accession_files <- df_DNAse %>%
        filter(donor_sex == sex[s],
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
        nm <- paste0(dir_out, tumors[t], "_", sex[s], "_", mark[m], ".bigWig")
        export.bw(gr, con = nm)
      }
    }
  }
}


#-------------------------------------------- Methylation
df_methyl <- read_tsv(file = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/suppl/filelist.txt")
tissues <- c("Gastric", "Breast", "Liver", "Prostate", "Keratinocyte", "Epidermal")
df_methyl <- df_methyl %>%
  filter(grepl("bigwig", Name),
         !grepl("hg38", Name),
         grepl(paste(tissues, collapse="|"), Name)) %>%
  filter(!grepl("Prostate-Smooth-Muscle", Name),
         !grepl("Liver-Macrophages", Name)) %>%
  mutate(`Cancer Type` = case_when(grepl("Gastric", Name)  ~ "Stomach-AdenoCA",
                                   grepl("Prostate", Name)  ~ "Prost-AdenoCA",
                                   grepl("Breast", Name)  ~ "Breast-Cancer",
                                   grepl("Keratinocyte", Name) | grepl("Epidermal", Name) ~ "Skin-Melanoma",
                                   grepl("Liver", Name)  ~ "Liver-HCC")) %>%
  as.data.frame()

tumors <- unique(df_methyl$`Cancer Type`)

dir_out <- "~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/"
dir_1kb <- "~/SigPoisProcess/data/ENCODE_pvalues/Methylation_data/Geo_methylation_bigWig/"
for(t in 1:length(tumors)){
  print(tumors[t])
  # Find all accession files
  accession_files <- df_methyl %>%
    filter(`Cancer Type`== tumors[t]) %>%
    pull(Name)
  if(length(accession_files) > 0){
    # Load first file
    gr <- rtracklayer::import(paste0(dir_1kb, accession_files[1]))
    # Iterate across all others
    if(length(accession_files) > 1){
      for(i in 2:length(accession_files)){
        gr_tmp <- rtracklayer::import(paste0(dir_1kb, accession_files[i]))
        if(length(gr_tmp)!=length(gr)){
          break
        }
        # running mean
        gr$score <- ((i - 1) * gr$score + gr_tmp$score)/i
      }
    }
    nm <- paste0(dir_out, tumors[t], "_Methylation.bigWig")
    export.bw(gr, con = nm)
  }
}


#-------------------------------------------- Baseline mutation rate from whole ICGC

# Load the genome and split it
genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X", "Y"))
df_clinical <- read_tsv("data/ICGC_rawdata/pcawg_donor_clinical_August2016_v9.tsv")

dir_mutations <- "data/ICGC_snp_mutations/"
files <- list.files(dir_mutations)

df_all <- as.data.frame(expand.grid(region = 1:length(genome_1kb), donor_sex = c("male", "female")))

for(f in files){
  print(f)
  gr <- readRDS(paste0(dir_mutations, f))
  df_tumor <- as.data.frame(gr) %>%
    left_join(df_clinical %>% dplyr::select(icgc_donor_id, donor_sex),
              by = c("sample" = "icgc_donor_id"))
  overlaps <- findOverlaps(gr, genome_1kb)
  df_tumor$region <- subjectHits(overlaps)

  # Find aggregate mutation rate
  muts <- df_tumor %>%
    group_by(region, donor_sex) %>%
    summarise(n = n()) %>%
    left_join(df_tumor %>%
                dplyr::select(sample, donor_sex) %>%
                distinct() %>%
                group_by(donor_sex) %>%
                summarize(nsex = n()), by = "donor_sex") %>%
    mutate(n = n/nsex) %>%
    dplyr::select(region, donor_sex, n)
  colnames(muts)[3] <- f
  df_all <- df_all %>%
    left_join(muts, by = c("region", "donor_sex"))
}

# Replace NA with 0
df_all[is.na(df_all)] <- 0
df_all$totals <- rowSums(df_all[, -c(1,2)])


gr_male <- genome_1kb
gr_male$score <- df_all$totals[df_all$donor_sex == "male"]
nm <- paste0("data/ENCODE_pvalues/Cancer_covariates/Baseline_male.bigWig")
export.bw(gr_male, con = nm)

gr_female <- genome_1kb
gr_female$score <- df_all$totals[df_all$donor_sex == "female"]
nm <- paste0("data/ENCODE_pvalues/Cancer_covariates/Baseline_female.bigWig")
export.bw(gr_female, con = nm)

















