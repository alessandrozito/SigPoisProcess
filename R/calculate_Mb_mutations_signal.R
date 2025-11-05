# This script merges covariates at the 1kb resolution with the patients.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

calculate_Mb_mutations_signal <- function(TumorData){
  # Load the genome 1MegaBase
  genome <- BSgenome.Hsapiens.UCSC.hg19
  chrom_lengths <- seqlengths(genome)[1:24]
  hg19_1Mb <- tileGenome(chrom_lengths, tilewidth = 1e6,
                         cut.last.tile.in.chrom = TRUE)
  # Overlap with Tumor
  overlaps_1Mb <- findOverlaps(TumorData$gr_tumor, hg19_1Mb)

  # Step 1 - Tensor of mutations at the megabase scale
  df_tumor <- as.data.frame(TumorData$gr_tumor)
  df_tumor$seqnames <- factor(df_tumor$seqnames, levels = names(chrom_lengths))
  df_tumor$region <- subjectHits(overlaps_1Mb)
  df_tumor$region <- factor(df_tumor$region, levels = 1:length(hg19_1Mb))
  df_tumor$sample <- as.factor(df_tumor$sample)

  # Calculate mutation rate at Mb scale
  MutationsMb <- as.array(table(df_tumor$region, df_tumor$sample, df_tumor$channel))

  covariates <- c("baseline", colnames(TumorData$df_files[, -c(1,2)]))
  patients <- colnames(MutationsMb)
  regions <- 1:length(hg19_1Mb)

  # Step 2 - Summing signal at the megabase scale

  # Initialize the quantity
  CovariatesMb <- array(0, dim = c(length(hg19_1Mb), ncol(MutationsMb), length(covariates)),
                        dimnames = list(regions, patients, covariates))
  # Initialize the quantity
  for(p in 1:length(covariates)){
    print(covariates[p])
    if(covariates[p] != "baseline"){
      files <- TumorData$df_files[, colnames(TumorData$df_file) == covariates[p]]
      unique_files <- unique(files)
      for(f in 1:length(unique_files)){
        file <- unique_files[f]
        gr <- rtracklayer::import(file)
        overlaps <- findOverlaps(gr, hg19_1Mb)
        gr_df <- as.data.frame(gr)
        gr_df$region <-  subjectHits(overlaps)
        df_score <- gr_df %>%
          group_by(region) %>%
          summarise(sumScore = sum(score * width))
        score_sc <- rep(0, length(regions))
        score_sc[df_score$region] <- df_score$sumScore
        # Append to all patients
        patients_sub <- TumorData$df_files$sample[files == file]
        for(j in 1:length(patients_sub)){
          CovariatesMb[, patients_sub[j], covariates[p]] <- score_sc
        }
      }
    } else {
      for(j in 1:length(patients)){
        CovariatesMb[, j, covariates[p]] <- width(hg19_1Mb)
      }
    }
  }
  return(list(MutationsMb = MutationsMb,
              CovariatesMb = CovariatesMb))
}

# Stomach-AdenoCA
print("Stomach-AdenoCA")
TumorData <- readRDS("data/ICGC_snp_merged/Stomach-AdenoCA_merged_1kb.rds.gzip")
Mutations_Covariates_Mb <- calculate_Mb_mutations_signal(TumorData)
saveRDS(Mutations_Covariates_Mb, file = "data/ICGC_snp_merged/Stomach-AdenoCA_Mutations_Covariates_Mb.rds.gzip", compress = "gzip")

# Liver-HCC
print("Liver-HCC")
TumorData <- readRDS("data/ICGC_snp_merged/Liver-HCC_merged_1kb.rds.gzip")
Mutations_Covariates_Mb <- calculate_Mb_mutations_signal(TumorData)
saveRDS(Mutations_Covariates_Mb, file = "data/ICGC_snp_merged/Liver-HCC_Mutations_Covariates_Mb.rds.gzip", compress = "gzip")

# Breast-AdenoCa
print("Breast-AdenoCa")
TumorData <- readRDS("data/ICGC_snp_merged/Breast-Cancer_merged_1kb.rds.gzip")
Mutations_Covariates_Mb <- calculate_Mb_mutations_signal(TumorData)
saveRDS(Mutations_Covariates_Mb, file = "data/ICGC_snp_merged/Breast-Cancer_Mutations_Covariates_Mb.rds.gzip", compress = "gzip")

# Prost-AdenoCa
print("Prost-AdenoCa")
TumorData <- readRDS("data/ICGC_snp_merged/Prost-AdenoCA_merged_1kb.rds.gzip")
Mutations_Covariates_Mb <- calculate_Mb_mutations_signal(TumorData)
saveRDS(Mutations_Covariates_Mb, file = "data/ICGC_snp_merged/Prost-AdenoCA_Mutations_Covariates_Mb.rds.gzip", compress = "gzip")

# Skin-melanoma
print("Skin-melanoma-AdenoCa")
TumorData <- readRDS("data/ICGC_snp_merged/Skin-Melanoma_merged_1kb.rds.gzip")
Mutations_Covariates_Mb <- calculate_Mb_mutations_signal(TumorData)
saveRDS(Mutations_Covariates_Mb, file = "data/ICGC_snp_merged/Skin-Melanoma_Mutations_Covariates_Mb.rds.gzip", compress = "gzip")






