library(rvest)
# Initialization of the SigPoisProcess method using standard NMF +
# Poisson regression with identity link.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
#library(SigPoisProcess)
library(rtracklayer)
library(foreach)
library(doParallel)
devtools::document()


make_scored_genome <- function(gr, genome="hg19", chroms = c(paste0("chr", 1:22), "chrX", "chrY")) {
  # Get chromosome info
  si <- Seqinfo(genome=genome)
  si <- si[chroms] # Only standard chroms

  # Prune & update gr to match chrom set, set seqlengths
  seqlevels(gr, pruning.mode="coarse") <- chroms
  seqlengths(gr) <- seqlengths(si)

  # FORCE UNSTRANDED
  strand(gr) <- "*"

  # Generate full coverage for each chrom
  cov <- coverage(gr, width = seqlengths(si)[chroms])

  # Convert coverage RleList to GRanges
  cov_gr <- as(cov, "GRanges")
  seqinfo(cov_gr) <- si

  # Keep only intervals with coverage 0 or 1+ (avoid negatives, NAs...)
  cov_gr <- cov_gr[as.character(seqnames(cov_gr)) %in% chroms] # only wanted chroms
  cov_gr <- cov_gr[start(cov_gr) > 0 & end(cov_gr) <= seqlengths(si)[as.character(seqnames(cov_gr))]]
  cov_gr$score <- ifelse(cov_gr$score > 0, 1L + 1e-18, 1e-18)

  cov_gr
}

# Average mutation rate at Kb
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:24]
gr_tumor <- readRDS(file = "~/SigPoisProcess/data/ICGC_snp_mutations/Liver-HCC_snp.rds.gzip")
hg19_1kb <- tileGenome(chrom_lengths, tilewidth = 1000, cut.last.tile.in.chrom = TRUE)

df_Kb <- as.data.frame(gr_tumor) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor, hg19_1kb))) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  tidyr::complete(region = 1:length(hg19_1kb), fill = list(n = 0))

plot(table(df_Kb$n))

# Open files
# CTCF
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_female_CTCF.bigWig")
df_Kb$CTCF <- gr$score
# DNase
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_female_DNAse-seq.bigWig")
df_Kb$DNase <- gr$score
# H3K9me3
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_female_H3K9me3.bigWig")
df_Kb$H3K9me3 <- gr$score
# H3K27me3
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_male_H3K27me3.bigWig")
df_Kb$H3K27me3 <- gr$score
# H3K4me3
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_male_H3K4me3.bigWig")
df_Kb$H3K4me3 <- gr$score
# H3K27ac
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_male_H3K27ac.bigWig")
df_Kb$H3K27ac <- gr$score
# Methyl
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_Methylation.bigWig")
df_Kb$Methyl <- gr$score

library(nnls)
sol <- nnls(A = as.matrix(df_Kb[, -c(1,2)]), b = df_Kb$n)
sol$residuals
fit <- sol$fitted[1:1000]
n <- df_Kb$n[1:1000]
plot(fit-n)
sum(df_Kb$n)
sum(sol$fitted)

plot(fit, n)
fit_pois <- glm(n ~. , family = quasipoisson(), data = df_Kb %>% dplyr::select(-region))
summary(fit_pois)
plot(fit_pois)


df_Kb_scaled <- df_Kb %>%
  mutate(across(
    .cols = where(is.numeric) & !c(n, region),  # Exclude non-covariates
    .fns = ~ as.numeric(scale(.x))))

fit_poi_scaled <- glm(n ~. , family = quasipoisson(), data = df_Kb_scaled %>% dplyr::select(-region))
summary(fit_poi_scaled)

plot(fit_poi_scaled)

luad <- TCGAmutations::tcga_load(study = "BRCA")
table(luad@maf.silent$Tumor_Sample_Barcode)



