# This script merges covariates at the 1kb resolution with the patients.
# Stomach-AdenoCA
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

source("~/SigPoisProcess_analysis/utils.R")

#------------------------------------------------------------------------
transform_scores <- function(scores){
  #sc <- scores + 1e-12
  #norms <- (sc - mean(sc))/sd(sc)
  #-log10(pnorm(1 - norms) + 1e-12)
  #norms + 1 + 1e-12
  (scores + 1e-12)/mean(scores)
}

#transform_scores <- function(scores, ...){
  #sc <- scores + 1e-12
  #norms <- (sc - mean(sc))/sd(sc)
  #-log10(pnorm(1 - norms) + 1e-12)
  #norms + 1 + 1e-12
 # sc_score <- (scores - mean(scores))/(2*sd(scores)) + 1
#  print(table(sc_score < 0))
#  sc_score[sc_score < 1e-12] <- 1e-12
#  c(sc_score)
  #sc_score <- exp(sc_score)
  #sc_score[]
  #(scores + 1e-12)/mean(scores)
#}


transform_gc <- function(scores, q = 0.265, square = FALSE){
  scores <- scores - q
  scores[scores <= 0] <- 0
  (scores + 1e-12)/mean(scores)
}


transform_gc <- function(scores, ...){
  q <- quantile(scores[scores>0], 0.01)
  scores[scores>0] <- scores[scores>0] - q
  scores[scores<1e-12] <- 1e-12
  sc_score <- (scores - mean(scores))/(2*sd(scores)) + 1
  print(table(sc_score < 0))
  sc_score[sc_score < 1e-12] <- 1e-12
  c(sc_score)
}

transform_scores <- function(scores, ...){
  q <- quantile(scores[scores>0], 0.01)
  scores[scores>0] <- scores[scores>0] - q
  scores[scores<1e-12] <- 1e-12
  scores/mean(scores)
}

transform_gc <- transform_scores

merge_with_patients <- function(gr_tumor, tumor){
  files_subset <- gr_files[grepl(tumor, gr_files)]

  donor_sex <- c("male", "female")
  mark <- unique(sub(".*_(.+)\\.bigWig$", "\\1", files_subset))

  patients <- gr_tumor$sample
  is_male <- which(patients %in% patiens_male)
  is_female <- which(patients %in% patiens_female)

  all(unique(patients) %in% df_clinical$icgc_donor_id)

  df_areas <- data.frame(sample = unique(patients)) %>%
    left_join(df_clinical %>% dplyr::select(icgc_donor_id, donor_sex),
              by = c("sample" = "icgc_donor_id"))


  ###########################################
  # Merge histone marks, POL2RA and CTCF

  for(i in 1:length(mark)){
    print(mark[i])
    if(mark[i] == "Methylation"){
      file_male <-  files_subset[grepl("Methylation", files_subset)]
      file_female <-  files_subset[grepl("Methylation", files_subset)]
    } else {
      file_male <-  files_subset[grepl("_male_", files_subset) & grepl(mark[i], files_subset)]
      file_female <-  files_subset[grepl("_female_", files_subset) & grepl(mark[i], files_subset)]
      print(file_male)
      print(file_female)
    }
    scores <- rep(1e-12, length(gr_tumor))
    # Load both files
    gr_male <-  rtracklayer::import(paste0(dir_covariates, file_male))
    gr_female <-  rtracklayer::import(paste0(dir_covariates, file_female))
    gr_female <- keepSeqlevels(gr_female, std_chrs[-24], pruning.mode = "coarse")

    # Merge them with male
    gr_male$score <- transform_scores(gr_male$score)
    scores_male <- rep(1e-12, length(gr_tumor))
    overlaps_male <- findOverlaps(gr_tumor, gr_male)
    scores_male[queryHits(overlaps_male)] <- gr_male[subjectHits(overlaps_male)]$score

    # Merge them with female
    gr_female$score <- transform_scores(gr_female$score)
    scores_female <- rep(1e-12, length(gr_tumor))
    overlaps_female <- findOverlaps(gr_tumor, gr_female)
    scores_female[queryHits(overlaps_female)] <- gr_female[subjectHits(overlaps_female)]$score

    # Find which is male and which is female
    scores[is_female] <- scores_female[is_female]
    scores[is_male] <- scores_male[is_male]

    # Append to new covariate
    gr_tumor$score <- scores
    colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- mark[i]

    # Calculate areas
    df_areas <- df_areas %>%
      left_join(data.frame("donor_sex" = c("female", "male"),
                           "area" = c(sum(width(gr_male) * gr_male$score),
                                      sum(width(gr_female) * gr_female$score))),
                by = "donor_sex")
    colnames(df_areas)[ncol(df_areas)] <- mark[i]

  }

  ###########################################
  print("Nucl Occupancy")
  # Merge Nucleosome occupancy (this is the only file related to nucleosome occupancy)
  gr_Nucl <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/NuclOccupancy/GSM920557_hg19_wgEncodeSydhNsomeK562Sig_1kb.bigWig")
  gr_Nucl$score <- transform_scores(gr_Nucl$score)
  overlaps <- findOverlaps(gr_tumor, gr_Nucl)
  scores <- rep(1e-12, length(gr_tumor))
  scores[queryHits(overlaps)] <- gr_Nucl[subjectHits(overlaps)]$score
  gr_tumor$score <- scores
  colnames(mcols(gr_tumor))[length(colnames(mcols(gr_tumor)))] <- "NuclOccup"
  df_areas$NuclOccup <- sum(width(gr_Nucl) * gr_Nucl$score)

  ###########################################
  print("Replication timing")
  # Merge Replication timing (Stomach is NHEK) only available for females...
  if(grepl("Breast", tumor)){
    tmp <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig")
  } else if (grepl("Liver", tumor)){
    tmp <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig")
  } else {
    tmp <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig")
  }
  tmp$score <- transform_scores(tmp$score)
  overlaps <- findOverlaps(gr_tumor, tmp)
  scores <- rep(1e-12, length(gr_tumor))
  scores[queryHits(overlaps)] <- tmp[subjectHits(overlaps)]$score
  gr_tumor$score <- scores
  colnames(mcols(gr_tumor))[length(colnames(mcols(gr_tumor)))] <- "RepliTime"
  df_areas$RepliTime <- sum(width(tmp) * tmp$score)

  ###########################################
  print("GC content")
  # Merge GC content
  gr_male <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig")
  gr_female <- keepSeqlevels(gr_male, std_chrs[-24], pruning.mode = "coarse")
  scores <- rep(1e-12, length(gr_tumor))

  # Merge them with male
  gr_male$score <- transform_gc(gr_male$score, q = 0.265)
  scores_male <- rep(1e-12, length(gr_tumor))
  overlaps_male <- findOverlaps(gr_tumor, gr_male)
  scores_male[queryHits(overlaps_male)] <- gr_male[subjectHits(overlaps_male)]$score

  # Merge them with female
  gr_female$score <- transform_gc(gr_female$score, q = 0.265)
  scores_female <- rep(1e-12, length(gr_tumor))
  overlaps_female <- findOverlaps(gr_tumor, gr_female)
  scores_female[queryHits(overlaps_female)] <- gr_female[subjectHits(overlaps_female)]$score

  # Find which is male and which is female
  scores[is_female] <- scores_female[is_female]
  scores[is_male] <- scores_male[is_male]

  # Append to new covariate
  gr_tumor$score <- scores
  colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- "GC"

  # Calculate areas
  df_areas <- df_areas %>%
    left_join(data.frame("donor_sex" = c("female", "male"),
                         "area" = c(sum(width(gr_male) * gr_male$score),
                                    sum(width(gr_female) * gr_female$score))),
              by = "donor_sex")
  colnames(df_areas)[ncol(df_areas)] <- "GC"

  ###########################################
  print("GC content square")
  # Merge GCsq content
  gr_male <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig")
  gr_female <- keepSeqlevels(gr_male, std_chrs[-24], pruning.mode = "coarse")
  scores <- rep(1e-12, length(gr_tumor))

  # Merge them with male
  gr_male$score <- transform_gc(gr_male$score^2, q = 0.1)
  scores_male <- rep(1e-12, length(gr_tumor))
  overlaps_male <- findOverlaps(gr_tumor, gr_male)
  scores_male[queryHits(overlaps_male)] <- gr_male[subjectHits(overlaps_male)]$score

  # Merge them with female
  gr_female$score <- transform_gc(gr_female$score^2, q = 0.1)
  scores_female <- rep(1e-12, length(gr_tumor))
  overlaps_female <- findOverlaps(gr_tumor, gr_female)
  scores_female[queryHits(overlaps_female)] <- gr_female[subjectHits(overlaps_female)]$score

  # Find which is male and which is female
  scores[is_female] <- scores_female[is_female]
  scores[is_male] <- scores_male[is_male]

  # Append to new covariate
  gr_tumor$score <- scores
  colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- "GCsq"

  # Calculate areas
  df_areas <- df_areas %>%
    left_join(data.frame("donor_sex" = c("female", "male"),
                         "area" = c(sum(width(gr_male) * gr_male$score),
                                    sum(width(gr_female) * gr_female$score))),
              by = "donor_sex")
  colnames(df_areas)[ncol(df_areas)] <- "GCsq"

  ########################################### SAVE FILE!
  TumorData <- list(gr_tumor = gr_tumor,
                    df_areas = df_areas)
  saveRDS(TumorData, paste0("data/ICGC_snp_merged/", tumor, "_merged_1kb.rds.gzip"), compress = "gzip")
}
#------------------------------------------------------------------------

# Chromosome names
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Load clinical data
df_clinical <- read_tsv("data/ICGC_rawdata/pcawg_donor_clinical_August2016_v9.tsv")
patiens_male <- df_clinical$icgc_donor_id[df_clinical$donor_sex == "male"]
patiens_female <- df_clinical$icgc_donor_id[df_clinical$donor_sex == "female"]

# Load covariate names
dir_covariates <- "data/ENCODE_pvalues/Cancer_covariates/"
gr_files <- list.files(dir_covariates)

##########################################
# Stomach-AdenoCA
##########################################
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")
tumor <- "Stomach-AdenoCA"
merge_with_patients(gr_tumor, tumor)

##########################################
# Breast-AdenoCa
##########################################
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Breast-AdenoCa_snp.rds.gzip")
tumor <- "Breast-Cancer"
merge_with_patients(gr_tumor, tumor)

##########################################
# Skin-Melanoma
##########################################
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Skin-Melanoma_snp.rds.gzip")
tumor <- "Skin-Melanoma"
merge_with_patients(gr_tumor, tumor)

##########################################
# Prost-AdenoCA
##########################################
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Prost-AdenoCA_snp.rds.gzip")
tumor <- "Prost-AdenoCA"
merge_with_patients(gr_tumor, tumor)

##########################################
# Liver-HCC
##########################################
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Liver-HCC_snp.rds.gzip")
tumor <- "Liver-HCC"
merge_with_patients(gr_tumor, tumor)


LiverTumor <- readRDS("data/ICGC_snp_merged/Liver-HCC_merged_1kb.rds.gzip")
# Merge also each patient with gene expression. These are available for livehcc only.
gr_tumor <- LiverTumor$gr_tumor
df_areas <- LiverTumor$df_areas
list_gene_expr <- list.files("data/ICGC_GeneExpr_bigWig/")
patients <- gsub(".bigWig", "", list_gene_expr)

patients <- unique(gr_tumor$sample)[unique(gr_tumor$sample) %in% patients]

gr_tumor_expr <- NULL
df_areas_expr <- data.frame()
for(i in 1:length(patients)){
  print(i)
  # Look at patient i
  gr_tmp <- gr_tumor[gr_tumor$sample == patients[i]]
  scores <- rep(1e-12, length(gr_tmp))
  # Load gene expression and find overlaps
  tmp_exp <- rtracklayer::import(paste0("data/ICGC_GeneExpr_bigWig/", patients[i], ".bigWig"))
  tmp_exp$score <- transform_scores(tmp_exp$score)
  overlaps <- findOverlaps(gr_tmp, tmp_exp)
  # Merge the score
  scores[queryHits(overlaps)] <- tmp_exp[subjectHits(overlaps)]$score
  gr_tmp$score <- scores
  colnames(mcols(gr_tmp))[ncol(mcols(gr_tmp))] <- "GeneExpr"
  gr_tumor_expr <- append(gr_tumor_expr, gr_tmp)
  # Add area
  area <- sum(width(tmp_exp) * tmp_exp$score)
  df_areas_expr <- rbind(df_areas_expr,
                         data.frame("sample" = patients[i],
                                    "GeneExpr" = area))
}

# Merge with old areas
df_areas_expr <- df_areas %>%
  left_join(df_areas_expr, by = "sample") %>%
  drop_na()

# Merge with old areas
TumorData <- list(gr_tumor = gr_tumor_expr,
                  df_areas = df_areas_expr)
saveRDS(TumorData, paste0("data/ICGC_snp_merged/Liver-HCC_GeneExpr_merged_1kb.rds.gzip"), compress = "gzip")



#
# ###########################################
# # Run the method
# # Load the genome
# library(SigPoisProcess)
# X <- as.matrix(mcols(gr_tumor)[, -c(1,2,3)])
# X_totals <- as.matrix(df_areas[, -c(1,2)])
# rownames(X_totals) <- df_areas$sample
# #X_totals2 <- cbind("baseline" = sum(1 * chrom_lengths), X_totals)
# Mutations <- as.data.frame(mcols(gr_tumor)[, c(2,3)])
# Mutations$sample <- factor(Mutations$sample, levels = rownames(X_totals))
# Mutations$channel <- as.factor(Mutations$channel)
#
# a <- 3
# alpha <- 1.01
# K <- 10
# set.seed(42)
# out_Proc_MAP <- SigPoisProcess(Mutations = Mutations,
#                                 X = X,
#                                 X_total = X_totals,
#                                 method = "map",
#                                 prior_params = list(a = a, alpha = alpha),
#                                 K = K,
#                                 controls = list(tol = 1e-6,
#                                                 maxiter = 1000,
#                                                 merge_move = TRUE))
#
# round(out_Proc_MAP$Mu, 5)
# CompressiveNMF::plot_SBS_signature(out_Proc_MAP$Signatures)
#
# ###########################################
# # Run the method
# # Load the genome
# genome <- BSgenome.Hsapiens.UCSC.hg19
# std_chrs <- paste0("chr", c(1:22, "X", "Y"))
# chrom_lengths <- seqlengths(genome)[1:24]
# X <- cbind("baseline" = 1, as.matrix(mcols(gr_tumor)[, -c(1,2,3, 4)]))
# X_totals <- as.matrix(df_areas[, -c(1,2,3)])
# rownames(X_totals) <- df_areas$sample
# X_totals2 <- cbind("baseline" = sum(1 * chrom_lengths), X_totals)
# Mutations <- as.data.frame(mcols(gr_tumor)[, c(2,3)])
# Mutations$sample <- factor(Mutations$sample, levels = rownames(X_totals))
# Mutations$channel <- as.factor(Mutations$channel)
# colnames(X_totals2)
# library(SigPoisProcess)
# a <- 3
# alpha <- 1.01
# K <- 12
# J <- length(unique(Mutations$sample))
# a0 <- a * J - 1
# b0 <- a * 0.001 * J
#
#
# set.seed(10)
# out_Proc_MAP2 <- SigPoisProcess(Mutations = Mutations,
#                                 X = X,
#                                 X_total = X_totals2,
#                                 method = "map",
#                                 prior_params = list(a = a, alpha = alpha,
#                                                     a0 = a0, b0 = b0),
#                                 K = K,
#                                 controls = list(tol = 1e-6,
#                                                 maxiter = 1200,
#                                                 merge_move = TRUE))
#
# round(out_Proc_MAP2$Sig_Cov_prob[1, ,], 5)
# j<- 3
# round(out_Proc_MAP2$Betas[,j,] * out_Proc_MAP2$Xbar[, j, ], 5)
#
# CompressiveNMF::plot_SBS_signature(out_Proc_MAP2$Signatures)
# round(out_Proc_MAP2$Mu, 5)
# match_to_RefSigs(out_Proc_MAP2$Signatures, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
# round(apply(out_Proc_MAP2$Betas * out_Proc_MAP2$Xbar, c(1,3), mean), 5)
#
#
# set.seed(10)
# out_Proc_MLE <- SigPoisProcess(Mutations = Mutations,
#                                 X = X,
#                                 X_total = X_totals2,
#                                 method = "mle",
#                                 prior_params = list(a = a, alpha = alpha),
#                                 K = K,
#                                 controls = list(tol = 1e-6,
#                                                 maxiter = 1000,
#                                                 merge_move = TRUE))
#
# p2 <- CompressiveNMF::plot_SBS_signature(out_Proc_MLE$Signatures)
# round(apply(out_Proc_MLE$Betas * out_Proc_MLE$Xbar, c(1,3), mean), 5)
# round(out_Proc_MLE$Mu, 5)
#
# ggpubr::ggarrange(p1, p2, ncol = 2)
#
# ###########################################
# # EDA
#
# df_tumor <- as.data.frame(gr_tumor)
# df_tumor %>%
#   #filter(sample != "DO218693") %>%
#   group_by(sample, channel) %>%
#   mutate(
#     mut_type = str_extract(channel, "(?<=\\[)[^\\]]+"),
#     across(
#     .cols = CTCF:GC,
#     .fns = mean,
#     .names = "{.col}"
#   )) %>%
#   distinct() %>%
#   dplyr::select(sample, channel, mut_type, DNAse.seq) %>%
#   gather(key = "mark", value = "n", -sample, -channel, -mut_type) %>%
#   #mutate(n = case_when(sqrt(n) > 6.5 ~ 6.5^2,
#   #                     TRUE ~ n)) %>%
#   ggplot() +
#   geom_tile(aes(x = channel, y = sample, fill = n)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
#         axis.text.y = element_text(size = 3),
#         panel.spacing.x = unit(0, "lines"),
#         panel.spacing.y = unit(0.1, "lines"))+
#   facet_grid(mark~mut_type, scales = "free") +
#   scale_fill_viridis_c(name = "signal", option="plasma", direction = -1)
#
#
#
# corrplot::corrplot(cor(as.matrix(mcols(gr_tumor)[, -c(1,2,3)])))
# X <- CompressiveNMF::Mutations_21breast
# pp <- CompressiveNMF::CompressiveNMF_map(X = CompressiveNMF::Mutations_21breast,
#                                          a = 1.0, alpha = 0.5)
# round(c(pp$mapOutput$Mu), 4)
# CompressiveNMF::plot_SBS_signature(pp$Signatures)
# ###########################################
# # Merge with Gene expression (only for liver unfortunately)
# df_igcg <- read_tsv("data/ICGC_rawdata/pcawg_specimen_histology_August2016_v9_donor.tsv")
# df_igcg2 <- read_tsv("data/ICGC_rawdata/pcawg_donor_clinical_August2016_v9.tsv")
# table(df_igcg$icgc_specimen_id %in% df_igcg2$icgc_donor_id)
#
# expr_files <- list.files("data/ICGC_GeneExpr_bigWig/")
# expr_files <- expr_files[gsub(".bigWig", "", expr_files) %in% patients]
#
# unique(patients) %in% df_igcg$icgc_specimen_id
#
# gr <- readRDS("data/ICGC_snp_mutations/Breast-AdenoCa_snp.rds.gzip")
#
# df_igcg %>%
#   filter(icgc_specimen_id %in% colnames(Mexpr2)) %>%
#   pull(histology_abbreviation) %>%
#   table()
#
# unique(patients) %in% colnames(Mexpr)
#
# pats <- unique(gr$sample)
# pats %in% colnames(Mexpr)
#
#
# Xena <- read_tsv(file = "~/October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS")
# sp <- unique(Xena$Sample)[unique(Xena$Sample) %in%colnames(Mexpr)]
#
# df_igcg %>%
#   filter(icgc_specimen_id %in% sp) %>%
#   pull(histology_abbreviation) %>%
#   table()
#
# df_igcg
#
#
#
#
# df_igcg$
# gr1$score == gr2$score
# head(as.matrix(mcols(gr_tumor)[-c(1,2,3)]))
# corrplot::corrplot(cor(as.matrix(mcols(gr_tumor)[-c(1,2,3)])))
#
# colMeans(as.matrix(mcols(gr_tumor)[-c(1,2,3)]))
#
# sc_male <- gr_male$score
# sc_female <- gr_female$score
#
# plot(sc_male, sc_female)
# gp <- keepSeqlevels(gr_male, "chrY", pruning.mode = "coarse")
# hist(gp$score, breaks = 40)
#
#
#
#
#
# ###########################################
# # Merge Baseline mutation rate
# gr_male <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Baseline_male.bigWig")
# gr_female <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Baseline_female.bigWig")
# gr_female <- keepSeqlevels(gr_male, std_chrs[-24], pruning.mode = "coarse")
# scores <- rep(1e-12, length(gr_tumor))
#
# # Merge them with male
# gr_male$score <- transform_scores(gr_male$score)
# scores_male <- rep(1e-12, length(gr_tumor))
# overlaps_male <- findOverlaps(gr_tumor, gr_male)
# scores_male[queryHits(overlaps_male)] <- gr_male[subjectHits(overlaps_male)]$score
#
# # Merge them with female
# gr_female$score <- transform_scores(gr_female$score)
# scores_female <- rep(1e-12, length(gr_tumor))
# overlaps_female <- findOverlaps(gr_tumor, gr_female)
# scores_female[queryHits(overlaps_female)] <- gr_female[subjectHits(overlaps_female)]$score
#
# # Find which is male and which is female
# scores[is_female] <- scores_female[is_female]
# scores[is_male] <- scores_male[is_male]
#
# # Append to new covariate
# gr_tumor$score <- scores
# colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- "baseline"
#
# # Calculate areas
# df_areas <- df_areas %>%
#   left_join(data.frame("donor_sex" = c("female", "male"),
#                        "area" = c(sum(width(gr_male) * gr_male$score),
#                                   sum(width(gr_female) * gr_female$score))),
#             by = "donor_sex")
# colnames(df_areas)[ncol(df_areas)] <- "baseline"
