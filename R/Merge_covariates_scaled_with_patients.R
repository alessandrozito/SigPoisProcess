# This script merges covariates at the 1kb resolution with the patients.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)


merge_with_patients <- function(gr_tumor, tumor, files_tumor){

  # Load the genome 1MegaBase
  genome <- BSgenome.Hsapiens.UCSC.hg19
  chrom_lengths <- seqlengths(genome)[1:24]

  # Remove mutations from blacklist
  blacklist <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/hg19-blacklist.v2.bed")
  overlaps <- findOverlaps(gr_tumor, blacklist)
  gr_tumor <- gr_tumor[-queryHits(overlaps)]

  # Find male and female in patients
  patients <- gr_tumor$sample
  is_male <- which(patients %in% patiens_male)
  is_female <- which(patients %in% patiens_female)

  # Names of the marks
  name_marks <- gsub("_1kb_standardized.bigWig.rds.gzip", "", files_tumor)
  name_marks <- gsub(paste0(tumor, "_"), "", name_marks)
  mark <- unique(sub("^(female|male)_", "", name_marks))

  # Handle sex of the patient/histone
  donor_sex <- c("male", "female")
  sex <- ifelse(grepl("female_", name_marks), "female",
                ifelse(grepl("male_", name_marks), "male", NA))

  # Areas data.frame
  df_areas <- data.frame(sample = unique(patients)) %>%
    left_join(df_clinical %>% dplyr::select(icgc_donor_id, donor_sex),
              by = c("sample" = "icgc_donor_id"))

  # Areas file data.frame
  df_files <- data.frame(sample = unique(patients)) %>%
    left_join(df_clinical %>% dplyr::select(icgc_donor_id, donor_sex),
              by = c("sample" = "icgc_donor_id"))

  ###########################################
  # Merge histones and CTCF
  covariates_to_save <- c("score_scaled", "score_plus", "score_minus")
  for(i in 1:length(mark)){
    print(mark[i])
    file_male <-  files_tumor[grepl(mark[i], files_tumor) & grepl("male", files_tumor) & !grepl("female", files_tumor)]
    file_female <-  files_tumor[grepl(mark[i], files_tumor) & grepl("female", files_tumor)]

    if(length(file_female) == 0 & length(file_male) == 1) {
        file_female <- file_male
    } else if (length(file_female) == 1 & length(file_male) == 0){
        file_male <- file_female
    }

    scores <- matrix(0, nrow = length(gr_tumor), ncol = 3)
    colnames(scores) <- covariates_to_save

    # Load both files
    gr_male <-  readRDS(paste0(dir_covariates, file_male))
    gr_female <- readRDS(paste0(dir_covariates, file_female))
    gr_female[gr_female@seqnames=="chrY"]$bin_weight <- 0

    # Merge them with male
    scores_male <- matrix(0, nrow = length(gr_tumor), ncol = 3)
    overlaps_male <- findOverlaps(gr_tumor, gr_male)
    scores_male[queryHits(overlaps_male), ] <- as.matrix(mcols(gr_male[subjectHits(overlaps_male)])[, covariates_to_save])

    # Merge them with female
    overlaps_female <- findOverlaps(gr_tumor, gr_female)
    scores_female <- matrix(0, nrow = length(gr_tumor), ncol = 3)
    scores_female[queryHits(overlaps_female), ] <- as.matrix(mcols(gr_female[subjectHits(overlaps_female)])[, covariates_to_save])

    # Find which is male and which is female
    scores[is_female, ] <- scores_female[is_female, ]
    scores[is_male, ] <- scores_male[is_male, ]

    # Append to gr
    colnames(scores) <- paste0(mark[i], "_", colnames(scores))
    mcols(gr_tumor) <- cbind(mcols(gr_tumor), scores)

    # Calculate areas
    df_tmp <- data.frame("donor_sex" = c("male", "female"),
                         "score_plus" = c(sum(gr_male$bin_weight * gr_male$score_plus),
                                          sum(gr_female$bin_weight * gr_female$score_plus)),
                         "score_minus" = c(sum(gr_male$bin_weight * gr_male$score_minus),
                                           sum(gr_female$bin_weight * gr_female$score_minus)))
    colnames(df_tmp)[-1] <-  paste0(mark[i], "_", colnames(df_tmp)[-1])
    df_areas <- df_areas %>%
      left_join(df_tmp, by = "donor_sex")


    # Append the files
    df_files <- df_files %>%
      left_join(data.frame("donor_sex" = c("female", "male"),
                           "file" = c(paste0(dir_covariates, file_female),
                                      paste0(dir_covariates, file_male))))
    colnames(df_files)[ncol(df_files)] <- mark[i]

  }


  ########################################### SAVE FILE!
  TumorData <- list(gr_tumor = gr_tumor,
                    df_areas = df_areas,
                    df_files = df_files)

  saveRDS(TumorData, paste0("~/SigPoisProcess/data/ICGC_snp_merged/", tumor, "_merged_1kb_standardized.rds.gzip"), compress = "gzip")
}


#--------- Load the cancer
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Load clinical data
df_clinical <- read_tsv("~/SigPoisProcess/data/ICGC_rawdata/pcawg_donor_clinical_August2016_v9.tsv")
patiens_male <- df_clinical$icgc_donor_id[df_clinical$donor_sex == "male"]
patiens_female <- df_clinical$icgc_donor_id[df_clinical$donor_sex == "female"]

# Load covariate names
dir_covariates <- "~/SigPoisProcess/data/Covariates_standardized/"
gr_files <- list.files(dir_covariates)


##########################################
# Stomach-AdenoCA
##########################################
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")
tumor <- "Stomach-AdenoCA"

files_shared <- c("male_GC_1kb_standardized.bigWig", "female_GC_1kb_standardized.bigWig",
                  "male_RepliTiming_1kb_standardized.bigWig", "female_RepliTiming_1kb_standardized.bigWig",
                  "male_NuclOccup_1kb_standardized.bigWig", "female_NuclOccup_1kb_standardized.bigWig")
files_shared <- paste0(files_shared, ".rds.gzip")
files_tumor <- gr_files[grepl(tumor, gr_files) | gr_files %in% files_shared]

merge_with_patients(gr_tumor, tumor, files_tumor)


TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Stomach-AdenoCA_merged_1kb_standardized.rds.gzip")


gr_tumor <- TumorData$gr_tumor
mcols(gr_tumor) <- mcols(gr_tumor)[, grepl("sample|channel|minus|plus", colnames(mcols(gr_tumor)))]
mcols(gr_tumor) <- mcols(gr_tumor)[, !grepl("DNAse", colnames(mcols(gr_tumor)))]
colnames(mcols(gr_tumor)) <- gsub("score_minus", "minus", colnames(mcols(gr_tumor)))
colnames(mcols(gr_tumor)) <- gsub("score_plus", "plus", colnames(mcols(gr_tumor)))

df_areas <- TumorData$df_areas
df_areas <- df_areas[, grepl("sample|minus|plus", colnames(df_areas))]
colnames(df_areas) <- gsub("score_minus", "minus", colnames(df_areas))
colnames(df_areas) <- gsub("score_plus", "plus", colnames(df_areas))
df_areas <- df_areas[, !grepl("DNAse", colnames(df_areas))]
colnames(df_areas)

set.seed(10)
out <- SigPoisProcess(gr_Mutations = gr_tumor,
                      df_areas = df_areas,
                      K = 12,
                      method = "map",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 2000, merge_move = TRUE, tol = 1e-5))
plot(out)

tmp <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_H3K36me3_1kb_standardized.bigWig.rds.gzip")
gr_tumor2 <- gr_tumor; df_areas2 <- df_areas
gr_tumor2$baseline <- 1
df_areas2$baseline <- sum(tmp$bin_weight)


set.seed(10)
out <- SigPoisProcess(gr_Mutations = gr_tumor2,
                      df_areas = df_areas2,
                      K = 10,
                      method = "map",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 2000, merge_move = FALSE, tol = 1e-5))
plot(out)

refSigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
match_to_RefSigs(out$Signatures, refSigs)


df <- as.data.frame(gr_tumor)

df_temp <- df[df$sample == "DO222286", ]
plot(df_temp$RepliTiming_score_minus)
lines(df_temp$H3K9me3_score_plus, col = "red")

sort(table(apply(df_temp[, -c(1:7)], 1, function(x) names(which.max(x)))))

round(out$Betas[, "DO222286",] * out$Xbar[, 1, ], 4)



### Try multiplicative effects now.

TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Stomach-AdenoCA_merged_1kb_standardized.rds.gzip")

gr_Mutations <- TumorData$gr_tumor
df_files <- TumorData$df_files

build_scoreTracks <- function(df_files, tilewidth = 5e5){

  genome_part <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                           tilewidth = tilewidth,
                           cut.last.tile.in.chrom = TRUE)

  samples <- df_files$sample
  covs <- colnames(df_files)[-c(1,2)]
  Files_mat <- as.matrix(df_files)[, -c(1,2)]

  J <- nrow(df_files)
  scoreTracks <- array(NA, c(length(genome_part), J, length(covs)),
                       dimnames = list(1:length(genome_part), samples, covs))
  Bins <- matrix(NA, nrow = length(genome_part), ncol = J,
                 dimnames = list(1:length(genome_part), samples))

  files_all <- unique(Files_mat)
  for(file in files_all){
    print(file)
    gr_cov <- readRDS(file)
    overlaps <- findOverlaps(gr_cov, genome_part)
    df_tmp <- as.data.frame(gr_cov) %>%
      mutate(region = subjectHits(overlaps)) %>%
      group_by(region) %>%
      summarise(score = sum(score_scaled * bin_weight)/sum(bin_weight),
                bin = sum(bin_weight),
                score = case_when(bin == 0 ~ 0,
                                  TRUE ~ score))
    which_has_file <- which(Files_mat == file, arr.ind = TRUE)
    for(i in 1:nrow(which_has_file)){
      j <- which_has_file[i, 1]
      l <- which_has_file[i, 2]
      Bins[, j] <- df_tmp$bin
      scoreTracks[, j, l] <- df_tmp$score
    }
  }
  return(list(scoreTracks = scoreTracks, Bins = Bins))
}

tracks_all <- build_scoreTracks(df_files)
dim(tracks_all$scoreTracks)

gr_Mutations <- TumorData$gr_tumor
mcols(gr_Mutations) <- mcols(gr_Mutations)[grepl("_score_scaled|sample|channel", colnames(mcols(gr_Mutations)))]
colnames(mcols(gr_Mutations)) <- gsub("_score_scaled","", colnames(mcols(gr_Mutations)))
MutMatrix <- getTotalMutations(gr_Mutations)

scoreTraks <- tracks_all$scoreTracks[, colnames(MutMatrix)  ,]
Bins <- tracks_all$Bins[, colnames(MutMatrix)]

refSigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
resMap <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = 3000 * refSigs + 1, K = 10)


R_start <- resMap$Signatures
Theta_start <- resMap$Theta[, colnames(MutMatrix)]/mean(colSums(Bins))
K <- ncol(R_start)
colnames(R_start) <- rownames(Theta_start) <- paste0("SigN", sprintf("%02d", 1:K))
init <- SigPoisProcess_mult.init(R_start = R_start, Theta_start = Theta_start)

scoreTracks <- apply(tracks_500Kb$scoreTracks, c(2,3), scale)
scoreTracks <- tracks_500Kb$scoreTracks
Bins <- tracks_500Kb$Bins
set.seed(10)
out_stomach <- SigPoisProcess_mult(gr_Mutations = gr_Mutations,
                                   scoreTracks = scoreTracks,
                                   Bins = Bins,
                                   K = K,
                                   method = "mle",
                                   init = init)
round(out_stomach$Betas, 2)

if(TRUE){
  X <- as.data.frame(GenomicRanges::mcols(gr_Mutations)) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()
  MutMatrix <- getTotalMutations(gr_Mutations)
  # Extract mutation data
  Mutations <- as.data.frame(GenomicRanges::mcols(gr_Mutations)) %>%
    dplyr::select(sample, channel)
  Mutations$sample <- factor(Mutations$sample, levels = unique(Mutations$sample))
  Mutations$channel <- as.factor(Mutations$channel)

  # Extract channel ID and mutation ID
  channel_id <- as.numeric(Mutations$channel) -1
  sample_id <- as.numeric(Mutations$sample) - 1
  Betas_start <- matrix(0, nrow = K, ncol = p)
}

Betas <- Update_Theta_Betas_loglinear(R = R_start,
                                     Theta = Theta_start,
                                     Betas_start = Betas_start,
                                     Xtotal = scoreTracks,
                                     Bins = t(Bins),
                                     X = X,
                                     channel_id = channel_id,
                                     sample_id = sample_id,
                                     iter = 20)
Betas
plot_96SBSsig(R_start)
plot_96SBSsig(out_stomach$Signatures)
round(out_stomach$Thetas * mean(colSums(Bins)), 5)

plot(R_start %*% (out_stomach$Thetas[, colnames(MutMatrix)] * mean(colSums(Bins))), MutMatrix)
colnames(out_stomach$Thetas) == colnames(MutMatrix)


tracks_500Kb <- build_scoreTracks(df_files = df_files, tilewidth = 500000)
tracks_50Kb <- build_scoreTracks(df_files = df_files, tilewidth = 50000)
tracks_10Kb <- build_scoreTracks(df_files = df_files, tilewidth = 10000)


corrplot(cor(tracks_50Kb$scoreTracks[, 1, ]), method = "number")

p <- dim(tracks_500Kb$scoreTracks)[3]

tmp <- tracks_500Kb$scoreTracks[, 1, ]

Theta_start
betas <- rnorm(p, sd = 0.2)
sum(0.001 * 2e-7 * tracks_50Kb$Bins[, 1] * exp(tracks_50Kb$scoreTracks[, 1, ] %*% betas))
sum(0.001 * 2e-7 * tracks_500Kb$Bins[, 1] * exp(tracks_500Kb$scoreTracks[, 1, ] %*% betas))

tmp <- 1 + tracks_500Kb$scoreTracks[, 1, ] %*% betas + c(0.5 * betas %*% crossprod(tracks_500Kb$scoreTracks[, 1, ]) %*% betas)
c(crossprod(tracks_500Kb$Bins[, 1], tmp))


DNase <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_DNAse_seq_1kb_standardized.bigWig.rds.gzip")
NuclOccup <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_NuclOccup_1kb_standardized.bigWig.rds.gzip")

subseq <- sample(2700000, size = 20000)
#subseq <- 1000:5000
cor(NuclOccup[subseq]$score_tr, DNase[subseq]$score_tr)
plot(NuclOccup[subseq]$score_tr, DNase[subseq]$score_tr)
subseq <- 5000:20000
par(mfrow = c(2,1))
plot(NuclOccup[subseq]$score_tr, type = "l"); plot(DNase[subseq]$score_tr, type = "l")
par(mfrow = c(1,1))

"~/SigPoisProcess/data/ENCODE_pvalues/NuclOccupancy/GSM920557_hg19_wgEncodeSydhNsomeK562Sig.bigWig"
chain <- import.chain("hg38ToHg19.over.chain.gz")
gg <- GRanges(seqnames = "chr1",
        ranges = IRanges(start = 1000000, end = 1500000))
NuclOccup_all <- import("~/SigPoisProcess/data/ENCODE_pvalues/NuclOccupancy/GSM920557_hg19_wgEncodeSydhNsomeK562Sig.bigWig",
                        which= gg)
DNAse_all <- unlist(liftOver(import("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/ENCFF865IXT.bigWig",
                which = gg), chain = chain))

overlaps <- findOverlaps(NuclOccup_all, DNAse_all)
no <- NuclOccup_all[queryHits(overlaps)]$score
dd <- DNAse_all[subjectHits(overlaps)]$score
base::plot(dd, no,  col = rgb(0, 0, 0, alpha = 0.1), pch = 20)
base::plot(log2(dd+0.1), log2(no + 0.1))
cor(dd, no)
ids <- sample(1:length(dd), size = 1000)
plot(dd[ids], no[ids], xlim = c(0, 0.3))
abline(lm(no ~ dd))
cor(log2(dd+0.1), log2(no + 0.1))
abline(lm(log2(no + 0.1) ~ log2(dd+0.1)))
plot(scale(dd), type = "l"); lines(scale(no), type = "l", col = rgb(1, 0, 0, alpha = 0.2))

#-------------------------------------------------------------------------------
gg <- GRanges(seqnames = "chr1",
              ranges = IRanges(start = 1000000, end = 1500000))
gr1_all <- import("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/ENCFF890TDY.bigWig",
                        which= gg)
gr2_all <- import("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/ENCFF862XUL.bigWig",
                                    which = gg)

overlaps <- findOverlaps(gr1_all, gr2_all)
no <- NuclOccup_all[queryHits(overlaps)]$score
dd <- H3K9me3_all[subjectHits(overlaps)]$score
ids <- sample(1:length(dd), size = 2000)
plot(dd[ids], no[ids])
abline(lm(asinh(no) ~ asinh(dd)))
means <- rep(gr2_all$score, times = width(gr2_all))
plot(means, type = "l")

cor(log2(dd+0.1), log2(no + 0.1))
abline(lm(log2(no + 0.1) ~ log2(dd+0.1)))
plot(scale(dd), type = "l"); lines(scale(no), type = "l", col = rgb(1, 0, 0, alpha = 0.2))

# Let's do it a 1kb
NuclOccup <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_NuclOccup_1kb_standardized.bigWig.rds.gzip")
DNAse <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_DNAse_seq_1kb_standardized.bigWig.rds.gzip")

ids <- sample(1:length(DNAse), size = 3000)
plot(DNAse[ids]$score_scaled, NuclOccup[ids]$score_scaled)
abline(lm(NuclOccup[ids]$score_scaled~DNAse[ids]$score_scaled))
cor(DNAse[ids]$score_scaled, NuclOccup[ids]$score_scaled)


gr1 <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_RepliTiming_1kb_standardized.bigWig.rds.gzip")
gr2 <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_GC_1kb_standardized.bigWig.rds.gzip")

ids <- sample(1:length(gr1), size = 5000)
ggplot(data = data.frame(x = gr1[ids]$score, y = gr2[ids]$score) %>% filter(x>0 & y > 0),
       aes(x=x, y = y)) +
  geom_point(alpha = 0.1) +
  geom_density_2d(alpha = 0.5) +
  theme_bw() +
  geom_smooth()

abline(lm(NuclOccup[ids]$score_scaled~DNAse[ids]$score_scaled))
cor(H3K9me3[ids]$score_scaled, CTCF[ids]$score_scaled)



# Let's calculate their univariate distribution and the distribution at the
# distribution at the mutations

# Replication timing
gr_tumor <- TumorData$gr_tumor
gr_cov <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_RepliTiming_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(gr_cov, gr_tumor)
scores <- gr_cov[queryHits(overlaps)]$score
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score
plot(density(sc_all), ylim = c(0, 0.034))
lines(density(scores), col = "red")


# GC content
gr_cov <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_GC_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(gr_cov, gr_tumor)
scores <- gr_cov[queryHits(overlaps)]$score
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score
plot(density(sc_all), ylim = c(0, 6.5))
lines(density(scores), col = "red")

# Methylation
gr_cov <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_Methylation_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(gr_cov, gr_tumor)
scores <- gr_cov[queryHits(overlaps)]$score
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score
plot(density(sc_all), ylim = c(0, 70))
lines(density(scores), col = "red")

# CTCF
gr_cov <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_CTCF_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(gr_cov, gr_tumor)
scores <- gr_cov[queryHits(overlaps)]$score
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score
plot(density(sc_all))
lines(density(scores), col = "red")
scores <- gr_cov[queryHits(overlaps)]$score_tr
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score_tr
plot(density(sc_all))
lines(density(scores), col = "red")


# Nucelosome occupancy
gr_cov <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_NuclOccup_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(gr_cov, gr_tumor)
scores <- gr_cov[queryHits(overlaps)]$score
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score
plot(density(sc_all))
lines(density(scores), col = "red")
scores <- gr_cov[queryHits(overlaps)]$score_tr
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score_tr
plot(density(sc_all))
lines(density(scores), col = "red")


# H3K9me3
gr_cov <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_H3K9me3_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(gr_cov, gr_tumor)
scores <- gr_cov[queryHits(overlaps)]$score
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score
plot(density(sc_all))
lines(density(scores), col = "red")
scores <- gr_cov[queryHits(overlaps)]$score_tr
sc_all <- gr_cov[gr_cov$bin_weight > 0]$score_tr
plot(density(sc_all))
lines(density(scores), col = "red")










