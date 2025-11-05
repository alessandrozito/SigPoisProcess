library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(SigPoisProcess)
library(rtracklayer)
library(foreach)
library(doParallel)


#registerDoParallel(20)
maxiter <- 20000
K <- 20

##########################################
# Stomach-AdenoCA
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Stomach-AdenoCA_merged_1kb.rds.gzip")

# Get initial estimate for the number of signatures
#MutMatrix <- SigPoisProcess::getTotalMutations(TumorData$gr_tumor)
# set.seed(10, kind = "L'Ecuyer-CMRG")
# MapSol_sim <- foreach(i = 1:20, .combine = "c") %dopar% {
#   res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 30, a = 1.1, alpha = 1.1)
#   ncol(res$Signatures)
# }
# K <- max(MapSol_sim) + 5
colMeans(as.matrix(mcols(TumorData$gr_tumor)[, -c(1,2,3)]))
corrplot::corrplot(cor(as.matrix(mcols(TumorData$gr_tumor)[, -c(1,2,3)])))

set.seed(42)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                       df_areas = TumorData$df_areas,
                       K = 10,
                       add_baseline = TRUE,
                       prior_params = list(a = 1.1, alpha = 1.1),
                       controls = SigPoisProcess.controls(maxiter = 1000, merge_move = FALSE))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Stomach-AdenoCA.rds.gzip", compress = "gzip")
#out <- readRDS("output/ICGC/Stomach-AdenoCA.rds.gzip")
plot(out)

# How well do we recover the muation loadings?
Covariates <- readRDS("data/ICGC_snp_merged/Stomach-AdenoCA_Mutations_Covariates_Mb.rds.gzip")

Covariates <- readRDS("data/ICGC_snp_merged/Skin-Melanoma_Mutations_Covariates_Mb.rds.gzip")
MutationsMb <- apply(Covariates$MutationsMb, c(2,3), function(x) x)[,colnames(out$Betas), ]
MutationsMb_aggr <- apply(Covariates$MutationsMb, c(2,3), cumsum)[,colnames(out$Betas), ]
Xbar_region <- Covariates$CovariatesMb[,colnames(out$Betas), ]
Xbar_aggr <- apply(Covariates$CovariatesMb, c(2,3), cumsum)[,colnames(out$Betas), ]

Xest <- estimateTotalCounts.SigPoisProcess(out)
Xtot <- getTotalMutations(TumorData$gr_tumor)
plot(Xtot, Xest[, colnames(Xtot)])

)
# Multiply by each patients
R <- out$Signatures
Betas <- out$Betas
TotalCounts <- RegionCounts <- array(NA, dim = dim(MutationsMb_aggr), dimnames = dimnames(MutationsMb_aggr))
for(j in 1:ncol(TotalCounts)){
  TotalCounts[, j, ] <- t(R %*% Betas[, j, ] %*% t(Xbar_aggr[, j, ]))
  RegionCounts[, j, ] <- t(R %*% Betas[, j, ] %*% t(Xbar_region[, j, ]))
}
i <- 3000
plot(t(TotalCounts[i, ,]), t(MutationsMb_aggr[i, ,]))
plot(t(RegionCounts[i, ,]), t(MutationsMb[i, ,]))

j <- 37
plot(rowSums(MutationsMb_aggr[, j,]), type = "l")
lines(rowSums(TotalCounts[, j,]), type = "l", col = "red")



set.seed(10)
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:24]
gr_baseline <- TumorData$gr_tumor
mcols(gr_baseline) <- mcols(gr_baseline)[, c(1,2,3)]
gr_baseline$baseline = 1
df_areas_baseline <- TumorData$df_areas[, c(1,2)]
df_areas_baseline$baseline <- sum(chrom_lengths)
set.seed(10)
out_base <- SigPoisProcess(gr_Mutations = gr_baseline,
                      df_areas = df_areas_baseline,
                      K = 10,
                      add_baseline = FALSE,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 5000, merge_move = FALSE))

plot(out_base)
R <- out_base$Signatures
Betas <- out_base$Betas
TotalCounts_base <- RegionCounts_base <- array(NA, dim = dim(MutationsMb_aggr), dimnames = dimnames(MutationsMb_aggr))
for(j in 1:ncol(TotalCounts_base)){
  TotalCounts_base[, j, ] <- t(R %*% Betas[, j, ] %*% t(Xbar_aggr[, j, 1]))
  RegionCounts_base[, j, ] <- t(R %*% Betas[, j, ] %*% t(Xbar_region[, j, 1]))
}
i <- 10
plot(t(TotalCounts_base[i, ,]), t(MutationsMb_aggr[i, ,]))
plot(t(RegionCounts_base[i, ,]), t(MutationsMb[i, ,]))

j <- 3
plot(rowSums(MutationsMb_aggr[, j,]), type = "l")
lines(rowSums(TotalCounts[, j,]), type = "l", col = "red")
lines(rowSums(TotalCounts_base[, j,]), type = "l", col = "blue")

j <- 3
plot(rowSums(MutationsMb[, j,]), type = "l")
lines(rowSums(RegionCounts[, j,]), type = "l", col = "red")
lines(rowSums(RegionCounts_base[, j,]), type = "l", col = "blue")



Covariates <- readRDS("data/ICGC_snp_merged/Skin-Melanoma_Mutations_Covariates_Mb.rds.gzip")
MutationsMb <- apply(Covariates$MutationsMb, c(2,3), function(x) x)
MutationsMb_aggr <- apply(Covariates$MutationsMb, c(2,3), cumsum)
Xbar_region <- Covariates$CovariatesMb
Xbar_aggr <- apply(Covariates$CovariatesMb, c(2,3), cumsum)

tmp <- apply(MutationsMb, c(1,3), sum)
tmp2 <- apply(Xbar_region, c(1,3), mean)
nMuts <- rowSums(tmp)
corrplot(cor(cbind("C>A" = rowSums(tmp[, grep("C>A", colnames(tmp))]),
                   "C>T" = rowSums(tmp[, grep("C>T", colnames(tmp))]),
                   "C>G" = rowSums(tmp[, grep("C>G", colnames(tmp))]),
                   "T>A" = rowSums(tmp[, grep("T>A", colnames(tmp))]),
                   "T>G" = rowSums(tmp[, grep("T>G", colnames(tmp))]),
                   "T>C" = rowSums(tmp[, grep("T>C", colnames(tmp))]),
                   tmp2)))

plot(nMuts, tmp2[, "CTCF"])
plot(rowSums(tmp[, grep("T>C", colnames(tmp))]), tmp2[, "H3K9me3"])


Covariates_skin <- readRDS("data/ICGC_snp_merged/Skin-Melanoma_Mutations_Covariates_Mb.rds.gzip")
Covariates_stomach <- readRDS("data/ICGC_snp_merged/Stomach-AdenoCA_Mutations_Covariates_Mb.rds.gzip")
Covariates_liver <- readRDS("data/ICGC_snp_merged/Liver-HCC_Mutations_Covariates_Mb.rds.gzip")
Covariates_breast <- readRDS("data/ICGC_snp_merged/Breast-Cancer_Mutations_Covariates_Mb.rds.gzip")
Covariates_Prost <- readRDS("data/ICGC_snp_merged/Prost-AdenoCA_Mutations_Covariates_Mb.rds.gzip")

sumMutations <- function(MutationsMb, name) {
  tmp <- apply(MutationsMb, c(1,3), sum)
  names_cols <- paste(c("C>A", "C>T", "C>G", "T>A", "T>G", "T>C", "all"), name, sep = "_")
  ress <- cbind(rowSums(tmp[, grep("C>A", colnames(tmp))]),
        rowSums(tmp[, grep("C>T", colnames(tmp))]),
        rowSums(tmp[, grep("C>G", colnames(tmp))]),
        rowSums(tmp[, grep("T>A", colnames(tmp))]),
        rowSums(tmp[, grep("T>G", colnames(tmp))]),
        rowSums(tmp[, grep("T>C", colnames(tmp))]),
        rowSums(tmp))
  colnames(ress) <- names_cols
  return(ress)
}


Muts_all <- cbind(sumMutations(Covariates_skin$MutationsMb, "skin"),
      sumMutations(Covariates_breast$MutationsMb, "breast"),
      sumMutations(Covariates_liver$MutationsMb, "liver"),
      sumMutations(Covariates_Prost$MutationsMb, "Prost"),
      sumMutations(Covariates_stomach$MutationsMb, "stomach"))

corrplot(cor(Muts_all[, grepl("all", colnames(Muts_all))]))

plot(Muts_all[, "all_skin"]/mean(Muts_all[, "all_skin"]), type = "l")
lines(Muts_all[, "all_breast"]/mean(Muts_all[, "all_breast"]), col = "red")

plot(Muts_all[, grepl("all", colnames(Muts_all))][, c(1,5)])


as.data.frame(Muts_all[, grepl("all", colnames(Muts_all))]) %>%
  gather(key = "key", value = "value") %>%
  group_by(key) %>%
  mutate(value = value/mean(value),
         region = row_number()) %>%
  ggplot()+
  geom_line(aes(x = region, y = value)) +
  facet_grid(key~.)

sqrt(sum((MutationsMb_aggr - TotalCounts)^2))
sqrt(sum((MutationsMb_aggr - TotalCounts_base)^2))

j <- 37
pp1 <- reshape2::melt((MutationsMb_aggr - TotalCounts)[, j, ])
pp1$type = "FullModel"
pp2 <- reshape2::melt((MutationsMb_aggr - TotalCounts_base)[, j, ])
pp2$type = "Baseline"
pp <- rbind(pp1, pp2)
pp$Var3 <- str_extract(pp$Var2, "(?<=\\[)[^\\]]+")
ggplot(pp) +
  geom_line(aes(x = Var1, y = value, color = Var3, linetype = Var2)) +
  scale_linetype_manual(values = rep("solid", 96)) +
  guides(linetype = "none") +
  scale_color_manual(values = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) +
  facet_grid(Var3 ~ type, scale = "free") +
  geom_hline(yintercept = 0, col = "blue", linetype = "dashed")+
  theme_bw()
pp %>%
  group_by(type) %>%
  summarise(rmse = sqrt(sum(value^2)))

match_to_RefSigs(sigs = out$Signatures, out_base$Signatures)


##########################################
# Breast-AdenoCa
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Breast-Cancer_merged_1kb.rds.gzip")

# Get initial estimate for the number of signatures
#MutMatrix <- SigPoisProcess::getTotalMutations(TumorData$gr_tumor)
# set.seed(10, kind = "L'Ecuyer-CMRG")
# MapSol_sim <- foreach(i = 1:20, .combine = "c") %dopar% {
#   res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 30, a = 1.1, alpha = 1.1)
#   ncol(res$Signatures)
# }
# K <- max(MapSol_sim) + 5

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = K,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = maxiter))
saveRDS(object = out, file = "output/ICGC/Breast-Cancer.rds.gzip", compress = "gzip")
out <- readRDS("output/ICGC/Breast-Cancer.rds.gzip")
plot(out)

##########################################
# Skin-Melanoma
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Skin-Melanoma_merged_1kb.rds.gzip")

# Get initial estimate for the number of signatures
#MutMatrix <- SigPoisProcess::getTotalMutations(TumorData$gr_tumor)
# set.seed(10, kind = "L'Ecuyer-CMRG")
# MapSol_sim <- foreach(i = 1:20, .combine = "c") %dopar% {
#   res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 30, a = 1.1, alpha = 1.1)
#   ncol(res$Signatures)
# }
# K <- max(MapSol_sim) + 5

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = 20,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 200, merge_move = FALSE))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Skin-Melanoma.rds.gzip", compress = "gzip")
out <- readRDS("output/ICGC/Skin-Melanoma.rds.gzip")
plot(out)

# How well do we recover the muation loadings?
Covariates <- readRDS("data/ICGC_snp_merged/Skin-Melanoma_Mutations_Covariates_Mb.rds.gzip")

MutationsMb_aggr <- apply(Covariates$MutationsMb, c(2,3), cumsum)[,colnames(out$Betas), ]
Xbar_aggr <- apply(Covariates$CovariatesMb, c(2,3), cumsum)[,colnames(out$Betas), ]

Xest <- estimateTotalCounts.SigPoisProcess(out)
Xtot <- getTotalMutations(TumorData$gr_tumor)
plot(Xtot, Xest[, colnames(Xtot)])

# Multiply by each patients
R <- out$Signatures
Betas <- out$Betas
TotalCounts <- array(NA, dim = dim(MutationsMb_aggr), dimnames = dimnames(MutationsMb_aggr))
for(j in 1:ncol(TotalCounts)){
  TotalCounts[, j, ] <- t(R %*% Betas[, j, ] %*% t(Xbar_aggr[, j, ]))
}
i <- 160
plot(t(TotalCounts[i, ,]), t(MutationsMb_aggr[i, ,]))


set.seed(10)
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:24]
gr_baseline <- TumorData$gr_tumor
mcols(gr_baseline) <- mcols(gr_baseline)[, c(1,2,3)]
gr_baseline$baseline = 1
df_areas_baseline <- TumorData$df_areas[, c(1,2)]
df_areas_baseline$baseline <- sum(chrom_lengths)

set.seed(10)
out_base <- SigPoisProcess(gr_Mutations = gr_baseline,
                           df_areas = df_areas_baseline,
                           K = 20,
                           add_baseline = FALSE,
                           prior_params = list(a = 1.1, alpha = 1.1),
                           controls = SigPoisProcess.controls(maxiter = 2000, merge_move = FALSE))

R <- out_base$Signatures
Betas <- out_base$Betas
TotalCounts_base <- array(NA, dim = dim(MutationsMb_aggr), dimnames = dimnames(MutationsMb_aggr))
for(j in 1:ncol(TotalCounts_base)){
  TotalCounts_base[, j, ] <- t(R %*% Betas[, j, ] %*% t(Xbar_aggr[, j, 1]))
}
i <- 10
plot(t(TotalCounts_base[i, ,]), t(MutationsMb_aggr[i, ,]))



##########################################
# Prost-AdenoCA
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Prost-AdenoCA_merged_1kb.rds.gzip")

# Get initial estimate for the number of signatures
#MutMatrix <- SigPoisProcess::getTotalMutations(TumorData$gr_tumor)
# set.seed(10, kind = "L'Ecuyer-CMRG")
# MapSol_sim <- foreach(i = 1:20, .combine = "c") %dopar% {
#   res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 30, a = 1.1, alpha = 1.1)
#   ncol(res$Signatures)
# }
# K <- max(MapSol_sim) + 5

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = K,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = maxiter))
saveRDS(object = out, file = "output/ICGC/Prost-AdenoCA.rds.gzip", compress = "gzip")
out <- readRDS("output/ICGC/Prost-AdenoCA.rds.gzip")
plot(out)
sigs <- out$Signatures[, rowMeans(out$Mu) > 0.001]
round(rowMeans(out$Mu), 5)
match_to_RefSigs(sigs, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
plot_96SBSsig(sigs)

round(out$Betas *out$Xbar, 5)[, ,9]
##########################################
# Liver-HCC
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Liver-HCC_merged_1kb.rds.gzip")

# Get initial estimate for the number of signatures
#MutMatrix <- SigPoisProcess::getTotalMutations(TumorData$gr_tumor)
# set.seed(10, kind = "L'Ecuyer-CMRG")
# MapSol_sim <- foreach(i = 1:20, .combine = "c") %dopar% {
#   res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 30, a = 1.1, alpha = 1.1)
#   ncol(res$Signatures)
# }
# K <- max(MapSol_sim) + 5

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = K,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))
saveRDS(object = out, file = "output/ICGC/Liver-HCC.rds.gzip", compress = "gzip")

out <- readRDS("output/ICGC/Liver-HCC.rds.gzip")
plot(out)
sigs <- out$Signatures[, rowMeans(out$Mu) > 0.001]
match_to_RefSigs(sigs, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
plot_96SBSsig(sigs)

j <- 50
plot(rowSums(MutationsMb_aggr[, j,]), type = "l")
lines(rowSums(TotalCounts[, j,]), type = "l", col = "red")
lines(rowSums(TotalCounts_base[, j,]), type = "l", col = "blue")

sqrt(sum((MutationsMb_aggr - TotalCounts)^2))
sqrt(sum((MutationsMb_aggr - TotalCounts_base)^2))

j <- 3
pp1 <- reshape2::melt((MutationsMb_aggr - TotalCounts)[, j, ])
pp1$type = "FullModel"
pp2 <- reshape2::melt((MutationsMb_aggr - TotalCounts_base)[, j, ])
pp2$type = "Baseline"
pp <- rbind(pp1, pp2)

pp$Var3 <- str_extract(pp$Var2, "(?<=\\[)[^\\]]+")
ggplot(pp) +
  geom_line(aes(x = Var1, y = value, color = Var3, linetype = Var2)) +
  scale_linetype_manual(values = rep("solid", 96)) +
  guides(linetype = "none") +
  scale_color_manual(values = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) +
  facet_grid(Var3 ~ type, scale = "free") +
  geom_hline(yintercept = 0, col = "blue", linetype = "dashed")+
  theme_bw()
pp %>%
  group_by(type) %>%
  summarise(rmse = sqrt(sum(value^2)))




round(MutationsMb_aggr[3113, 2, ] - TotalCounts[3113, 2, ], 50)

nrow(MutationsMb_aggr)






