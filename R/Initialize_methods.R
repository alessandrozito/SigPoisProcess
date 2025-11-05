

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
#registerDoParallel(20)
#maxiter <- 10000
#K <- 20

Sigs <-CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
SigPrior <- 3000 * CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37 + 1
#registerDoParallel(4)
##########################################
# CompressiveNMF
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Liver-HCC_merged_1kb.rds.gzip")
MutMatrix <- getTotalMutations(TumorData$gr_tumor, df_areas = TumorData$df_areas)

sol <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 30, a = 1.1, alpha = 1.1)

plot_96SBSsig(sol$Signatures)

# Now, lets calculate the average value of the covariates for each patient and mutation
R <- sol$Signatures
Theta <- sol$Theta

df_tumor <- as.data.frame(TumorData$gr_tumor)

channels <- rownames(R)
samples <- colnames(Theta)
# Pick a single signature
sig_attr <- rep(NA, nrow(df_tumor))
for(s in 1:nrow(df_tumor)){
  if(s%%100000 == 0){
    print(s)
  }
  i <- which(channels == df_tumor$channel[s])
  j <- which(samples == df_tumor$sample[s])
  sig_attr[s] <- which.max((R[i, ] * Theta[, j])/sum(R[i, ] * Theta[, j]))
}

df_tumor$sig_attr <- sig_attr

df_aggr <- df_tumor %>%
  group_by(channel, sample, sig_attr) %>%
  summarise(
    n = n(),
    across(.cols = baseline:GCsq, .fns = mean, .names = "{.col}"))


k <- 5
j <- which(samples == "DO23508")

tmp <- df_aggr %>%
  filter(sample == "DO23508", sig_attr == 5)

yX <- matrix(0, nrow = 96, ncol = 16)
colnames(yX) <- colnames(tmp)[-c(1:3)]
rownames(yX) <- rownames(R)
yX[tmp$channel, ] <- as.matrix(tmp[, -c(1:3)])
yX[, -1] <- yX[, -1] + matrix(rgamma(length(yX[, -1]), shape = 1, rate = 1000), nrow = nrow(yX))

glm(n ~. -1, family = poisson(link = "identity"),
    data = as.data.frame(yX),
    start = rep(1, ncol(yX) - 1))

library(nnls)
x <- nnls(A = yX[, -1] * R[, rep(k, ncol(yX[, -1]))], b = yX[, 1])$x
names(x) <- colnames(yX[, -1])
x

TumorData$gr_tumor$baseline

Theta[k, j]





































