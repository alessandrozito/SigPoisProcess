##########################################
# Stomach-AdenoCA
##########################################
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
#library(SigPoisProcess)
library(rtracklayer)
library(foreach)
library(doParallel)

devtools::document()
registerDoParallel(20)
maxiter <- 1000
K <- 20

Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
SigPrior <- 3000 * Sigs + 1
#registerDoParallel(4)
##########################################
# CompressiveNMF
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Stomach-AdenoCA_merged_1kb.rds.gzip")

# Get initial estimate for the number of signatures
MutMatrix <- SigPoisProcess::getTotalMutations(TumorData$gr_tumor)
set.seed(10, kind = "L'Ecuyer-CMRG")
MapSol <- foreach(i = 1:10) %dopar% {
  res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = SigPrior, K = 20, a = 1.1, alpha = 1.1)
  res
}

# Use consensus to determine which signatures we should use in the prior
matches <- lapply(MapSol, function(x) {
  match <- match_to_RefSigs(x$Signatures, ref = Sigs)
  sort(match$best[match$cosine > 0.85])
})
sig_names <- sort(unique(unlist(matches)))

# I will initialize the solution at the CompressiveNMF one.
set.seed(10)
MapSol_final <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = 1e8 * Sigs[, sig_names] + 1, K = 0)

CompressiveNMF::plot_weights(MapSol_final$Theta)

# Try a better starting point
K <- ncol(MapSol_final$Signatures)
J <- ncol(MapSol_final$Theta)
p <- ncol(mcols(TumorData$gr_tumor)[-c(1,2,3)])
R_start <- as.matrix(MapSol_final$Signatures)
Xtotals <- as.matrix(TumorData$df_areas[, -c(1,2)])
Betas_start <- array(NA, dim = c(K, J, p))
for(i in 1:p){
  Betas_start[, , i] <- t(apply(MapSol_final$Theta, 1, function(x) x/Xtotals[, i]))
}
Mu_start <- matrix(rowMeans(MapSol_final$Theta))[, rep(1, p)]

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = K,
                      method = "map",
                      init =  SigPoisProcess.init(R_start = R_start,
                                                  Betas_start = Betas_start,
                                                  Mu_start = Mu_start),
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 1000, merge_move = FALSE))
#saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MAP.rds.gzip", compress = "gzip")
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MAP_init.rds.gzip", compress = "gzip")

out <- readRDS("~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MAP.rds.gzip")
plot(out)


##########################################
# Standard Poisson-Process factorization
##########################################

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = 12,
                      method = "map",
                      compressType = "sig",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))

saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MAP_sig.rds.gzip", compress = "gzip")
plot(out)

##########################################
# Standard Poisson-Process factorization
##########################################

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = 12,
                      method = "mle",
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MLE.rds.gzip", compress = "gzip")
plot(out)

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = K,
                      method = "mle",
                      init =  SigPoisProcess.init(R_start = R_start,
                                                  Betas_start = Betas_start,
                                                  Mu_start = Mu_start),
                      controls = SigPoisProcess.controls(maxiter = 1000, merge_move = FALSE))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MLE_init.rds.gzip", compress = "gzip")


outMAP <- readRDS("~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MAP.rds.gzip")
plot(out)
out$mle$logLik_trace[-1]/out$mle$logLik_trace[-19] - 1
match_to_RefSigs(out$Signatures, Sigs)

outMLE <- readRDS("~/SigPoisProcess/output/ICGC/Stomach-AdenoCA_MLE.rds.gzip")
plot(outMLE)





















