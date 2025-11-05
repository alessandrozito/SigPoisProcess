##########################################
# Liver-HCC
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

re_initialize <- FALSE
if(re_initialize) {
  # Get initial estimate for the number of signatures
  set.seed(10, kind = "L'Ecuyer-CMRG")
  MapSol <- foreach(i = 1:20) %dopar% {
    res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = SigPrior, K = 20, a = 1.1, alpha = 1.1)
    res
  }

  plot(MutMatrix, MapSol[[13]]$Signatures %*% MapSol[[13]]$Theta)

  # Use consensus to determine which signatures we should use in the prior
  matches <- lapply(MapSol, function(x) {
    match <- match_to_RefSigs(x$Signatures, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37 )
    sort(match$best[match$cosine > 0.9])
  })
  sig_names <- sort(unique(unlist(matches)))

  ncols <- unlist(lapply(MapSol, function(x) ncol(x$Signatures)))
  logPost <- unlist(lapply(MapSol, function(x) x$logPosterior))
  sig_names <- sort(unique(unlist(matches)))
  which.min(logPost)

  sig_names <- c(names((sort(table(unlist(matches)), decreasing = TRUE)/20)[1:14]), "SBS17a", "SBS17b")
  sig_names <- sig_names[sig_names!="SBS44"]
  sig_names[sig_names=="SBS8"] <- "SBS9"
  sig_names[sig_names=="SBS31"] <- "SBS29"

  sig_names <- sig_names[order(as.numeric(sub("SBS([0-9]+).*", "\\1", sig_names)), sub("SBS[0-9]+", "", sig_names))]
  # sig_names <- c("SBS1", "SBS4", "SBS5", "SBS6", "SBS9", "SBS12", "SBS16", "SBS17b", "SBS19", "SBS22a", "SBS29", "SBS35", "SBS40a")
  ## Starting point based on NNLS (compressive NMF)
  set.seed(10)
  MapSol_final <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = 1e8 * Sigs[, sig_names] + 1, K = 0)
  CompressiveNMF::plot_weights(MapSol_final$Theta)
  plot_96SBSsig(MapSol_final$Signatures)

  plot(MutMatrix, MapSol_final$Signatures %*% MapSol_final$Theta)

  saveRDS(MapSol_final, "~/SigPoisProcess/output/ICGC/StartPoint_LiverHcc_v2.rds.gzip")

}

#MapSol_final <- readRDS("~/SigPoisProcess/output/ICGC/StartPoint_LiverHcc.rds.gzip")
MapSol_final <- readRDS("~/SigPoisProcess/output/ICGC/StartPoint_LiverHcc_v2.rds.gzip")

# Starting point
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

#########################################################
# Compressive Poisson point process factorization
#########################################################
# Distribute now the starting point
maxiter <- 5000

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = K,
                      method = "map",
                      init =  SigPoisProcess.init(R_start = R_start,
                                                  Betas_start = Betas_start,
                                                  Mu_start = Mu_start),
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE, tol = 1e-6))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Liver-HCC_MAP_init_v3.rds.gzip", compress = "gzip")

#########################################################
# Standard Poisson-Process factorization
#########################################################

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = K,
                      method = "mle",
                      init =  SigPoisProcess.init(R_start = R_start,
                                                  Betas_start = Betas_start,
                                                  Mu_start = Mu_start),
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE, tol = 1e-6))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Liver-HCC_MLE_init_v3.rds.gzip", compress = "gzip")


outMAP <- readRDS("~/SigPoisProcess/output/ICGC/Liver-HCC_MAP_init_v3.rds.gzip")
plot(outMAP)
ggsave("~/SigPoisProcess/figures/LiverHcc_out_map.pdf", width = 7.47, height= 6.73)
outMLE <- readRDS("~/SigPoisProcess/output/ICGC/Liver-HCC_MLE_init_v3.rds.gzip")
plot(outMLE)
ggsave("~/SigPoisProcess/figures/LiverHcc_out_mle.pdf", width = 7.47, height= 6.73)

# plot(outMAP$map$logPost_trace)
# plot(outMLE$mle$logLik_trace)
# #
# match_to_RefSigs(outMAP$Signatures, ref = Sigs)
#
# plot(estimateTotalCounts.SigPoisProcess(outMLE), MutMatrix)
# plot(estimateTotalCounts.SigPoisProcess(outMAP), MutMatrix)
#
# estimateTotalCounts_v2 <- function(R, Betas, X_total){
#   MatOut <-sapply(1:ncol(Betas), function(j) R %*% Betas[, j, ] %*% X_total[j, ], simplify = TRUE)
#   colnames(MatOut) <- colnames(Betas)
#   rownames(MatOut) <- rownames(R)
#   return(MatOut)
# }
# plot(estimateTotalCounts_v2(R_start, Betas_start, as.matrix(TumorData$df_areas[, -c(1,2)])), MutMatrix)


#
# set.seed(10)
# out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
#                       df_areas = TumorData$df_areas,
#                       K = 15,
#                       method = "map",
#                       prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
#                       controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))
# saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Liver-HCC_MAP.rds.gzip", compress = "gzip")
#
#
# ##########################################
# # Standard Poisson-Process factorization
# ##########################################
#
# set.seed(10)
# out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
#                       df_areas = TumorData$df_areas,
#                       K = 15,
#                       method = "map",
#                       compressType = "sig",
#                       prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
#                       controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))
#
# saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Liver-HCC_MAP_sig.rds.gzip", compress = "gzip")
#
#
# ##########################################
# # Standard Poisson-Process factorization
# ##########################################
#
# set.seed(10)
# out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
#                       df_areas = TumorData$df_areas,
#                       K = 15,
#                       method = "mle",
#                       controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))
# saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Liver-HCC_MLE.rds.gzip", compress = "gzip")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
