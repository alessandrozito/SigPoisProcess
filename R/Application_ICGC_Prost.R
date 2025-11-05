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
#registerDoParallel(20)
maxiter <- 10000
K <- 20

SigPrior <- 3000 * CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37 + 1
#registerDoParallel(4)
##########################################
# CompressiveNMF
##########################################
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Prost-AdenoCA_merged_1kb.rds.gzip")

# # Get initial estimate for the number of signatures
# MutMatrix <- SigPoisProcess::getTotalMutations(TumorData$gr_tumor)
# set.seed(10, kind = "L'Ecuyer-CMRG")
# MapSol <- foreach(i = 1:4) %dopar% {
#   res <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = SigPrior, K = 20, a = 1.1, alpha = 1.1)
#   res
# }
# lapply(MapSol, function(x) ncol(x$Signatures))
# # Use consensus to determine which signatures we should use in the prior
# matches <- lapply(MapSol, function(x) {
#   match <- match_to_RefSigs(x$Signatures, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37 )
#   sort(match$best[match$cosine > 0.85])
# })
# sig_names <- sort(unique(unlist(matches)))


set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = 10,
                      method = "map",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Prost_MAP.rds.gzip", compress = "gzip")


##########################################
# Standard Poisson-Process factorization
##########################################

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = 10,
                      method = "map",
                      compressType = "sig",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))

saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Prost_MAP_sig.rds.gzip", compress = "gzip")


##########################################
# Standard Poisson-Process factorization
##########################################

set.seed(10)
out <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                      df_areas = TumorData$df_areas,
                      K = 10,
                      method = "mle",
                      controls = SigPoisProcess.controls(maxiter = maxiter, merge_move = FALSE))
saveRDS(object = out, file = "~/SigPoisProcess/output/ICGC/Prost_MLE.rds.gzip", compress = "gzip")




outMAP <- readRDS("~/SigPoisProcess/output/ICGC/Prost_MAP.rds.gzip")
plot(outMAP)


outMLE <- readRDS("~/SigPoisProcess/output/ICGC/Prost_MLE.rds.gzip")
plot(outMLE)




















