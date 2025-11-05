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
#------------------------------------------------------------
# Chain file URL
chain_url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz"

# Download and decompress
chain_file <- tempfile(fileext = ".chain.gz")
download.file(chain_url, chain_file)
chain_unzipped <- R.utils::gunzip(chain_file, remove = FALSE, temporary = TRUE)

# Import the chain
chain <- import.chain(chain_unzipped)
#--------------------------------------------------------------

devtools::document()
url <- "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak/"
page <- read_html(url)

# Extract file names from all <a> tags
file_names <- page %>%
  html_elements("a") %>%
  html_attr("href")

file_names <- file_names[grepl("E118", file_names)]
file_names <- file_names[!grepl(".tbi", file_names)]

signals <- c("DNase", "H2A.Z", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me3", "H3K27ac",
           "H3K9ac", "H3K9me3", "H3K36me3", "H3K27me3", "H4K20me1")

file_names[sapply(signals , function(x) which(grepl(x, file_names)))]


# Very useful function
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


# Download narrowPeaks now
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Liver-HCC_snp.rds.gzip")

fs <- list.files("~/SigPoisProcess/data/Roadmap/narrowPeaks/")[-1]
main_dir <- "~/SigPoisProcess/data/Roadmap/narrowPeaks/"

totals <- c()
for(i in 1:length(fs)){

  # Load file
  df_Pk <- read_tsv(paste0(main_dir, fs[i]), col_names = FALSE)
  gr_covariate <- GRanges(seqnames = df_Pk$X1,
                          IRanges(start = df_Pk$X2, end = df_Pk$X3))
  gr_covariate <- make_scored_genome(gr_covariate)
  score <- rep(1e-18, length(gr_tumor))
  overlaps <- findOverlaps(gr_tumor, gr_covariate)
  gr_tumor$score <- gr_covariate$score[subjectHits(overlaps)]

  # Append
  nm <- gsub(".imputed.narrowPeak.bed.nPk.gz", "", gsub("E118-", "", fs[i]))
  colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- nm

  # Calculate Area
  totals <- c(totals, sum(width(gr_covariate) * gr_covariate$score))
  names(totals)[length(totals)] <- nm
}

# Merge CTCF
std_chrs <- paste0("chr", c(1:22, "X", "Y"))
gr_ctcf <- import("~/ENCFF757EKU.bigBed")
# LiftOver
genome(gr_ctcf) <- "hg38"
gr_ctcf <- keepSeqlevels(unlist(liftOver(gr_ctcf, chain)), std_chrs, pruning.mode = "coarse")
genome(gr_ctcf) <- "hg19"
# Merge with covariate
gr_covariate <- make_scored_genome(gr_ctcf)
score <- rep(1e-18, length(gr_tumor))
overlaps <- findOverlaps(gr_tumor, gr_covariate)
gr_tumor$score <- gr_covariate$score[subjectHits(overlaps)]
colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- "CTCF"
# Calculate Area
totals <- c(totals, sum(width(gr_covariate) * gr_covariate$score))
names(totals)[length(totals)] <- "CTCF"

# Merge GC content
gr_gc <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig")
scores <- rep(1e-18, length(gr_tumor))
overlaps <- findOverlaps(gr_tumor, gr_gc)
scores[queryHits(overlaps)] <- gr_gc[subjectHits(overlaps)]$score
gr_tumor$score <- scores
colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- "GC"
totals <- c(totals, sum(width(gr_gc) * gr_gc$score))
names(totals)[length(totals)] <- "GC"

# Merge RepliTiming
gr_repli <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig")
sc <- gr_repli$score
sc_unif <- 1 - (sc - min(sc))/(max(sc) - min(sc))
gr_repli$score <- sc_unif
scores <- rep(1e-18, length(gr_tumor))
overlaps <- findOverlaps(gr_tumor, gr_repli)
scores[queryHits(overlaps)] <- gr_repli[subjectHits(overlaps)]$score
gr_tumor$score <- scores
colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- "RepliTiming"
totals <- c(totals, sum(width(gr_repli) * gr_repli$score))
names(totals)[length(totals)] <- "RepliTiming"

# Add Methylation
gr_methyl <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_Methylation.bigWig")
gr_methyl$score <- gr_methyl$score/mean(gr_methyl$score)
scores <- rep(1e-18, length(gr_tumor))
overlaps <- findOverlaps(gr_tumor, gr_methyl)
scores[queryHits(overlaps)] <- gr_methyl[subjectHits(overlaps)]$score
gr_tumor$score <- scores
colnames(mcols(gr_tumor))[ncol(mcols(gr_tumor))] <- "Methyl"
totals <- c(totals, sum(width(gr_methyl) * gr_methyl$score))
names(totals)[length(totals)] <- "Methyl"

# Create Areas dataset
df_areas <- data.frame(sample = unique(gr_tumor$sample),
                       t(totals))

gr_tumor_expr <- NULL
df_areas_expr <- data.frame()
for(i in 1:length(patients)){
  print(i)
  # Look at patient i
  gr_tmp <- gr_tumor[gr_tumor$sample == patients[i]]
  scores <- rep(1e-18, length(gr_tmp))
  # Load gene expression and find overlaps
  tmp_exp <- rtracklayer::import(paste0("~/SigPoisProcess/data/ICGC_GeneExpr_bigWig/", patients[i], ".bigWig"))
  tmp_exp$score <- (tmp_exp$score)/mean(tmp_exp$score)
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

df_areas_expr_merged <- df_areas %>%
  left_join(df_areas_expr, by = "sample") %>%
  drop_na()


set.seed(10)
sp <- sample(unique(gr_tumor$sample), 30)

gr_tumor_red <- gr_tumor[gr_tumor$sample%in% sp]
df_areas_red <- df_areas[df_areas$sample%in% sp, ]
rownames(df_areas_red) <- NULL

MutMatrix <- getTotalMutations(gr_tumor_red, df_areas = df_areas_red)
resMap <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 20, a = 1, alpha = 1)
plot_96SBSsig(resMap$Signatures)

set.seed(42)
out_map1 <- SigPoisProcess(gr_Mutations = gr_tumor_expr,
                      df_areas = df_areas_expr_merged,
                      K = 10,
                      method = "map",
                      compressType = "sig_cov",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 500,
                                                         merge_move = FALSE, tol = 1e-6))
plot(out_map1)

set.seed(42)
out_map2 <- SigPoisProcess(gr_Mutations = gr_tumor_expr,
                      df_areas = df_areas_expr_merged,
                      K = 10,
                      method = "map",
                      compressType = "sig",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 500,
                                                         merge_move = FALSE, tol = 1e-6))
plot(out_map2)

set.seed(42)
out_mle <- SigPoisProcess(gr_Mutations = gr_tumor_expr,
                       df_areas = df_areas_expr_merged,
                       K = 10,
                       method = "mle",
                       compressType = "sig",
                       prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                       controls = SigPoisProcess.controls(maxiter = 500,
                                                          merge_move = FALSE, tol = 1e-6))
plot(out_mle)

set.seed(10)
out <- SigPoisProcess(gr_Mutations = gr_tumor_red,
                      df_areas = df_areas_red,
                      K = 8,
                      method = "mle",
                      compressType = "sig",
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1, a0 = 11, b0 = 0.01 * 10),
                      controls = SigPoisProcess.controls(maxiter = 1000,
                                                         merge_move = FALSE, tol = 1e-5))

plot(out)
round((out$Betas * out$Xbar)[, 1,], 3)


list_gene_expr <- list.files("~/SigPoisProcess/data/ICGC_GeneExpr_bigWig/")
patients <- gsub(".bigWig", "", list_gene_expr)

patients <- unique(gr_tumor$sample)[unique(gr_tumor$sample) %in% patients]


MutMatrix <- getTotalMutations(gr_tumor_expr, df_areas = df_areas_expr_merged)
set.seed(10)
MapSol_final <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 20, a = 1.1, alpha = 1.1)
plot_96SBSsig(res$Signatures)

# Starting point
K <- ncol(MapSol_final$Signatures)
J <- ncol(MapSol_final$Theta)
p <- ncol(mcols(gr_tumor_expr)[-c(1,2,3)])
R_start <- as.matrix(MapSol_final$Signatures)
Xtotals <- as.matrix(df_areas_expr_merged[, -1])
Betas_start <- array(NA, dim = c(K, J, p))
for(i in 1:p){
  Betas_start[, , i] <- t(apply(MapSol_final$Theta, 1, function(x) x/Xtotals[, i]))
}
Mu_start <- matrix(rowMeans(MapSol_final$Theta))[, rep(1, p)]



set.seed(42)
out2 <- SigPoisProcess(gr_Mutations = gr_tumor_expr,
                      df_areas = df_areas_expr_merged,
                      K = ncol(R_start),
                      method = "map",
                      compressType = "sig",
                      init =  SigPoisProcess.init(R_start = R_start,
                                                  Betas_start = Betas_start,
                                                  Mu_start = Mu_start),
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1, a0 =11, b0 =  0.01 * 10),
                      controls = SigPoisProcess.controls(maxiter = 1000,
                                                         merge_move = FALSE, tol = 1e-5))

plot(out)
plot(out2)

set.seed(10)
out_gibbs <- SigPoisProcess(gr_Mutations = gr_tumor_expr,
                            df_areas = df_areas_expr_merged,
                            K = ncol(R_start),
                            method = "gibbs",
                            compressType = "sig",
                            init =  SigPoisProcess.init(R_start = out2$Signatures,
                                                        Betas_start = out2$Betas,
                                                        Mu_start = out2$Mu),
                            prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1, a0 =11, b0 =  0.01 * 10),
                            controls = SigPoisProcess.controls(maxiter = 1000,
                                                               merge_move = FALSE, tol = 1e-5,
                                                               nsamples = 500, burnin = 500))

plot(out_gibbs)
out_gibbs$gibbs$MU
match_to_RefSigs(out$Signatures, RefSigs)



RefSigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
match_to_RefSigs(out$Signatures, RefSigs)
# Merge with old areas
df_areas_expr <- df_areas %>%
  left_join(df_areas_expr, by = "sample") %>%
  drop_na()



Probs <- compute_SignatureCovariate_probs(out, gr_tumor_red)
round(Probs[10, , ], 4)
round(apply(Probs, c(2,3), sum), 5)

mltools::mcc(gr_tumor_red$H4K20me1, gr_tumor_red$H2A.Z)
mltools::mcc(gr_tumor_red$H3K4me1, gr_tumor_red$H3K27ac)
table(gr_tumor_red$H3K4me1, gr_tumor_red$H3K27ac)





gr_gc <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig")


gr_gc$score[gr_gc$score > 0]
hist(gr_gc$score)


sc <- gr_gc$score
sc2 <- (sc - mean(sc[sc > 0]))/sd(sc[sc > 0])
sc3 <- sc2[5000:10000]

plot(pmax(sc3, 0), type = "l")
lines(pmax(-sc3, 0), col = "red")
plot(sc3, type = "l")
abline(h = 0)

gr_ctcf <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/CTCF_histones_1kb/ENCFF094CYN.bigWig")

hist(asinh(gr_ctcf$score))
hist(gr_ctcf$score)

sc <- asinh(gr_ctcf$score)
sc2 <- (sc - median(sc))/sd(sc)
sc3 <- sc2[1:1000]
plot(sc3, type = "l")
plot(pmax(sc3, 0), type = "l")
lines(pmax(-sc3, 0), col = "red")




# Take 30 patients
set.seed(10)
pats <- sample(size = 30, unique(gr_tumor$sample))

gr_tumor_red <- gr_tumor[gr_tumor$sample %in% pats]
totals <- c()

# CTCF
gr_tmp <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_female_CTCF.bigWig")
sc <- asinh(gr_tmp$score)
gr_tmp$score <- (sc - mean(sc))/sd(sc)
gr_tmp$score_rich <- pmax(0, gr_tmp$score)
gr_tmp$score_depl <- pmax(0, -gr_tmp$score)
score_rich <- score_depl <- rep(0, length(gr_tumor_red))
overlaps <- findOverlaps(gr_tumor_red, gr_tmp)
score_rich[queryHits(overlaps)] <- gr_tmp[subjectHits(overlaps)]$score_rich
score_depl[queryHits(overlaps)] <- gr_tmp[subjectHits(overlaps)]$score_depl

gr_tumor_red$score_rich <- score_rich
gr_tumor_red$score_depl <- score_depl
col <- ncol(mcols(gr_tumor_red))
colnames(mcols(gr_tumor_red))[c(col - 1, col)] <- paste0("CTCF", c("_rich", "_depl"))

totals <- c(totals, c(sum(width(gr_tmp) * gr_tmp$score_rich),
                      sum(width(gr_tmp) * gr_tmp$score_depl)))
names(totals)[c(length(totals) - 1, length(totals))] <- paste0("CTCF", c("_rich", "_depl"))

# H3K9me3
gr_tmp <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_female_H3K9me3.bigWig")
sc <- asinh(gr_tmp$score)
gr_tmp$score <- (sc - mean(sc))/sd(sc)
gr_tmp$score_rich <- pmax(0, gr_tmp$score)
gr_tmp$score_depl <- pmax(0, -gr_tmp$score)
score_rich <- score_depl <- rep(0, length(gr_tumor_red))
overlaps <- findOverlaps(gr_tumor_red, gr_tmp)
score_rich[queryHits(overlaps)] <- gr_tmp[subjectHits(overlaps)]$score_rich
score_depl[queryHits(overlaps)] <- gr_tmp[subjectHits(overlaps)]$score_depl

gr_tumor_red$score_rich <- score_rich
gr_tumor_red$score_depl <- score_depl
col <- ncol(mcols(gr_tumor_red))
colnames(mcols(gr_tumor_red))[c(col - 1, col)] <- paste0("H3K9me3", c("_rich", "_depl"))

totals <- c(totals, c(sum(width(gr_tmp) * gr_tmp$score_rich),
                      sum(width(gr_tmp) * gr_tmp$score_depl)))
names(totals)[c(length(totals) - 1, length(totals))] <- paste0("H3K9me3", c("_rich", "_depl"))

# GC content
gr_tmp <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/gc_content/gc_content_1kb.bigWig")
when_zero <- gr_tmp$score > 0
sc <- gr_tmp$score[when_zero]
gr_tmp$score[when_zero] <- (sc - mean(sc))/sd(sc)
gr_tmp$score_rich <- pmax(0, gr_tmp$score)
gr_tmp$score_depl <- pmax(0, -gr_tmp$score)
score_rich <- score_depl <- rep(0, length(gr_tumor_red))
overlaps <- findOverlaps(gr_tumor_red, gr_tmp)
score_rich[queryHits(overlaps)] <- gr_tmp[subjectHits(overlaps)]$score_rich
score_depl[queryHits(overlaps)] <- gr_tmp[subjectHits(overlaps)]$score_depl

gr_tumor_red$score_rich <- score_rich
gr_tumor_red$score_depl <- score_depl
col <- ncol(mcols(gr_tumor_red))
colnames(mcols(gr_tumor_red))[c(col - 1, col)] <- paste0("GC", c("_rich", "_depl"))

totals <- c(totals, c(sum(width(gr_tmp) * gr_tmp$score_rich),
                      sum(width(gr_tmp) * gr_tmp$score_depl)))
names(totals)[c(length(totals) - 1, length(totals))] <- paste0("GC", c("_rich", "_depl"))

# Areas
df_areas_red <- data.frame(sample = unique(gr_tumor_red$sample),
                       t(totals))

set.seed(42)
out_map1 <- SigPoisProcess(gr_Mutations = gr_tumor_red,
                           df_areas = df_areas_red,
                           K = 10,
                           method = "map",
                           compressType = "sig_cov",
                           prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                           controls = SigPoisProcess.controls(maxiter = 4000,
                                                              merge_move = FALSE, tol = 1e-5))
match_to_RefSigs(out_map1$Signatures, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)
plot(out_map1)


