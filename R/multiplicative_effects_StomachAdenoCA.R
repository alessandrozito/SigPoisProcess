# Multiplicative effects on Stomach-AdenoCA
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

genome_5kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 5000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

add_bin_weights <- function(gr){
  # Import blacklist and regions with gaps
  blacklist <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/hg19-blacklist.v2.bed")
  gap_gr <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/gaps_hg19.bed")

  # Adjust bin sizes to remove blacklist and gaps
  BinWeight <- width(gr)
  score <- gr$score
  # Remove gaps
  overlapGaps <- findOverlaps(gr, gap_gr)
  rangesGaps <- gr[queryHits(overlapGaps)]
  # Find the number of ACGT in the gap regions, and adjust the bin sizes to remove
  # the NNNNN sequences
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rangesGaps)
  countACGT <- rowSums((Biostrings::alphabetFrequency(seqs))[, c("A", "C", "G", "T")])
  BinWeight[queryHits(overlapGaps)] <- countACGT
  # Exlude now the blacklist as well
  overlap_blackList <- findOverlaps(gr, blacklist)
  pairs <- Pairs(gr[queryHits(overlap_blackList)], blacklist[subjectHits(overlap_blackList)])
  grIntersect <- pintersect(pairs)
  BinWeight[queryHits(overlap_blackList)] <- pmax(BinWeight[queryHits(overlap_blackList)] - width(grIntersect), 0)
  score[BinWeight == 0] <- 0 # Replace it with 0, knowing that we will remove it afterwards

  gr$score <- score
  gr$bin_weight <- BinWeight
  return(gr)
}

# Load tumor data
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")
# Remove mutations from blacklist
blacklist <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/hg19-blacklist.v2.bed")
overlaps <- findOverlaps(gr_tumor, blacklist)
gr_tumor <- gr_tumor[-queryHits(overlaps)]
gr_tumor <- gr_tumor[seqnames(gr_tumor) != "chrY"]

#---- Load Methylation covariate
gr_methyl <- import("~/SigPoisProcess/data/ENCODE_pvalues/old_data/Stomach_AdenoCA/methylation_bigWig_GRCh38/Stomach_methylation_male_hg19_1kb.bigWig")
gr_methyl <- binnedAverage(bins = genome_5kb,
                              numvar = coverage(gr_methyl, weight = "score"),
                              varname = "score")
gr_methyl <- add_bin_weights(gr_methyl)
gr_methyl <- gr_methyl[seqnames(gr_methyl) != "chrY"]

sc <- gr_methyl$score[gr_methyl$bin_weight > 0]
gr_methyl$score2 <- 0
sc <- pmin(sc, quantile(sc, 0.999))
gr_methyl$score2[gr_methyl$bin_weight > 0] <- (sc - min(sc))/(max(sc) - min(sc))

#---- Load replication timing
gr_repli <- import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig")
gr_repli <- binnedAverage(bins = keepSeqlevels(genome_5kb, std_chrs[-24], pruning.mode = "coarse"),
                           numvar = coverage(gr_repli, weight = "score"),
                           varname = "score")
gr_repli <- add_bin_weights(gr_repli)
sc <- gr_repli$score[gr_repli$bin_weight > 0]
gr_repli$score2 <- 0
gr_repli$score2[gr_repli$bin_weight > 0] <- pmin(sc, quantile(sc, 0.999))/100

#---- Load H3K9me3
gr_H3K9me3 <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K9me3_1kb.bigWig")
gr_H3K9me3 <- binnedAverage(bins = keepSeqlevels(genome_5kb, std_chrs, pruning.mode = "coarse"),
                          numvar = coverage(gr_H3K9me3, weight = "score"),
                          varname = "score")
gr_H3K9me3 <- gr_H3K9me3[seqnames(gr_H3K9me3) != "chrY"]
gr_H3K9me3 <- add_bin_weights(gr_H3K9me3)
sc <- gr_H3K9me3$score[gr_H3K9me3$bin_weight > 0]
sc <- pmin(sc, quantile(sc, 0.999))
gr_H3K9me3$score2 <- 0
gr_H3K9me3$score2[gr_H3K9me3$bin_weight > 0] <- (sc - min(sc))/(max(sc) - min(sc))

#---- Load CTCF
gr_CTCF <- import("../data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_CTCF_1kb.bigWig")
gr_CTCF <- binnedAverage(bins = keepSeqlevels(genome_5kb, std_chrs, pruning.mode = "coarse"),
                            numvar = coverage(gr_CTCF, weight = "score"),
                            varname = "score")
gr_CTCF <- gr_CTCF[seqnames(gr_CTCF) != "chrY"]
gr_CTCF <- add_bin_weights(gr_CTCF)
sc <- gr_CTCF$score[gr_CTCF$bin_weight > 0]
sc <- pmin(sc, quantile(sc, 0.999))
gr_CTCF$score2 <- 0
gr_CTCF$score2[gr_CTCF$bin_weight > 0] <- (sc - min(sc))/(max(sc) - min(sc))


#---- Load H3K36me3
gr_H3K36me3 <- import("../data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K36me3_1kb.bigWig")
gr_H3K36me3 <- binnedAverage(bins = keepSeqlevels(genome_5kb, std_chrs, pruning.mode = "coarse"),
                         numvar = coverage(gr_H3K36me3, weight = "score"),
                         varname = "score")
gr_H3K36me3 <- gr_H3K36me3[seqnames(gr_H3K36me3) != "chrY"]
gr_H3K36me3 <- add_bin_weights(gr_H3K36me3)
sc <- gr_H3K36me3$score[gr_H3K36me3$bin_weight > 0]
sc <- pmin(sc, quantile(sc, 0.999))
gr_H3K36me3$score2 <- 0
gr_H3K36me3$score2[gr_H3K36me3$bin_weight > 0] <- (sc - min(sc))/(max(sc) - min(sc))



# Merge with patients
overlaps_methyl <- findOverlaps(gr_tumor, gr_methyl)
gr_tumor$Methyl <- gr_methyl[subjectHits(overlaps_methyl)]$score2
overlaps_repli <- findOverlaps(gr_tumor, gr_repli)
gr_tumor$RepliTime <- gr_repli[subjectHits(overlaps_repli)]$score2
overlaps_H3K9me3 <- findOverlaps(gr_tumor, gr_H3K9me3)
gr_tumor$H3K9me3  <- gr_H3K9me3[subjectHits(overlaps_H3K9me3)]$score2
overlaps_CTCF <- findOverlaps(gr_tumor, gr_CTCF)
gr_tumor$CTCF  <- gr_CTCF[subjectHits(overlaps_CTCF)]$score2
overlaps_H3K36me3 <- findOverlaps(gr_tumor, gr_H3K36me3)
gr_tumor$H3K36me3  <- gr_H3K36me3[subjectHits(overlaps_H3K36me3)]$score2


# Model paramters
SignalTrack <- cbind("Methyl" = gr_methyl$score2,
                     "RepliTime" = gr_repli$score2,
                     "H3K9me3" = gr_H3K9me3$score2,
                     "CTCF" = gr_CTCF$score2,
                     "H3K36me3" = gr_H3K36me3$score2)

bin_weight <- gr_methyl$bin_weight



#----------------------------------- Run the model from random starting point
set.seed(10)
out <- SigPoisProcess_mult(gr_Mutations = gr_tumor,
                           SignalTrack = SignalTrack,
                           bin_weight = bin_weight,
                           K = 12,
                           controls = SigPoisProcess_mult.controls(maxiter = 100))

out$Betas
out$Thetas
out$Mu

# Start from a matrix of known signatures, and compare it against the
# Fixed-case solution
Sigs_to_include <- c("SBS1", "SBS2", "SBS3", "SBS5","SBS13", "SBS17a",
                     "SBS17b", "SBS18", "SBS21", "SBS22a", "SBS26", "SBS44")
MutMatrix <- getTotalMutations(gr_Mutations = gr_tumor)
Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
SigPrior <- 1e6 * Sigs[, Sigs_to_include] + 1
startSol <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = SigPrior, K = 0, a = 1.1, alpha = 1.1)

# Now, run the model with fixed signatures.
set.seed(10)
out <- SigPoisProcess_mult(gr_Mutations = gr_tumor,
                           SignalTrack = SignalTrack,
                           bin_weight = bin_weight,
                           K = 12,
                           controls = SigPoisProcess_mult.controls(maxiter = 100),
                           update_Signatures = FALSE,
                           init = SigPoisProcess_mult.init(R_start = startSol$Signatures))




CompressiveNMF::plot_SBS_signature(out$Signatures)

sol <- SigPoisProcess_mult(gr_tumor)


X <- as.matrix(mcols(gr_tumor)[, -c(1,2,3)])

gr_tumor$sample <- as.factor(gr_tumor$sample)
channel_id = as.numeric(gr_tumor$channel) - 1
sample_id = as.numeric(gr_tumor$sample) - 1

Sigs_to_include <- c("SBS1", "SBS2", "SBS3", "SBS5","SBS13", "SBS17a",
                     "SBS17b", "SBS18", "SBS26", "SBS44")
# Starting point
MutMatrix <- getTotalMutations(gr_Mutations = gr_tumor)
Sigs <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
SigPrior <- 3000 * Sigs + 1
startSol <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = SigPrior, K = 20, a = 1.1, alpha = 1.1)

set.seed(42)
J <- length(unique(sample_id))
K <- ncol(startSol$Signatures)
R_start <- startSol$Signatures #sample_signatures(matrix(1, 96, K))
Theta_start <- startSol$Theta/sum(bin_weight)

set.seed(10)
R_start <- Sigs[, Sigs_to_include]#sample_signatures(matrix(1, 96, K))
K <- ncol(R_start)
Theta_start <- apply(MutMatrix, 2, function(x) nnls::nnls(b = x, A = R_start)$x + 1)/sum(bin_weight)
#Theta_start <- matrix(rgamma(K * J, 1,1)/sum(bin_weight), nrow = K)
Betas_start <- matrix(0, nrow = ncol(X), ncol = K)
SigPrior <- matrix(1.1, 96, K)
#Betas_start <- matrix(rnorm(p * K, sd = 2), nrow = p)
sol <- PPF_loglinear(R_start = R_start,
                     Theta_start = Theta_start,
                     Betas_start = Betas_start,
                     X = X,
                     SignalTrack = SignalTrack,
                     bin_weight = bin_weight,
                     channel_id = channel_id,
                     sample_id = sample_id,
                     SigPrior = SigPrior,
                     a = 1.1, a0 = 1.1 * J + 1, b0 = 1.1 * J * 0.01,
                     update_R = TRUE,
                     update_Theta = TRUE,
                     update_Betas = TRUE,
                     method = "map",
                     maxiter = 1000,
                     n_iter_betas = 2, tol = 1e-7)

sol$Mu
rownames(sol$Betas) <- colnames(X)
rownames(sol$R) <- levels(gr_tumor$channel)
#colnames(sol$R)[match_id] <- paste0("Sig", 1:ncol(sol$R))
colnames(sol$R) <- colnames(sol$Mu) <- colnames(sol$Betas) <- paste0("Sig", 1:ncol(sol$R))
SigPoisProcess::plot_96SBSsig(sol$R)
plot_96SBSsig(R_start)


match_to_RefSigs(sol$R, Sigs)
round(t(sol$Betas), 3)


tots <- colSums(bin_weight * exp(as.matrix(SignalTrack %*% sol$Betas)))
plot(sol$R %*% (sol$Theta * matrix(tots)[, rep(1, J)]), MutMatrix)

hist(exp(as.matrix(SignalTrack %*% sol$Betas)))

# Calculate the signatures
channel_id = as.numeric(gr_tumor$channel) - 1
sample_id = as.numeric(gr_tumor$sample) - 1
Probs_un <- sol$R[gr_tumor$channel, ] * t(sol$Theta[, gr_tumor$sample]) * exp(X %*% sol$Betas)
Probs <- Probs_un/rowSums(Probs_un)
W <- sample_Categorical(Probs)

pp <- sample_R(sol$R, W, SigPrior, channel_id)
dimnames(pp) <- dimnames(sol$R)
plot_96SBSsig(pp)

tt <- sample_Theta(t(sol$Theta), sol$Betas, SignalTrack, bin_weight,
                   W = W, sol$Mu, a = 1.1, sample_id= sample_id)
plot(tt, t(sol$Theta))


mm <- sample_Mu(t(sol$Theta), sol$Betas, a = 1.1, a0 = 1.1 * J + 1, b0 = 1.1 * J * 0.01,
                sum(bin_weight))
sol$Mu


sol <- .PoissonProcess_optim(R_start = R_start,
                     Theta_start = Theta_start,
                     Betas_start = Betas_start,
                     X = X,
                     SignalTrack = SignalTrack,
                     bin_weight = bin_weight,
                     channel_id = channel_id,
                     sample_id = sample_id,
                     SigPrior = SigPrior,
                     a = 1.1, a0 = 1.1 * J + 1, b0 = 1.1 * J * 0.01,
                     update_R = TRUE,
                     update_Theta = TRUE,
                     update_Betas = TRUE,
                     method = "map",
                     maxiter = 500,
                     n_iter_betas = 2, tol = 1e-7)

plot_96SBSsig(sol$R)
sol$Betas
# Relationship between risoluzione e  e

gr_H3K9me3_2 <- import("~/SigPoisProcess/data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K9me3_1kb.bigWig")
gr_H3K9me3 <- binnedAverage(bins = keepSeqlevels(genome_5kb, std_chrs, pruning.mode = "coarse"),
                            numvar = coverage(gr_H3K9me3, weight = "score"),
                            varname = "score")

gr_H3K36me3_2<- import("../data/ENCODE_foldchange/Cancer_covariates/Stomach-AdenoCA_male_H3K36me3_1kb.bigWig")
gr_H3K36me3 <- binnedAverage(bins = keepSeqlevels(genome_5kb, std_chrs, pruning.mode = "coarse"),
                             numvar = coverage(gr_H3K36me3, weight = "score"),
                             varname = "score")
k <- 1
sample_Betas(W[, k], X, SignalTrack, bin_weight,
             sol$Mu[k], sum(sol$Theta[k, ]))


set.seed(10)
id <- sort(sample(1:length(gr_H3K9me3_2), 10000))
plot(gr_H3K9me3_2[id]$score, gr_H3K36me3_2[id]$score)

Beta_old <-  sol$Betas[, 1]
Beta_prop <- sol$Betas[, 1] + sqrt(0.0001) * rnorm(5)




run_adaptive <- function(beta, burn_in = 200, nsamples = 200) {
  epsilon <- 1e-6
  p <- length(beta)
  out <- matrix(0, nsamples, p)
  S <- diag(epsilon, p)
  Sigma_r <- diag(0, p)
  mu_r <- beta
  pb <- txtProgressBar(style=3)
  for (r in 1:(burn_in + nsamples)) {
    setTxtProgressBar(pb, r/(nsamples + burn_in))
    # Updating the covariance matrix
    if(r > 1){
      Sigma_r <- (r - 2) / (r - 1) * Sigma_r + tcrossprod(beta - mu_r) / r
      mu_r <- (r - 1) / r * mu_r + beta / r
      S <- 2.38^2 * Sigma_r / p + diag(epsilon, p)
    }

    # Eigen-decomposition
    eig <- eigen(S, symmetric = TRUE)
    A1 <- t(eig$vectors) * sqrt(eig$values)

    beta_new <- beta + c(matrix(rnorm(p), 1, p) %*% A1)

    new <- - sum(sol$Theta[k, ]) * bin_weight %*% exp(SignalTrack %*% beta_new) +
      W[, k] %*% X %*% beta_new +
      sum(dnorm(beta_new, sd = sqrt(sol$Mu[k]), log = TRUE))

    old <- - sum(sol$Theta[k, ]) * bin_weight %*% exp(SignalTrack %*% beta) +
      W[, k] %*% X %*% beta +
      sum(dnorm(beta, sd = sqrt(sol$Mu[k]), log = TRUE))

    if(new - old > log(runif(1))){
      beta <- beta_new
    }

    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- beta
    }
  }
  close(pb)
  out
}



H <- sum(sol$Theta[k, ]) * crossprod(SignalTrack * c(exp(SignalTrack %*% sol$Betas[, k])),SignalTrack)
solve(H)


set.seed(10)
k <- 2
tmp <- run_adaptive(sol$Betas[, k], burn_in = 1000, nsamples = 3000)
hist(tmp[,3])
coda::effectiveSize(tmp)
plot(tmp[, 1], type = "l")
plot(tmp[, 2], type = "l")
plot(tmp[, 3], type = "l")
plot(tmp[, 4], type = "l")
plot(tmp[, 5], type = "l")

acf(tmp[, 5])

Probs[, 1]
S <- crossprod(X * Probs_un[, 1], X)



##### Let's make a cool plot now.
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

# Define parameters
chr <- "chr1"
region_start <- 1e6     # 1,000,000
region_end   <- 1.5e6   # 1,500,000
bin_size     <- 5000    # Tunable: change this as needed

# Create a GRanges for the region of interest
roi <- GRanges(chr, IRanges(region_start, region_end))

# Calculate bin starts and ends
bin_starts <- seq(region_start, region_end, by=bin_size)
bin_ends   <- pmin(bin_starts + bin_size - 1, region_end)

# Create a GRanges object with all bins
bins <- GRanges(seqnames = chr, ranges = IRanges(start=bin_starts, end=bin_ends))

# Inspect the first few bins
tilewidth <- 5000
genome_5kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                         tilewidth = 5000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X"))


# Model paramters
SignalTrack <- cbind("Methyl" = gr_methyl$score2,
                     "RepliTime" = gr_repli$score2,
                     "H3K9me3" = gr_H3K9me3$score2,
                     "CTCF" = gr_CTCF$score2,
                     "H3K36me3" = gr_H3K36me3$score2)

gr_SignalTrack <- granges(gr_methyl)
gr_SignalTrack$bin_weight <- gr_methyl$bin_weight
mcols(gr_SignalTrack) <- cbind(mcols(gr_SignalTrack), SignalTrack)


# Find overlaps with current gr_SignalTrack
overlaps <- findOverlaps(gr_select, gr_SignalTrack)
gr_subset <- gr_SignalTrack[subjectHits(overlaps)]

# For each region, extract the sequence of nucleotides (be mindful of chromosome
# edges, and include the parts on both sides)

start(gr_subset) <- pmax(start(gr_subset) - 1, 1)
end(gr_subset) <- end(gr_subset) + 1

# Find sequences and their reverse complement
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, granges(gr_subset))
seqs_rc <- reverseComplement(seqs)

# Take a C>T mutation and find all possible triplets that can mutate
seq <- as.character(seqs[1]);
triplets <- str_sub(seq, 1:(str_length(seq)-2), 3:str_length(seq))

# Take the centers.
centers <- str_sub(triplets, 2, 2)

triplets_C <- triplets[centers == "C"]
subs_C <- c("C>A", "C>G", "C>T")



triplets_T <- triplets[centers == "T"]
subs_T <- c("T>A", "T>C", "T>G")

Betas <- out$Betas
Sigs <- out$Signatures
Theta <- out$Thetas

SigTheta <- t(Sigs["A[C>T]G", ] * Theta[, 3])
bigProd <- SigTheta[rep(1, nrow(ExpSolBetas)), ] * ExpSolBetas
bigProd/rowSums(bigProd)

# Subset the signalTrack into relevant regions
SignalTrack_subset <- as.matrix(mcols(gr_subset)[, -1])
ExpSolBetas <- exp(SignalTrack_subset %*% Betas)

# Extract the mutation opportunity for each triplet
motif_labels <- rownames(Sigs)
patterns <-  paste0(substr(motif_labels, 1, 1), substr(motif_labels, 3, 3), substr(motif_labels, 7, 7))
counts_mat <- t(vcountPDict(DNAStringSet(patterns), seqs, fixed=TRUE))
colnames(counts_mat) <- motif_labels
prob_mat <- counts_mat/rowSums(counts_mat)

# Find patient 1
j <- "DO218068" #"DO217814"
tmp <- array(NA, dim = c(nrow(ExpSolBetas), ncol(ExpSolBetas), I),
             dimnames = list(NULL, colnames(ExpSolBetas), rownames(Sigs)))
for(i in 1:96){
  SigTheta <- t(Sigs[i, ] * Theta[, j])
  bigProd <- SigTheta[rep(1, nrow(ExpSolBetas)), ] * ExpSolBetas
  sig_probs <-  bigProd/rowSums(bigProd)
  #sig_probs_adj <- sig_probs * prob_mat[, rep(i, ncol(sig_probs))]
  tmp[, , i] <- sig_probs#sig_probs_adj/rowSums(sig_probs_adj)
}

i <- "A[C>A]C"
df_mod <- data.frame(tmp[, , i])
colnames(df_mod) <- colnames(startSol$Signatures)

pr_const <- startSol$Signatures[i, ] * startSol$Theta[, j]
df_const <- data.frame(t(pr_const/sum(pr_const))[rep(1,nrow(tmp)), ]) %>%
  mutate(region = 1:nrow(tmp),
         method = "Standard NMF")

window <- 4
df_mod %>%
  mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
                                                   align = "center"))) %>%
  mutate(region = 1:nrow(tmp),
         method = "PoissonProcess") %>%
  bind_rows(df_const) %>%
  gather(key = "Sig", value = "Prob", - region, -method) %>%
  filter(Sig %in% head(names(sort(pr_const, decreasing = TRUE)), 4))%>%
  ggplot() +
  geom_line(aes(x = region, y = Prob, color = Sig, linetype = method)) +
  scale_color_manual(values = c("black", "darkorange", "red", "darkblue")) +
  ylab("Signature attribution probability") +
  xlab("Genomic regions in chr1:18000001:22000000 (5kb)") +
  facet_wrap(~paste0("Mutational channel ", i)) +
  theme_bw()


plot(tmp[, 1,i], type = "l", ylim = c(0,1))
lines(tmp[, 2,i], col = "red")
lines(tmp[, 3,i], col = "blue")
lines(tmp[, 4,i], col = "green")
lines(tmp[, 5,i], col = "gray")




tmp2 <- array(NA, dim = dim(tmp))
for(k in 1:12) tmp2[, k ,] <- counts_mat

tmp2 <- counts_mat_arr * tmp

# Renormalize

vcountPattern

patterns <-  paste0(substr(motif_labels, 1, 1), substr(motif_labels, 3, 3), substr(motif_labels, 7, 7))



library(Biostrings)
vcountPattern("ACA", seqs, fixed=TRUE)



# Two key things to do by tonight

# 1) Test different aggregations (2kb vs 5kb)
# 2) Test standardized vs transformed between 0 and 1















