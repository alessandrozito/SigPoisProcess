# Let's try chormatin states
library(SigPoisProcess)
library(rtracklayer)

std_chrs <- paste0("chr", c(1:22, "X", "Y"))

gr_tumor <- readRDS("data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")

# Load chromatin states
chromStates <- rtracklayer::import("~/E110_15_coreMarks_dense.bed")
chromStates <- keepSeqlevels(chromStates, std_chrs, pruning.mode = "coarse")

# Reshape the matrix
df_states <- as.data.frame(chromStates)
df_states$name <- as.factor(df_states$name)
mm <- model.matrix(start ~ name - 1, data = df_states) + 1e-18
colnames(mm) <- sub("^[^_]+_", "", colnames(mm))

# Merge it with gr
mcols(chromStates) <- mm

# Merge it with data on tumors
overlaps <- findOverlaps(gr_tumor, chromStates)
mcols(gr_tumor) = cbind(mcols(gr_tumor), mcols(chromStates[subjectHits(overlaps)]))

# Calculate total areas
patients <- unique(gr_tumor$sample)
J <- length(patients)
df_areas <- as.data.frame(crossprod(width(chromStates), mm)[rep(1, J), ])
df_areas$sample <- patients

# Run the method
MutMatrix <- getTotalMutations(gr_tumor)
SigPrior <- 2000 * CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37 + 1
set.seed(10)
resMap <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = SigPrior, K = 10, a = 1.1, alpha = 1.1)
CompressiveNMF::plot_SBS_signature(resMap$Signatures)
match <- match_to_RefSigs(resMap$Signatures, CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)

SigPrior_map <- SigPrior[, sort(unique(match$best))]
set.seed(42)
outComp <- SigPoisProcess(gr_Mutations = gr_tumor,
                      df_areas = df_areas,
                      add_baseline = FALSE,
                      K = 5,
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1,
                                                                SigPrior = SigPrior_map),
                      controls = SigPoisProcess.controls(maxiter = 5000, tol = 1e-5))
plot(outComp)
round(outComp$Betas * outComp$Xbar, 5)

match_to_RefSigs(outComp$Signatures[, rowMeans(outComp$Mu)>0.01], CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)

set.seed(42)
out <- SigPoisProcess(gr_Mutations = gr_tumor,
                      df_areas = df_areas,
                      add_baseline = FALSE,
                      method = "mle",
                      K = 12,
                      prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1,
                                                                SigPrior = SigPrior_map),
                      controls = SigPoisProcess.controls(maxiter = 5000))

round(out$Betas * out$Xbar, 5)
plot(out)
log10(10)

corrplot(cor())

colSums(as.matrix(mcols(gr_tumor)[, -c(1,2,3)]))

table(gr_tumor$TssA)

gr_tumor2 <- gr_tumor
mcols(gr_tumor2) <- mcols(gr_tumor2)[, colnames(mcols(gr_tumor2)) != "Quies"]
df_areas2 <- df_areas[, colnames(df_areas) != "Quies"]

set.seed(42)
outComp <- SigPoisProcess(gr_Mutations = gr_tumor2,
                          df_areas = df_areas2,
                          add_baseline = TRUE,
                          K = 10,
                          prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                          controls = SigPoisProcess.controls(maxiter = 1000, tol = 1e-5))
plot(outComp)

set.seed(42)
outComp <- SigPoisProcess(gr_Mutations = gr_tumor,
                          df_areas = df_areas,
                          add_baseline = FALSE,
                          K = 10,
                          prior_params = SigPoisProcess.PriorParams(a = 1.1, alpha = 1.1),
                          controls = SigPoisProcess.controls(maxiter = 1000, tol = 1e-5))
plot(outComp)


mcols(gr_tumor)xs




