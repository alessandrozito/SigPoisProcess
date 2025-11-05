library(tidyverse)
library(BSgenome)
library(Biostrings)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

source("~/SigPoisProcess_analysis/utils.R")

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
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Liver-HCC_merged_1kb.rds.gzip")

library(SigPoisProcess)
out_Proc_MAP <- SigPoisProcess(gr_Mutations = TumorData$gr_tumor,
                               df_areas = TumorData$df_areas,
                               method = "map",
                               prior_params = list(a = 1.1, alpha = 1.1),
                               K = 10,
                               controls = list(tol = 1e-6,
                                               maxiter = 200,
                                               merge_move = TRUE))
plot(out_Proc_MAP)


TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Skin-Melanoma_merged_1kb.rds.gzip")
TumorData$gr_tumor

Xall <- getTotalMutations(TumorData$gr_tumor)


genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
std_chrs <- paste0("chr", c(1:22, "X", "Y"))
chrom_lengths <- seqlengths(genome)[1:24]
X <- cbind("baseline" = 1, as.matrix(mcols(TumorData$gr_tumor)[, -c(1,2,3)]))
X_totals <- as.matrix(TumorData$df_areas[, -c(1,2)])
rownames(X_totals) <- TumorData$df_areas$sample
X_totals2 <- cbind("baseline" = sum(1 * chrom_lengths), X_totals)
Mutations <- as.data.frame(mcols(TumorData$gr_tumor)[, c(2,3)])
Mutations$sample <- factor(Mutations$sample, levels = rownames(X_totals))
Mutations$channel <- as.factor(Mutations$channel)

library(SigPoisProcess)

# Part 1 - Use CompressiveNMF to get a first idea
MutMatrix <- getTotalMutations(Mutations)
dim(MutMatrix)

set.seed(10)
S <- 1e6 * CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37 + 1
resMap <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 20)
CompressiveNMF::plot_SBS_signature(resMap$Signatures)
plot(resMap$Signatures %*% resMap$Theta, MutMatrix)
match_to_RefSigs(resMap$Signatures, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)

library(SigPoisProcess)
a <- 2
alpha <- 1.01
K <- 12
J <- length(unique(Mutations$sample))
a0 <- a * J - 1
b0 <- a * 0.001 * J
set.seed(10)
out_Proc_MAP <- SigPoisProcess(Mutations = Mutations,
                                X = X,
                                X_total = X_totals2,
                                method = "map",
                                prior_params = list(a = a, alpha = alpha,
                                                    a0 = a0, b0 = b0),
                                K = K,
                                controls = list(tol = 1e-6,
                                                maxiter = 500,
                                                merge_move = TRUE))

plot.SigPoisProcess(out_Proc_MAP)

K <- 12
set.seed(10)
out_Proc_MLE <- SigPoisProcess(Mutations = Mutations,
                               X = X,
                               X_total = X_totals2,
                               method = "mle",
                               K = K,
                               controls = list(tol = 1e-6,
                                               maxiter = 500,
                                               merge_move = TRUE))

plot(out_Proc_MLE)

object <- out_Proc_MAP
plot_avgLoads(Betas = out_Proc_MAP$Betas, Xbar = out_Proc_MAP$Xbar, AttrSummary = out_Proc_MAP$AttrSummary)

Betas <- object$Betas
Xbar <- object$Xbar
AttrSummary <- object$AttrSummary


out_Proc_MAP$AttrSummary
Sigs <- out_Proc_MAP$Signatures
Loads <- out_Proc_MAP$Betas

Probs <- out_Proc_MAP$Sig_Cov_prob$
dimnames(Probs)[[2]] <- colnames(Sigs)
dimnames(Probs)[[3]] <- dimnames(Loads)[[3]]


MaxProbIds <- as.data.frame(out_Proc_MAP$Sig_Cov_prob$MaxProbIds) %>%
  group_by(V1, V2) %>%
  summarise(Fraction = n()/nrow(Mutations)) %>%
  as.data.frame()
MaxProbIds$V1 <- colnames(Sigs)[MaxProbIds$V1]
MaxProbIds$V2 <- dimnames(Loads)[[3]][MaxProbIds$V2]
colnames(MaxProbIds) <- c("Signature", "Covariate", "Fraction")


Xest <- estimateTotalCounts(R = out_Proc_MAP$Signatures,
                    Betas = out_Proc_MAP$Betas, X_total = out_Proc_MAP$Xbar[1, ,])

plot(Xest, MutMatrix)

round(out_Proc_MAP$Mu, 5)

Sigs <- out_Proc_MAP$Signatures[, rowMeans(out_Proc_MAP$Mu) > 0.001]
orders <- match_to_RefSigs(Sigs, ref = CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37) %>%
  arrange(best) %>% pull(sigs)
Sigs <- Sigs[, orders]
colnames(Sigs) <- paste0("Sig", LETTERS[1:ncol(Sigs)])

p1 <- plot_96SBSsig(Sigs, returnList = FALSE)


AvgBeta <- round(apply(out_Proc_MAP$Betas * out_Proc_MAP$Xbar, c(1,3), mean)[orders, ], 5)
rownames(AvgBeta) <- colnames(Sigs)
dfMu <- reshape2::melt(as.matrix(AvgBeta))
colnames(dfMu) <- c("Signature", "Feature", "Value")


out_Proc_MAP$Sig_Cov_prob[1, ,]
dimnames(out_Proc_MAP$Sig_Cov_prob) <- list(NULL,
                                            rownames(out_Proc_MAP$Betas),
                                            dimnames(out_Proc_MAP$Betas)[[3]])

out_Proc_MAP$Sig_Cov_prob <- out_Proc_MAP$Sig_Cov_prob[,orders , ]
dimnames(out_Proc_MAP$Sig_Cov_prob)[[2]] <- colnames(Sigs)

dim(out_Proc_MAP$Sig_Cov_prob)

mutProbs <- sapply(1:nrow(out_Proc_MAP$Sig_Cov_prob), function(i) {
  id <- which(out_Proc_MAP$Sig_Cov_prob[i, ,] == max(out_Proc_MAP$Sig_Cov_prob[i, ,]),
              arr.ind = TRUE)
  t(c(colnames(Sigs)[id[1,1]], dimnames(out_Proc_MAP$Sig_Cov_prob)[[3]][id[1, 2]]))
})

mutProbs <- t(mutProbs)

colnames(mutProbs) <- c("Signature", "Feature")
dfProbs <- data.frame(mutProbs) %>%
  group_by(Signature, Feature) %>%
  summarise(Prob = n()/nrow(mutProbs))


names(out_Proc_MAP$Sig_Cov_prob[2, id[1,1], id[1, 2]])

mutProbs <- apply(out_Proc_MAP$Sig_Cov_prob, 1, function(x) {
  id <- which(x == max(x), arr.ind = TRUE)
  c(colnames(Sigs)[id[1,1]], dimnames(out_Proc_MAP$Sig_Cov_prob)[[3]][id[1, 2]])
})


plot(mutProbs,AvgBeta)
summary(lm(c(mutProbs) ~ c(AvgBeta)))

Mu <- round(out_Proc_MAP$Mu[orders, ], 5)
rownames(Mu) <- colnames(Sigs)
dfMu <- reshape2::melt(as.matrix(Mu))
colnames(dfMu) <- c("Signature", "Feature", "Value")


corrplot::corrplot(cor(X))

dfMu2 <- dfMu %>%
  mutate(
    Feature = factor(Feature, levels=unique(Feature)), # ensure order
    AltCol = ifelse((as.numeric(Feature) %% 2) == 0, "gray93", "gray97") # or "gray90" etc
  )

p2 <- ggplot(dfMu2 %>%
               left_join(dfProbs, by = c("Signature", "Feature")) %>%
               mutate(Value = dplyr::case_when(Value < 0.01 ~ NA_real_,
                                         TRUE ~ Value)),
       aes(x = Feature, y = Signature)) +
  # background tiles: use AltCol (discrete!)
  geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
  scale_fill_manual(
    values = c("gray93" = "gray93", "gray97" = "gray97"),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  # data points: use Value (continuous)
  geom_point(aes(fill = Value, size = Prob),
             shape = 21, color = "grey30", stroke = 0.2,
             na.rm = TRUE) +
  scale_fill_gradientn(name = "Average\nloading",
                       colours =c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                                  "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
  scale_y_discrete(limits = rev(levels(dfMu$Signature))) +
  scale_size_continuous(name = "Fraction \nof total\nmutations")+
  scale_x_discrete(position = "top") +
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.ticks.x = element_blank()
  )

p1 + p2

tmp <- dfMu2 %>%
  left_join(dfProbs, by = c("Signature", "Feature")) %>%
  mutate(Value = dplyr::case_when(Value < 0.01 ~ NA_real_,
                                  TRUE ~ Value))

library(patchwork)
plist[[1]] / plist[[2]]
plist[[1]] / plist[[2]] / plist[[3]] +
  plot_layout(heights = c(1, 1, 1), )

round(out_Proc_MAP$Mu, 5)
Mu








set.seed(10)
object <- SigPoisProcess(Mutations = Mutations,
                                X = X,
                                X_total = X_totals2,
                                method = "map",
                                prior_params = list(a = 1.1, alpha = 1.1),
                                K = 10,
                                controls = list(tol = 1e-6,
                                                maxiter = 500,
                                                merge_move = TRUE))






# EDA

df_tumor <- as.data.frame(TumorData$gr_tumor)

MM <- as.matrix(df_tumor %>% dplyr::select(CTCF:GCsq))

bestMark <- apply(MM, 1, function(x) names(which.max(x)))
plot(table(bestMark))

df_tumor %>%
  group_by(sample, channel) %>%
  mutate(
    mut_type = str_extract(channel, "(?<=\\[)[^\\]]+"),
    across(
    .cols = CTCF:GC,
    .fns = mean,
    .names = "{.col}"
  )) %>%
  distinct() %>%
  dplyr::select(sample, channel, mut_type, CTCF, DNAse.seq, H3K9me3) %>%
  gather(key = "mark", value = "n", -sample, -channel, -mut_type) %>%
  #mutate(n = case_when(sqrt(n) > 6.5 ~ 6.5^2,
  #                     TRUE ~ n)) %>%
  ggplot() +
  geom_tile(aes(x = channel, y = sample, fill = n)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
        axis.text.y = element_text(size = 3),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0.1, "lines"))+
  facet_grid(mark~mut_type, scales = "free") +
  scale_fill_viridis_c(name = "signal", option="plasma", direction = -1)














