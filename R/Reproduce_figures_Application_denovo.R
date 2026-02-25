################################################################################
# This file makes the plots for the denovo application
################################################################################
#---- Load packages
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(foreach)
library(sigminer)
library(patchwork)
library(doParallel)
library(GenomicRanges)
library(AnnotationDbi)

# Load the package
library(SigPoisProcess)

######################################################################
# Plotting functions

get_PosteriorEffectiveSize <- function(chain, burnin = 1500){
  if(is.null(dim(chain))){
    coda::effectiveSize(chain[-c(1:burnin)])
  } else if(length(dim(chain)) == 2) {
    apply(chain[-c(1:burnin), ], 2, function(x) coda::effectiveSize(x))
  } else if (length(dim(chain)) == 3){
    apply(chain[-c(1:burnin), ,], c(2,3), function(x) coda::effectiveSize(x))
  }
}


plot_sigs <- function(signatures,
                      lowCI = NULL,
                      highCI = NULL,
                      palette = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) {

  signatures <- as.matrix(signatures)
  names_sig <- rownames(signatures)
  df_plot <- data.frame(signatures) %>%
    dplyr::mutate(
      Channel = names_sig,
      Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                      function(x) paste0(x[c(1,3,7)], collapse = "")),
      Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                       function(x) paste0(x[c(3,4,5)], collapse = "")),
      Mutation = as.factor(Mutation)
    ) %>%
    tidyr::gather(key = "Sig", value = "Prob", -Channel, -Triplet, -Mutation)

  if(!is.null(lowCI) & !is.null(highCI)){
    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         tidyr::gather(key = "Sig", value = "lowCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         tidyr::gather(key = "Sig", value = "highCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig"))
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Triplet, y = Prob, fill = Mutation))+
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::facet_grid(Sig~Mutation, scales = "free", switch = "y") +
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_manual(values = palette)+
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, color = "gray35",
                                          vjust = .5, size = 6.5, margin = ggplot2::margin(t = -4)),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(0, "lines"),
      panel.spacing.y = ggplot2::unit(0, "lines"),
      # The next line makes y facet labels horizontal
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1)
    )

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey65")
  }

  return(p)

}

plot_vector_facets <- function(df_Assign) {
  require(ggplot2)

  df_Assign <- df_Assign %>%
    dplyr::mutate(
      best_sig = factor(best_sig, levels = best_sig)
    )

  ggplot(df_Assign, aes(x = 1, y = 1, size = mu, fill = m)) +
    geom_point(shape = 21, color = "black", stroke = 0.2) +
    facet_wrap(~ best_sig, ncol = 1, strip.position = "left") +
    scale_size(name = expression(mu[k]), range = c(3, 12)) +
    ggplot2::scale_fill_gradientn(
      name = "N. mutations",
      colours = c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                  "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      strip.background = element_rect(fill = 'white', color = NA),  # WHITE BACKGROUND FOR STRIP
      legend.position = "right",
      panel.spacing = grid::unit(0.03, "lines")
    )
}

plot_betas_coef <- function(Betas_sol, lowCI, highCI) {

  df_betas <- as.data.frame(Betas_sol) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "Beta", -Covariate)

  df_low <- as.data.frame(lowCI) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "lowCI", -Covariate)

  df_high <- as.data.frame(highCI) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "highCI", -Covariate)

  df_betas <- df_betas %>%
    left_join(df_low, by = c("Covariate", "Signature")) %>%
    left_join(df_high, by = c("Covariate", "Signature")) %>%
    mutate(
      Covariate = factor(Covariate, levels = unique(Covariate)),
      Signature = as.factor(Signature),
      Beta = ifelse(lowCI < 0 & highCI > 0, NA, Beta)
    )

  max_abs_beta <- max(abs(df_betas$Beta), na.rm = TRUE)

  oldopt <- options()
  options(ggplot2.continuous.colour="ignore", ggplot2.continuous.fill="ignore")

  plot_out <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta)) +
    geom_tile(color = "black", width = 1, height = 1, na.rm = TRUE) +
    scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"),
      limits = c(-max_abs_beta, max_abs_beta),
      na.value = "gray90"
    ) +
    geom_text(aes(label = ifelse(is.na(Beta), "", round(Beta, 2))), size = 3, na.rm = TRUE) +
    scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 10,
                                 colour = "black"),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks = element_blank(), panel.grid = element_blank())

  options(oldopt)

  return(plot_out)
}

plot_baseline_weights <- function(Wmat, groups = NULL, clust = NULL, nclust = 10) {
  Wmat <- t(Wmat)
  if(is.null(clust)){
    clust <- stats::cutree(stats::hclust(dist(Wmat)), k = nclust)
  } else {
    if(length(clust) != nrow(Wmat)){
      stop("Length of clust must be equal to ncol(Wmat)")
    }
  }
  Wmat <- Wmat[order(clust), ]
  pW <- Wmat %>%
    as.data.frame() %>%
    dplyr::mutate(patient = rownames(Wmat),
                  group =1,
                  clust = clust[order(clust)]) %>%
    tidyr::gather(key = "Signature", value = "weight", -patient, -group, -clust) %>%
    dplyr::mutate(Signature = as.factor(Signature),
                  patient = as.factor(patient))%>%
    dplyr::mutate(patient = ordered(patient, unique(patient))) %>%
    ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::geom_bar(ggplot2::aes(x=patient, y = weight, fill = Signature), stat = "identity") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(axis.text.x =  element_blank(),
                   axis.title.x = element_blank(),
                   legend.position = "right")+
    ggplot2::ylab("Baseline activities") +
    facet_wrap(~clust, nrow = 1, scales = "free", space = "free_x")#+
    #ggforce::facet_row(~clust, scales = "free")

  if(!is.null(groups)){
    pW <- pW + ggforce::facet_row(~group, scales = "free", space = "free_x")
  }

  return(pW)
}

################################################################################
# Step 1 - Load the data and the output
################################################################################

# We source a separate script to load the data
source("~/SigPoisProcess/R/Load_ICGC_BreastAdenoCA_data.R")

# Load the output
resultsMCMC_denovo <- readRDS("~/SigPoisProcess/output/Application/Breast_DeNovo/resultsMCMC_denovo.rds.gzip")
postBurnChain_MCMC_denovo <- readRDS("~/SigPoisProcess/output/Application/Breast_DeNovo/postBurnin_MCMC_denovo.rds.gzip")

################################################################################
# Step 2 - Make plots
################################################################################

#-- Extract parameters

#-------------- Mu
mu <- resultsMCMC_denovo$Mu$mean
mu_filter <- mu[mu > 0.01]
sigs_order <- order(mu_filter, decreasing = TRUE)
mu_order <- mu_filter[sigs_order]
names(mu_order) <- paste0("Sig", 1:length(sigs_order))

#-------------- Signatures
SigsMean <- resultsMCMC_denovo$Signatures$mean[, mu > 0.01][, sigs_order]
SigsLowCI <- resultsMCMC_denovo$Signatures$lowCI[, mu > 0.01][, sigs_order]
SigsHighCI <- resultsMCMC_denovo$Signatures$highCI[, mu > 0.01][, sigs_order]
colnames(SigsMean) <- colnames(SigsLowCI) <- colnames(SigsHighCI) <- paste0("Sig", 1:length(sigs_order))

#-------------- Baselines
ThetaMean <- resultsMCMC_denovo$Baselines$mean[mu > 0.01, ][sigs_order, ]
ThetaLowCI <- resultsMCMC_denovo$Baselines$lowCI[mu > 0.01, ][sigs_order, ]
ThetaHighCI <- resultsMCMC_denovo$Baselines$highCI[mu > 0.01, ][sigs_order, ]
rownames(ThetaMean) <- rownames(ThetaLowCI) <- rownames(ThetaHighCI) <- paste0("Sig", 1:length(sigs_order))

#-------------- Beta coefficients
BetasMean <- resultsMCMC_denovo$Betas$mean[, mu > 0.01][, sigs_order]
BetasLowCI <- resultsMCMC_denovo$Betas$lowCI[, mu > 0.01][, sigs_order]
BetasHighCI <- resultsMCMC_denovo$Betas$highCI[, mu > 0.01][, sigs_order]
colnames(BetasMean) <- colnames(BetasLowCI) <- colnames(BetasHighCI) <- paste0("Sig", 1:length(sigs_order))

#-------------- Mutation assignments
bigProd <- SigsMean[data$gr_Mutations$channel, ] * t(ThetaMean[, data$gr_Mutations$sample]) *
  exp(as.matrix(mcols(data$gr_Mutations)[, -c(1:3)]) %*% BetasMean)
mm <- apply(bigProd, 1, which.max)
df_Assign <- data.frame(data$gr_Mutations) %>%
  mutate(best_sig = mm,
         best_sig = paste0("Sig", best_sig)) %>%
  group_by(sample, best_sig) %>%
  summarize(m = n()) %>%
  ungroup() %>%
  group_by(best_sig) %>%
  summarize(m = sum(m)) %>%
  left_join(data.frame(best_sig = names(mu_order), mu = mu_order), by = "best_sig")


########################################### Figure 4

#---------------------------------------------------- Figure 4 panel a, b and d

plot_sigs(SigsMean, SigsLowCI, SigsHighCI) +
  plot_vector_facets(df_Assign) +
  plot_spacer() +
  plot_betas_coef(BetasMean, BetasLowCI, BetasHighCI)  +
  plot_layout(widths = c(1.4, 0.2, 0.05, 1.5))

ggsave("~/SigPoisProcess/figures/Figure4_a_b_c_denovoPars.pdf", width = 12.58, height=5.82)

#---------------------------------------------------- Figure 4 panel d

# Cluster the baseline weights
ThetaCopy <- ThetaMean * t(colSums(data$CopyTrack))[rep(1, 9), ]

Thetas <- apply(ThetaMean, 2, function(x) x/sum(x))
Thetas <- Thetas[, names(sort(colSums(ThetaCopy)))]
Wmat <- ThetaCopy[, names(sort(colSums(ThetaCopy)))]

factoextra::fviz_nbclust(t(Thetas), cluster::pam, method = "silhouette")
clust <- cluster::pam(t(Thetas), k = 3, stand = FALSE)$clustering

#---- Make the plot and save it
col_values <- c("darkblue", "#4959C7", "skyblue1", "#F2AB67", "#F6C866", "#ED7470",
                "brown", "#999999", "black")

plot_baseline_weights(Wmat, clust = clust) +
  scale_fill_manual(values = col_values)
ggsave(filename = "~/SigPoisProcess/figures/Figure4_d_denovoBaselines.pdf", width = 12.34, height=2.82)


########################################### Figure 3

#----------- Panel a : full prediction at Mb scale
LambdaPred <- SigPoisProcess:::Reconstruct_Lambda(SignalTrack = data$SignalTrack,
                                           CopyTrack = data$CopyTrack,
                                           R = SigsMean,
                                           Theta = ThetaMean,
                                           Betas = BetasMean, verbose = TRUE)
LambdaPred <- rowSums(LambdaPred)
LambdaTrack <- data$gr_SignalTrack
LambdaTrack$LambdaPred <- LambdaPred

# Aggregate SignalTrack
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
over <- findOverlaps(LambdaTrack, hg19_Mb)
df_track <- as.data.frame(LambdaTrack) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred))

# Get chromosome for each region (tile)
df_chrom <- as.data.frame(hg19_Mb) %>%
  mutate(region = 1:n())

# Get Number of points in each bin, and join with chromosomes
df_tumor_all <- as.data.frame(data$gr_Mutations) %>%
  mutate(region = subjectHits(findOverlaps(data$gr_Mutations, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  left_join(df_track, by = "region") %>%
  dplyr::left_join(df_chrom %>%
                     dplyr::select(region, seqnames), by = "region") %>%
  dplyr::mutate(
    chr_num = as.integer(gsub("chr", "", seqnames)),
    chrom_color_group = chr_num %% 2,
    chrom_color_group = case_when(is.na(chrom_color_group) ~ 1,
                                  TRUE ~ chrom_color_group),
    point_color_group = as.factor(chrom_color_group),
    line_color_group  = as.factor(chrom_color_group)) %>%
  dplyr::mutate(
    rolling_lambda = zoo::rollmean(Lambda, k = 1, fill = NA, align = "center"))

# Make the plot
p_Track <- ggplot() +
  geom_rect(data = as.data.frame(hg19_Mb) %>%
              mutate(region = 1:length(hg19_Mb)) %>%
              group_by(seqnames) %>%
              summarise(xmin = min(region), xmax = max(region) + 1) %>%
              arrange(seqnames) %>%
              mutate(fill_color = factor(seq_along(seqnames) %% 2)),
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = fill_color),
            alpha = 0.1) +
  scale_fill_manual(values = c("#4682B4", "antiquewhite"))+
  geom_point(data = df_tumor_all, size = 0.66,
             aes(x = region, y = n, color = point_color_group)) +
  scale_color_manual(values = c("#000D8B" ,"#93BAF1"), guide = "none") +
  ggnewscale::new_scale_color() +
  geom_line(data = df_tumor_all, linewidth = 0.7, alpha = 0.7,
            aes(x = region, y = Lambda, color = line_color_group, group = seqnames)) +
  scale_color_manual(values = c("#CD2626", "tomato"), guide = "none") +
  theme_bw() +
  ylab("Number of mutations") +
  xlab("Genomic region (Mb)") +
  scale_x_continuous(expand = c(0.01, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  theme(legend.position = "none")

p_reconstr <- ggplot(df_tumor_all) +
  geom_point(aes(x = Lambda, y = n, color = point_color_group),
             alpha =0.25, size = 0.5, shape = 20) +
  scale_color_manual(values = c("#000D8B" ,"#93BAF1"), guide = "none") +
  theme_bw() +
  theme(aspect.ratio =  1) +
  xlim(c(0, 610)) +
  ylim(c(0, 610)) +
  geom_abline(slope = 1, intercept = 0, color = "#CD2626", linewidth = 0.7) +
  ylab("Observed mutations") +
  xlab("Predicted mutations")


#----------- Panel b: model with only copy numbers

outNoCovs <- readRDS(file = "~/SigPoisProcess/output/Application_Breast/out_FreeSig_noCovariates.rds.gzip")

# Plot the number of samples at the MB scale
LambdaPredNoCovs <- SigPoisProcess:::Reconstruct_Lambda(SignalTrack = as.matrix(data$SignalTrack[, 1]),
                                       CopyTrack = data$CopyTrack,
                                       R = outNoCovs$Signatures,
                                       Theta = outNoCovs$Thetas,
                                       Betas = outNoCovs$Betas, verbose = TRUE)
LambdaPredNoCovs <- rowSums(LambdaPredNoCovs)
LambdaTrackNoCovs <- data$gr_SignalTrack
LambdaTrackNoCovs$LambdaPred <- LambdaPredNoCovs

df_trackNoCovs <- as.data.frame(LambdaTrackNoCovs) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred))

# Get Number of points in each bin, and join with chromosomes
df_tumor_allNoCovs <- as.data.frame(data$gr_Mutations) %>%
  mutate(region = subjectHits(findOverlaps(data$gr_Mutations, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  left_join(df_trackNoCovs, by = "region") %>%
  dplyr::left_join(df_chrom %>%
                     dplyr::select(region, seqnames), by = "region") %>%
  dplyr::mutate(
    chr_num = as.integer(gsub("chr", "", seqnames)),
    chrom_color_group = chr_num %% 2,
    chrom_color_group = case_when(is.na(chrom_color_group) ~ 1,
                                  TRUE ~ chrom_color_group),
    point_color_group = as.factor(chrom_color_group),
    line_color_group  = as.factor(chrom_color_group)) %>%
  dplyr::mutate(
    rolling_lambda = zoo::rollmean(Lambda, k = 2, fill = NA, align = "center"))

p_TrackNoCovs <- ggplot() +
  geom_rect(data = as.data.frame(hg19_Mb) %>%
              mutate(region = 1:length(hg19_Mb)) %>%
              group_by(seqnames) %>%
              summarise(xmin = min(region), xmax = max(region) + 1) %>%
              arrange(seqnames) %>%
              mutate(fill_color = factor(seq_along(seqnames) %% 2)),
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = fill_color),
            alpha = 0.1) +
  scale_fill_manual(values = c("#4682B4", "antiquewhite"))+
  geom_point(data = df_tumor_allNoCovs, size =0.66,
             aes(x = region, y = n, color = point_color_group)) +
  scale_color_manual(values = c("#000D8B" ,"#93BAF1"), guide = "none") +
  ggnewscale::new_scale_color() +
  geom_line(data = df_tumor_allNoCovs, linewidth = 0.8,alpha =0.7,
            aes(x = region, y = Lambda, color = line_color_group, group = seqnames)) +
  scale_color_manual(values = c( "#CD2626", "tomato"), guide = "none") +
  theme_bw() +
  ylab("Number of mutations")+
  xlab("Genomic region (Mb)") +
  scale_x_continuous(expand = c(0.01, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  theme(legend.position = "none")

p_reconstrNoCovs <- ggplot(df_tumor_allNoCovs) +
  geom_point(aes(x = Lambda, y = n, color = point_color_group),
             alpha =0.25, size = 0.5, shape = 20) +
  scale_color_manual(values = c("#000D8B" ,"#93BAF1"), guide = "none") +
  theme_bw() +
  theme(aspect.ratio =  1) +
  xlim(c(0, 610)) +
  ylim(c(0, 610)) +
  geom_abline(slope = 1, intercept = 0,color = "#CD2626", linewidth = 0.7) +
  ylab("Observed mutations") +
  xlab("Predicted mutations")

(p_Track + p_reconstr) / (p_TrackNoCovs + p_reconstrNoCovs)
ggsave("../figures/Figure3_recontructed_Mbscale.pdf", width = 9.25, height=4.65)


########################################### Additional results in the supplement

#---- Figure S2 - correlation across covariates
# Correlation an the whole track
corrplot::corrplot(
  cor(data$SignalTrack),
  type = "lower",
  addCoef.col = "gray25",
  number.cex = 0.7,
  diag = FALSE
)

# Correlation an the mutations
corrplot::corrplot(
  cor(as.matrix(mcols(data$gr_Mutations[, -c(1,2,3)]))),
  type = "lower",
  addCoef.col = "gray25",
  number.cex = 0.7,
  diag = FALSE
)

################################################ Figure S3
# Confrontation between CompNMF, SignatureAnalyzer, and SigPoisProcess
mutMatrix <- getTotalMutations(data$gr_Mutations)

if(rerun_competitors) {
  set.seed(10)
  BaselineNMF <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                                    K = 30,
                                                    alpha = 1.01, a = 1.01)

  set.seed(10)
  out_SigAnalyzerL1KL <- sigminer::sig_auto_extract(nmf_matrix = t(mutMatrix), method = "L1KL", cores = 20)
  set.seed(10)
  out_SigAnalyzerL1W.L2H <- sigminer::sig_auto_extract(nmf_matrix = t(mutMatrix), method = "L1W.L2H", cores = 20)

  saveRDS(BaselineNMF, file = "~/SigPoisProcess/output/Application/Breast_DeNovo/BaselineNMF/out_CompNMF.rds")
  saveRDS(out_SigAnalyzerL1KL, file = "~/SigPoisProcess/output/Application/Breast_DeNovo/BaselineNMF/out_SigAnalyzerL1KL.rds")
  saveRDS(out_SigAnalyzerL1W.L2H, file = "~/SigPoisProcess/output/Application/Breast_DeNovo/BaselineNMF/out_SigAnalyzerL1WL2H.rds")

}

BaselineNMF <- readRDS(file = "~/SigPoisProcess/output/Application/Breast_DeNovo/out_CompNMF.rds")
out_SigAnalyzerL1KL <- readRDS(file = "~/SigPoisProcess/output/Application/Breast_DeNovo/out_SigAnalyzerL1KL.rds")
out_SigAnalyzerL1W.L2H <- readRDS(file = "~/SigPoisProcess/output/Application/Breast_DeNovo/out_SigAnalyzerL1WL2H.rds")

plot_sigs(BaselineNMF$Signatures)
match_to_RefSigs(BaselineNMF$Signatures, SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37)
plot_sigs(out_SigAnalyzerL1KL$Signature.norm)
match_to_RefSigs(out_SigAnalyzerL1KL$Signature.norm, SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37)
plot_sigs(out_SigAnalyzerL1W.L2H$Signature.norm)

MatEst <- SigPoisProcess:::Reconstruct_CountMatrix(SignalTrack = data$SignalTrack,
                                                   CopyTrack = data$CopyTrack,
                                                   R = SigsMean, Theta = ThetaMean,
                                                   Betas = BetasMean)
# Reconstructon PPF
plot(MatEst, mutMatrix)
rmse_1 <- sqrt(mean((MatEst - mutMatrix)^2))
# Reconstructon CompNMF
plot(BaselineNMF$Signatures %*% BaselineNMF$Theta, mutMatrix)
rmse_2 <- sqrt(mean((BaselineNMF$Signatures %*% BaselineNMF$Theta - mutMatrix)^2))
# Reconstructon SignatureAnalyzer
plot(out_SigAnalyzerL1KL$Signature.norm %*% out_SigAnalyzerL1KL$Exposure, mutMatrix)
rmse_3 <- sqrt(mean((out_SigAnalyzerL1KL$Signature.norm %*% out_SigAnalyzerL1KL$Exposure - mutMatrix)^2))

round(c(rmse_1, rmse_2, rmse_3), 2)

# Match with baseline
match_w_baseline <- match_MutSign(R_true = SigsMean,
                                  R_hat = BaselineNMF$Signatures)
colnames(match_w_baseline$R_hat) <- colnames(match_w_baseline$R_true)


# Match with SignatureAnalyzer
match_w_sigAnalyzer <- match_MutSign(R_true = SigsMean,
                                     R_hat = out_SigAnalyzerL1KL$Signature.norm)
colnames(match_w_sigAnalyzer$R_hat) <- colnames(match_w_sigAnalyzer$R_true)
plot_sigs(match_w_baseline$R_true) + plot_sigs(match_w_baseline$R_hat) + plot_sigs(match_w_sigAnalyzer$R_hat)
ggsave("~/SigPoisProcess/figures/Breast_suppl_Signatures_comparison.pdf",
       width = 13.42, height = 5.64)


###################################################### Figure S4
# Calculate effective sample sizes
mat_allCopy <- t(colSums(data$CopyTrack))[rep(1, sum(mu > 0.01)), ]
ThetaChain_adj <- postBurnChain_MCMC_denovo$Baselines[, mu > 0.01, ]
array_AllCopy <- array(NA, dim = c(3000, 9, 113))
for(s in 1:3000) ThetaChain_adj[s,,] <- mat_allCopy * ThetaChain_adj[s, , ]

REffects <- get_PosteriorEffectiveSize(postBurnChain_MCMC_denovo$Signatures[, , mu > 0.01], burnin = 0)
ThetaEffects <- get_PosteriorEffectiveSize(ThetaChain_adj, burnin = 0)
BetasEffects <- get_PosteriorEffectiveSize(postBurnChain_MCMC_denovo$Betas[, , mu > 0.01], burnin = 0)
Sigma2Effects <- get_PosteriorEffectiveSize(postBurnChain_MCMC_denovo$sigma2[, mu > 0.01], burnin = 0)
MuEffects <- get_PosteriorEffectiveSize(postBurnChain_MCMC_denovo$Mu[, mu > 0.01], burnin = 0)
lpEffects <- get_PosteriorEffectiveSize(postBurnChain_MCMC_denovo$logPost, burnin = 0)

effects_df <- tibble(
  value = c(REffects, ThetaEffects, BetasEffects, MuEffects),
  quantity = factor(rep(
    c("R", "Theta", "B", "mu"),
    times = c(length(REffects), length(ThetaEffects), length(BetasEffects),
              length(MuEffects))
  ),
  levels = c("R", "Theta", "B", "mu"))
)

p_effective <- ggplot(effects_df, aes(x = "", y = value)) +
  geom_boxplot() +
  facet_wrap(~ quantity, nrow = 1, scales = "free") +
  labs(x = "", y = "Effective Sample Size") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_logPost <- ggplot(data = data.frame(iter = 1:3000,
                                      logPosterior = postBurnChain_MCMC_denovo$logPost))+
  geom_line(aes(x = iter, y = logPosterior)) +
  theme_bw()+
  facet_wrap(~ "logposterior trace", nrow = 1) +
  xlab("post-burnin iteration")

p_ess <- p_effective + p_logPost + plot_layout(nrow = 1, widths =  c(2, 1))
p_ess

###################################################### Additional results

#----  Summary of the MAP outputs
outMAP1 <- readRDS("~/SigPoisProcess/output/Application/Breast_DeNovo/MapSolutions/Replicate01/MAPSolution.rds.gzip")
outMAP2 <- readRDS("~/SigPoisProcess/output/Application/Breast_DeNovo/MapSolutions/Replicate02/MAPSolution.rds.gzip")
outMAP3 <- readRDS("~/SigPoisProcess/output/Application/Breast_DeNovo/MapSolutions/Replicate03/MAPSolution.rds.gzip")

times_all <- rbind(c("iter" = outMAP1$sol$iter, "time" = as.numeric(outMAP1$time, unit = "secs")),
                   c("iter" = outMAP2$sol$iter, "time" = as.numeric(outMAP2$time, unit = "secs")),
                   c("iter" = outMAP3$sol$iter, "time" = as.numeric(outMAP3$time, unit = "secs")))

# Plot the times to run the algoithms
(times_all[, 2]/times_all[, 1]) * 60

# How many signatures were selected in each case?
sum(outMAP1$Mu > 0.001)
sum(outMAP2$Mu > 0.001)
sum(outMAP3$Mu > 0.001)


match12 <- match_MutSign(outMAP1$Signatures, outMAP2$Signatures)
mean(diag(cosine(match12$R_hat, match12$R_true)))

match13 <- match_MutSign(match12$R_true, outMAP3$Signatures)
mean(diag(cosine(match13$R_hat, match13$R_true)))

match23 <- match_MutSign(outMAP2$Signatures, outMAP3$Signatures)
mean(diag(cosine(match23$R_hat, match23$R_true)))

plot_96SBSsig(match12$R_hat) + plot_96SBSsig(match12$R_true) + plot_96SBSsig(match23$R_hat)

#---- Introduction: difference between gr for chromosome 1p and 1q
chr1_arms_gr <- GRanges(
  seqnames = c("chr1", "chr1"),
  ranges = IRanges(
    start = c(1, 124535434),
    end = c(121535433, 249250621)
  ),
  arm = c("1p", "1q")
)

df_chrom1pq <- as.data.frame(data$gr_Mutations)
df_chrom1pq$region <- NA
over1pq <- findOverlaps(data$gr_Mutations, chr1_arms_gr)
df_chrom1pq$region[queryHits(over1pq)] <- subjectHits(over1pq)

df_copy1pq <- as.data.frame(data$gr_CopyTrack)
df_copy1pq$region <- NA
overCopy1pq <- findOverlaps(data$gr_CopyTrack, chr1_arms_gr)
df_copy1pq$region[queryHits(overCopy1pq)] <- subjectHits(overCopy1pq)

df_chrom1pq %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  left_join(df_copy1pq %>%
              dplyr::select(region, DO1076:DO52562) %>%
              drop_na() %>%
              gather(key = "key", value = "value", -region) %>%
              group_by(region) %>%
              summarise(Copies = 2 * mean(value)), by = "region")




