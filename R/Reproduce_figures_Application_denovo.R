################################################################################
# This file plots Figure 3 and 4 in the main paper
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











