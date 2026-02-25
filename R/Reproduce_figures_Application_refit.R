################################################################################
# This file plots Figure 5 in the main paper
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

plot_betas_capped <- function(Betas_sol, lowCI, highCI) {
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
      Signature = factor(Signature, levels = colnames(Betas_sol)),
      Beta = ifelse(lowCI < 0 & highCI > 0, NA, Beta)
    ) %>%
    # Add capped value for colors only
    mutate(
      Beta_capped = case_when(
        is.na(Beta) ~ NA_real_,
        Beta < -1 ~ -1,
        Beta >  1 ~  1,
        TRUE      ~ Beta
      )
    )

  oldopt <- options()
  options(ggplot2.continuous.colour="ignore", ggplot2.continuous.fill="ignore")

  plot_out <- ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta_capped)) +
    geom_tile(color = "black", width = 1, height = 1, na.rm = TRUE) +
    scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"),
      values = scales::rescale(c(-1, 0, 1)),
      limits = c(-1, 1),
      na.value = "gray90"
    ) +
    geom_text(aes(label = ifelse(is.na(Beta), "", round(Beta, 2))),
              size = 3, na.rm = TRUE) +
    scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 10, colour = "black"),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks = element_blank(),
      panel.grid = element_blank())

  options(oldopt)
  return(plot_out)
}

plot_vector_facets_x <- function(df_Assign, levs) {
  require(ggplot2)
  require(dplyr)

  df_Assign <- df_Assign %>%
    dplyr::mutate(
      best_sig = factor(best_sig, levels = levs),
      compressed = ifelse(m == 0 | mu < 0.05, "compressed", "not compressed")
    )

  gg <- ggplot(df_Assign, aes(x = 1, y = 1, size = mu, fill = m, shape = compressed)) +
    geom_point(color = "black", stroke = 0.7) +
    facet_wrap(~ best_sig, ncol = 1, strip.position = "left") +
    scale_size(name = expression(mu[k]), range = c(3, 12)) +
    scale_fill_gradientn(
      name = "N. mutations",
      colours = c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                  "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")
    ) +
    scale_shape_manual(
      name = "Compressed",
      values = c("compressed" = 4, "not compressed" = 21)
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      strip.background = element_rect(fill = 'white', color = NA),
      legend.position = "left",
      panel.spacing = grid::unit(0.03, "lines")
    )

  return(gg)
}

#--- Useful functions
Compute_mutation_Probs <- function(gr_Mutations, R, Theta, Betas){
  X <- as.matrix(mcols(gr_Mutations)[, -c(1:3)])
  Probs <- R[gr_Mutations$channel, ] * t(Theta[, gr_Mutations$sample]) * exp(X %*% Betas)
  Probs <- t(apply(Probs, 1, function(x) x/sum(x)))
  colnames(Probs) <- colnames(Betas)
  return(Probs)
}

get_signal_regions <- function(chr, region_start, region_end, gr_SignalTrack) {
  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, gr_SignalTrack)
  gr_subset <- gr_SignalTrack[subjectHits(overlaps)]
  return(gr_subset)
}

find_mutations_in_region <- function(data,
                                     i, j,
                                     chr,
                                     region_start,
                                     region_end,
                                     R,
                                     Theta,
                                     Betas){
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  gr_mut_filter <- data$gr_Mutations[data$gr_Mutations$sample == j & data$gr_Mutations$channel == i]
  overlaps <- findOverlaps(gr_region, gr_mut_filter)
  positions <- start(gr_mut_filter[subjectHits(overlaps)])

  # Now, find the value of the signal at the mutations
  X <- as.matrix(mcols(get_signal_regions(chr, positions, positions + 1, data$gr_SignalTrack)))[, -1]
  intensity <- t(R[i, ] * Theta[, j])[rep(1, length(positions)), ] * exp(X %*% Betas)
  max_int <- apply(intensity, 1, which.max)
  df_pos <- data.frame("positions" = positions, "Sig" = colnames(Betas)[max_int])
  return(df_pos)
}

calculate_ChannelProbs_region <- function(chr,
                                          region_start, region_end,
                                          R,
                                          Theta,
                                          Betas,
                                          j = 1,
                                          gr_SignalTrack,
                                          CopyTrack) {


  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, gr_SignalTrack)
  gr_subset <- gr_SignalTrack[subjectHits(overlaps)]
  # Estimate the loadings in each subset
  SignalTrack_subset <- as.matrix(mcols(gr_subset)[, -1])
  CopyTrack_subset <- CopyTrack[subjectHits(overlaps), j]
  ExpSolBetas <- exp(SignalTrack_subset %*% Betas)
  # Big Storage
  ChannelProbs <- array(NA, dim = c(nrow(ExpSolBetas), ncol(ExpSolBetas), 96),
                        dimnames = list(NULL, colnames(ExpSolBetas), rownames(R)))

  ExpSolBetas_copy <- apply(ExpSolBetas, 2, function(x) x * CopyTrack_subset)
  for(i in 1:96){
    SigTheta <- t(R[i, ] * Theta[, j])
    bigProd <- SigTheta[rep(1, nrow(ExpSolBetas_copy)), ] * ExpSolBetas_copy
    sig_probs <-  bigProd#/rowSums(bigProd)
    ChannelProbs[, , i] <- sig_probs
  }
  ChannelProbs
}

################################################################################
# Step 1 - Load the data and the output
################################################################################

# We source a separate script to load the data
source("~/SigPoisProcess/R/Load_ICGC_BreastAdenoCA_data.R")

# Load the output
resultsMCMC_refit <- readRDS("~/SigPoisProcess/output/Application/Breast_refit/resultsMCMC_refit.rds.gzip")

#-- Signatures
Sigs_to_use <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8", "SBS20",
                 "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")
CosmicSigs <- SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37[, Sigs_to_use]

#----------------------------------- Figure 5 panel a
# Regression coefficients
p_Beats_sigs_fixed <- plot_betas_capped(resultsMCMC_refit$Betas$mean,
                                        resultsMCMC_refit$Betas$lowCI,
                                        resultsMCMC_refit$Betas$highCI) + theme(legend.position = "none")

# Relevance weights
bigProd <- resultsMCMC_refit$Signatures[data$gr_Mutations$channel, ] * t(resultsMCMC_refit$Baselines$mean[, data$gr_Mutations$sample]) *
  exp(as.matrix(mcols(data$gr_Mutations)[, -c(1:3)]) %*% resultsMCMC_refit$Betas$mean)
mm <- apply(bigProd, 1, which.max)
df_AssignFix <- data.frame(data$gr_Mutations) %>%
  mutate(best_sig = colnames(CosmicSigs)[mm],
         #best_sig = paste0("SigN", sprintf("%02d", best_sig)),
         best_sig = factor(best_sig, sort(colnames(CosmicSigs)))) %>%
  group_by(sample, best_sig) %>%
  summarize(m = n(), .groups = "drop") %>%
  tidyr::complete(sample, best_sig, fill = list(m = 0)) %>%
  ungroup() %>%
  group_by(best_sig) %>%
  summarize(m = sum(m)) %>%
  left_join(data.frame(best_sig = names(resultsMCMC_refit$Mu$mean),
                       mu = resultsMCMC_refit$Mu$mean), by = "best_sig")
p_mu_sigs_fixed <- plot_vector_facets_x(df_AssignFix,
                                        levs = colnames(resultsMCMC_refit$Betas$mean))

#----------------------------------- Figure 5 panel b

#--- Run the baseline CompressiveNMF
mutMatrix <- getTotalMutations(data$gr_Mutations)
BaseNMF <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                              K = 0,
                                              S = 1e6 * CosmicSigs + 1,
                                              alpha = 1.01, a = 1.01)
dimnames(BaseNMF$mapOutput$R) <- dimnames(CosmicSigs)
colnames(BaseNMF$mapOutput$Theta) <- colnames(mutMatrix)
rownames(BaseNMF$mapOutput$Theta) <- colnames(CosmicSigs)


# Probabilities for each mutation in PPF
MutProbs <- Compute_mutation_Probs(data$gr_Mutations,
                                   R = resultsMCMC_refit$Signatures,
                                   Theta = resultsMCMC_refit$Baselines$mean,
                                   Betas = resultsMCMC_refit$Betas$mean)

# Probabilities for each mutation in NMF
Betas_zero <- matrix(0,
                     nrow = nrow(resultsMCMC_refit$Betas$mean),
                     ncol = ncol(resultsMCMC_refit$Betas$mean),
                     dimnames = dimnames(resultsMCMC_refit$Betas$mean))
MutProbs_fixed <- Compute_mutation_Probs(data$gr_Mutations,
                                         R = BaseNMF$mapOutput$R,
                                         Theta = BaseNMF$mapOutput$Theta,
                                         Betas = Betas_zero)

# Signature with the highest probability in both cases
best_sigs_ppf <- Sigs_to_use[apply(MutProbs, 1, which.max)]
best_sigs_ppf <- factor(best_sigs_ppf, levels = Sigs_to_use)

best_sigs_nmf <- Sigs_to_use[apply(MutProbs_fixed, 1, which.max)]
best_sigs_nmf <- factor(best_sigs_nmf, levels = Sigs_to_use)

df_probs <- data.frame("Best_ppf" = best_sigs_ppf,
                       "Prob_ppf" = apply(MutProbs, 1, max),
                       "Best_nmf" = best_sigs_nmf,
                       "Prob_nmf" = apply(MutProbs_fixed, 1, max))

# Make the plot
p_confusion <- df_probs %>%
  group_by(Best_ppf, Best_nmf) %>%
  summarise(
    count = n(),
    mean_prob_ppf = mean(Prob_ppf),
    .groups = "drop"
  ) %>%
  complete(Best_ppf, Best_nmf) %>%
  mutate(AltCol = ifelse((as.numeric(Best_ppf) %% 2) == 0, "gray93", "gray97")) %>%
  ggplot(aes(x = Best_ppf, y = Best_nmf)) +
  geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
  scale_fill_manual(values = c("gray93", "gray97"), guide = "none") +
  geom_point(aes(size = count, color = mean_prob_ppf), na.rm = TRUE) +
  scale_size_binned(
    range = c(0.5, 10),
    breaks = c(10, 100, 1000, 5000, 10000, 40000, 100000, 140000),
    name = "Count"
  ) +
  ggplot2::scale_color_gradientn(
    name = "Average\nprobability",
    colours = rev(c("#fde725", "#b5de2b", "#6ece58", "#35b779", "#1f9e89",
                    "#26828e", "#31688e", "#3e4989", "#482878", "#440154"))) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(df_probs$Best_nmf))) +
  labs(x = "best_sigs_ppf", y = "best_sigs_nmf", size = "Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    legend.position = "right",        # <--- move legend top
    legend.box = "horizontal"       # <--- stack legends horizontally
  )


#----- Save the outputs
p_mu_sigs_fixed
ggsave("~/SigPoisProcess/figures/Figure5_a_mu.pdf", width = 8.42, height = 5.56)

(p_Beats_sigs_fixed + theme(aspect.ratio = 1)) + p_confusion
ggsave("~/SigPoisProcess/figures/Figure5_a_b.pdf",
       width = 11.68, height = 5.05)


#----------------------------------- Figure 5 panel c
i <- "A[C>A]C"
j <- "DO220823"
chr <- "chr1"
region_start <- 150e6
region_end   <- 160e6

ChannelProbs <- calculate_ChannelProbs_region(chr,
                                              region_start,
                                              region_end,
                                              resultsMCMC_refit$Signatures,
                                              resultsMCMC_refit$Baselines$mean,
                                              resultsMCMC_refit$Betas$mean,
                                              j = j,
                                              data$gr_SignalTrack,
                                              data$CopyTrack)

pr_const <- 2000 *BaseNMF$mapOutput$R[i, ] * BaseNMF$mapOutput$Theta[, j]/(sum(data$gr_SignalTrack$bin_weight))
df_const <- data.frame(t(pr_const)[rep(1,nrow(ChannelProbs)), ]) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1)) * 2000 + region_start,
         method = "Standard NMF")

df_mod <- data.frame(ChannelProbs[, , i])
colnames(df_mod) <- colnames(BaseNMF$mapOutput$R)

df_pos <- find_mutations_in_region(data, i, j, chr,
                                   region_start, region_end,
                                   resultsMCMC_refit$Signatures,
                                   resultsMCMC_refit$Baselines$mean,
                                   resultsMCMC_refit$Betas$mean)

#------- Top panel
window <- 10
options(scipen=99)
p_intensity <- df_mod %>%
  mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
                                                   align = "center"))) %>%
  mutate(region = (0:(nrow(df_mod) - 1))*2000 + region_start,
         method = "PoissonProcess") %>%
  bind_rows(df_const) %>%
  gather(key = "Sig", value = "Prob", -region, -method) %>%
  filter(Sig %in% c("SBS3", "SBS5", "SBS8")) %>%
  ggplot() +
  geom_line(aes(x = region, y = Prob, color = Sig, linetype = method, alpha = method),
            linewidth = 0.6) +
  scale_color_manual(name = "Signature",
                     values = c("red", "antiquewhite3", "blue")) +
  scale_alpha_manual(values = c(1, 1))+
  ylab("Expected mutations in 2kb region") +
  facet_wrap(~paste0("Mutation type ", i, " - Patient ", j)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  geom_rug(data = df_pos,
           aes(x = positions, color = Sig),
           sides = "b", # "b" for bottom
           inherit.aes = FALSE,
           length = grid::unit(0.035, "npc"),
           linewidth = 1.2)+
  theme(axis.title.x = element_blank())
p_intensity

#------- Covariates
subset_tracks <- get_signal_regions(chr, region_start, region_end, data$gr_SignalTrack)
subset_tracks$bin_weight <- NULL
df_tracks <- as.data.frame(mcols(subset_tracks)) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start) %>%
  gather(key = "covariate", value = "score", -region)

p_track <- ggplot(df_tracks) +
  geom_raster(aes(x = region + 1000, y = covariate, fill = score)) +
  #facet_wrap(.~ pos2, nrow = 1) +
  scale_fill_gradient2(name = "Signal track\n(standardized)",
                       mid = "white", low = "#2166AC", high = "#B2182B", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish)+
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  #theme(axis.title = element_blank())
  ylab("Genomic covariate")+
  xlab(paste0("Genomic regions in ",  chr, " (2kb)"))
p_panels <- p_intensity/p_track + plot_layout(heights = c(1.6, 1))
p_panels
ggsave(plot = p_panels, filename = "~/SigPoisProcess/figures/Figure5_c_.png",
       width = 11.68, height = 5.91)
