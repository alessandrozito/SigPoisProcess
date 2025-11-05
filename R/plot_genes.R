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
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(org.Hs.eg.db); #BiocManager::install("org.Hs.eg.db")
devtools::document()

create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}


#----- Useful functions
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
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1, size = 10,
                                 colour = "black"),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks = element_blank(), panel.grid = element_blank())

  options(oldopt)

  return(plot_out)
}

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
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1, size = 10, colour = "black"),
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

#------

#----- The tumor data
file_data <- "~/SigPoisProcess/data/ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip"
data <- readRDS(file_data)

#------ PART 1 - Plots for the Fully Unsupervised case

# Load the output from the MCMC
outMCMC <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate03/MCMCSolution.rds.gzip")

# Extract posterior quantities and confidence intervals
burnin <- c(1:7000)
sigma2 <- colMeans(outMCMC$chain$SIGMA2chain[-burnin, ])
mu <- colMeans(outMCMC$chain$MUchain[-burnin, ])

SigsCis <- get_PosteriorCI(outMCMC$chain$SIGSchain[-burnin, , mu > 0.01])
ThetaCis <- get_PosteriorCI(outMCMC$chain$THETAchain[-burnin, mu > 0.01, ])
BetasCis <- get_PosteriorCI(outMCMC$chain$BETASchain[-burnin, , mu > 0.01])


# Plot the model now. First, sort the signatures by their relevance weigth.
mu_filter <- mu[mu > 0.01]
sigs_order <- order(mu_filter, decreasing = TRUE)
mu_order <- mu_filter[sigs_order]
names(mu_order) <- paste0("Sig", 1:length(sigs_order))

#-------------- Assignment probabilities
# Fraction of mutations attributed to the signature
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

p_mu <- plot_vector_facets(df_Assign)

#-------------- Signatures
SigsMean <- SigsCis$mean[, sigs_order]
SigsLowCI <- SigsCis$lowCI[, sigs_order]
SigsHighCI <- SigsCis$highCI[, sigs_order]
colnames(SigsMean) <- colnames(SigsLowCI) <- colnames(SigsHighCI) <- paste0("Sig", 1:length(sigs_order))

#-------------- Baselines
ThetaMean <- ThetaCis$mean[sigs_order, ]
ThetaLowCI <- ThetaCis$lowCI[sigs_order, ]
ThetaHighCI <- ThetaCis$highCI[sigs_order, ]
rownames(ThetaMean) <- rownames(ThetaLowCI) <- rownames(ThetaHighCI) <- paste0("Sig", 1:length(sigs_order))

#-------------- Beta coefficients
BetasMean <- BetasCis$mean[, sigs_order]
BetasLowCI <- BetasCis$lowCI[, sigs_order]
BetasHighCI <- BetasCis$highCI[, sigs_order]
colnames(BetasMean) <- colnames(BetasLowCI) <- colnames(BetasHighCI) <- paste0("Sig", 1:length(sigs_order))

plot_sigs(SigsMean, SigsLowCI, SigsHighCI) +
  plot_spacer() +                # <--- Separator
  plot_betas_coef(BetasMean, BetasLowCI, BetasHighCI)  +
  p_mu +
  plot_layout(widths = c(1.4, 0.05, 1.5, 0.2))
ggsave("../figures/Breast_Output_Unsup_MCMC2.pdf", width = 12.58, height=5.82)


# Match to reference signatures
match_to_RefSigs(SigsMean, CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)

# Plot the model parameters
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
    ggplot2::ylab("Baseline activities")+
    ggforce::facet_row(~clust, scales = "free", space = "free")

  if(!is.null(groups)){
    pW <- pW + ggforce::facet_row(~group, scales = "free", space = "free")
  }

  return(pW)
}
col_values <- c("darkblue", "#4959C7", "skyblue1", "#F2AB67", "#F6C866", "#ED7470",
                "brown", "#999999", "black")
ThetaCopy <- ThetaMean * t(colSums(data$CopyTrack))[rep(1, 9), ]

Thetas <- apply(ThetaMean, 2, function(x) x/sum(x))
Thetas <- Thetas[, names(sort(colSums(ThetaCopy)))]
Wmat <- ThetaCopy[, names(sort(colSums(ThetaCopy)))]

# groups <- 1*(colnames(Wmat) %in% names(sort(colSums(Wmat)))[1:109])
# groups <- ifelse(groups == 1, "1.Normal", "2.Hypermutated")
factoextra::fviz_nbclust(t(Thetas), cluster::pam, method = "silhouette")
clust <- cluster::pam(t(Thetas), k = 3, stand = FALSE)$clustering

pup <- plot_baseline_weights(Thetas, clust = clust) +
  scale_fill_manual(values = col_values)
pdown <- plot_baseline_weights(Wmat, clust = clust) +
  scale_fill_manual(values = col_values)
pup/pdown

pdown
ggsave(plot = pdown, filename = "../figures/Breast_Output_Unsup_MCMC_Activities_clust.pdf", width = 12.34, height=2.82)

#------- Plot Mutation rate at the megabase scale
LambdaPred <- Reconstruct_Lambda(SignalTrack = data$SignalTrack,
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
p_Track

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

#ED7470

# ggpubr::ggarrange(p_Theta, p_Track, widths = c(1, 1.5))
# p_Theta +
# p_Track +
#   plot_layout(widths = c(0.9, 1.2))
# ggsave("../figures/Breast_Output_Unsup_MCMC.pdf", width = 12.34, height=5.82)

# Finally, do the same plot for the contant mutation rate
outNoCovs <- readRDS(file = "~/SigPoisProcess/output/Application_Breast/out_FreeSig_noCovariates.rds.gzip")

# Plot the number of samples at the MB scale
LambdaPredNoCovs <- Reconstruct_Lambda(SignalTrack = as.matrix(data$SignalTrack[, 1]),
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

#CD2626 tomato #000D8B  #4682B4 #antiquewhite #93BAF1

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

p_Track + p_TrackNoCovs
ggsave("../figures/Breast_Output_Unsup_MCMC_Track.pdf", width = 14.5, height=3.0)

(p_Track + p_reconstr) / (p_TrackNoCovs + p_reconstrNoCovs)
ggsave("../figures/Breast_Output_Unsup_MCMC_Track_points2.pdf", width = 11.32, height=4.15)
ggsave("../figures/Breast_Output_Unsup_MCMC_Track_points2.pdf", width = 9.25, height=4.65)

sqrt(mean((df_tumor_allNoCovs$n - df_tumor_allNoCovs$Lambda)^2))
sqrt(mean((df_tumor_all$n - df_tumor_all$Lambda)^2))

#############################################
# EDA for the introduction
#############################################
df_copy <- read_tsv("~/20170119_final_consensus_copynumber_donor") %>%
  filter(sampleID %in% colnames(data$CopyTrack))

gr_copy <- GRanges(
  seqnames = paste0("chr",df_copy$chr),
  ranges = IRanges(start = df_copy$start, end = df_copy$end),
  strand = "*",
  sample = df_copy$sampleID,
  score = df_copy$value)

tmp <- as.data.frame(gr_copy) %>%
  group_by(sample) %>%
  filter(is.na(score)) %>%
  summarize(w = sum(width))

gr_copy$score[is.na(gr_copy$score)] <- 2

gr_CopyTrack <- data$gr_CopyTrack
gr_CopyTrack$bin_weight <- NULL
mm <- rowMeans(as.matrix(mcols(gr_CopyTrack)) * 2)
mcols(gr_CopyTrack) <- mm

df_copy_agg <- as.data.frame(data$gr_Mutations) %>%
  mutate(region = subjectHits(findOverlaps(data$gr_Mutations, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  left_join(as.data.frame(gr_CopyTrack) %>%
              mutate(region = subjectHits(findOverlaps(gr_CopyTrack, hg19_Mb))) %>%
              group_by(region) %>%
              summarize(Copies = mean(X)), by = "region")

cor(df_copy_agg$n, df_copy_agg$Copies)
cor.test(df_copy_agg$n, df_copy_agg$Copies)

plot(df_copy_agg$n, df_copy_agg$Copies)

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

df_copy1pq <- as.data.frame(gr_CopyTrack)
df_copy1pq$region <- NA
overCopy1pq <- findOverlaps(gr_CopyTrack, chr1_arms_gr)
df_copy1pq$region[queryHits(overCopy1pq)] <- subjectHits(overCopy1pq)

df_chrom1pq %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  left_join(df_copy1pq %>%
              group_by(region) %>%
              summarize(Copies = mean(X)), by = "region")


#------------------------------------------------------------------
# I need to plot
# 1) Correlation among covariates both along the track and at the mutation points

set.seed(10)
subset_points <- sample(nrow(data$SignalTrack), size = 20000)
VarsTemp <- data$SignalTrack[subset_points, ]

# Color function: blue-white-red, correlation [-1, 1]
corr_col_fun <- scales::col_numeric(
  palette = c("blue", "white", "red"),
  domain = c(-0.75, 0.75)
)

library(GGally)
# Custom correlation panel with colored correlation values
# Custom correlation panel with colored correlation values and smaller font
cor_fun_col <- function(data, mapping, ...){
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  corval <- cor(x, y, use = "complete.obs")
  col <- corr_col_fun(corval) # from previous code
  GGally::ggally_cor(data, mapping, ..., size = 5, color = col)
}

# ggpairs plot
p <- GGally::ggpairs(
  VarsTemp[rowSums(VarsTemp > 3) ==0,],
  lower =  list(continuous = GGally::wrap("density", alpha = 0.7)),
  upper = list(continuous = cor_fun_col),
  diag = list(continuous = "densityDiag")
)

print(p)

# I simply will plot the correlation between all variables, and then plot the density
# at the mutation points

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

#------ Density at the mutation points
df_allData <- data.frame(data$SignalTrack) %>%
  mutate(across(where(is.numeric), ~pmin(.x, quantile(.x, 0.99, na.rm = TRUE))))
df_allData$type <- "All SignalTrack"

df_MutData <- as.data.frame(data$gr_Mutations) %>%
  mutate(type =  str_extract(channel, "(?<=\\[)[^\\]]+(?=\\])")) %>%
  dplyr::select(GC:type) %>%
  mutate(across(where(is.numeric), ~pmin(.x, quantile(.x, 0.99, na.rm = TRUE))))

rbind(df_allData, df_MutData) %>%
  gather(key = "Variable", value = "value", -type) %>%
  ggplot() +
  geom_density(aes(x = value, color = type, fill = type)) +
  facet_wrap(.~Variable, scales = "free") +
  scale_fill_manual(values = c("magenta", NA, NA, NA, NA, NA, NA))+
  scale_color_manual(values = c("magenta", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) +
  theme_bw()

ggplot() +
  geom_density(data = df_allData %>%

                 gather(key = "Variable", value = "value", -type),
               aes(x = value), color = "darkblue", linetype = "dotted") +
  geom_density(data = df_MutData %>%
                 gather(key = "Variable", value = "value", -type),
               aes(x = value, color = type)) +
  facet_wrap(.~Variable, scales = "free") +
  scale_color_manual(values = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) +
  theme_bw()


### Extract signatures from the MAP
outMAP1 <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate01/MAPSolution.rds.gzip")
outMAP2 <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate02/MAPSolution.rds.gzip")
outMAP3 <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate03/MAPSolution.rds.gzip")

times_all <- rbind(c("iter" = outMAP1$sol$iter, "time" = as.numeric(outMAP1$time, unit = "secs")),
      c("iter" = outMAP2$sol$iter, "time" = as.numeric(outMAP2$time, unit = "secs")),
      c("iter" = outMAP3$sol$iter, "time" = as.numeric(outMAP3$time, unit = "secs")))

(times_all[, 2]/times_all[, 1]) * 60

outMAP2$time
outMAP2$sol$iter

outMAP3$time
outMAP3$sol$iter

# Mu values
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


plot_96SBSsig(match12$R_hat) + plot_96SBSsig(match12$R_true)


plot_96SBSsig(match12$R_hat) + plot_96SBSsig(match12$R_true)

################################################
# Confrontation between CompNMF, SignatureAnalyzer, and SigPoisProcess
mutMatrix <- getTotalMutations(data$gr_Mutations)

set.seed(10)
BaselineNMF <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                              K = 30,
                                              alpha = 1.01, a = 1.01)

set.seed(10)
out_SigAnalyzerL1KL <- sigminer::sig_auto_extract(nmf_matrix = t(mutMatrix), method = "L1KL", cores = 20)
set.seed(10)
out_SigAnalyzerL1W.L2H <- sigminer::sig_auto_extract(nmf_matrix = t(mutMatrix), method = "L1W.L2H", cores = 20)

plot_sigs(BaselineNMF$Signatures)
match_to_RefSigs(BaselineNMF$Signatures, SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37)
plot_sigs(out_SigAnalyzerL1KL$Signature.norm)
match_to_RefSigs(out_SigAnalyzerL1KL$Signature.norm, SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37)
plot_sigs(out_SigAnalyzerL1W.L2H$Signature.norm)

saveRDS(BaselineNMF, file = "~/SigPoisProcess/output/Application_Breast/BaselineNMF/out_CompNMF.rds")
saveRDS(out_SigAnalyzerL1KL, file = "~/SigPoisProcess/output/Application_Breast/BaselineNMF/out_SigAnalyzerL1KL.rds")
saveRDS(out_SigAnalyzerL1W.L2H, file = "~/SigPoisProcess/output/Application_Breast/BaselineNMF/out_SigAnalyzerL1WL2H.rds")

devtools::document()
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


match_to_RefSigs(SigsMean, SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37)
match_to_RefSigs(match_w_baseline$R_hat + 1e-10, SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37)
match_to_RefSigs(match_w_sigAnalyzer$R_hat, SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37)


plot_96SBSsig(match12$R_hat) + plot_96SBSsig(match12$R_true) + plot_96SBSsig(match23$R_hat)
plot_96SBSsig(match12$R_hat) + plot_96SBSsig(match12$R_true)
plot_96SBSsig(match12$R_hat) + plot_96SBSsig(match12$R_true)



################################################
# Supervised analysis

#------------------------------------------------------------------
# Finally, plot the last application.
Sigs_to_use <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8", "SBS20",
                 "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30")
CosmicSigs <- SigPoisProcess::COSMIC_v3.4_SBS96_GRCh37[, Sigs_to_use]

# Load the estimates
outMCMCFixed <- readRDS("~/SigPoisProcess/output/Application_Breast/SemiSupervised/MCMCSolution_SemiSupervised_noMap.rds.gzip")
SigsCisFix <- get_PosteriorCI(outMCMCFixed$chain$SIGSchain[-c(1:2000), , ])
ThetaCisFix <- get_PosteriorCI(outMCMCFixed$chain$THETAchain[-c(1:2000), , ])
BetasCisFix <- get_PosteriorCI(outMCMCFixed$chain$BETASchain[-c(1:2000), , ])
MuFix <- colMeans(outMCMCFixed$chain$MUchain[-c(1:2000), ])
names(MuFix)  <- colnames(CosmicSigs)
colnames(BetasCisFix$mean) <- colnames(BetasCisFix$lowCI) <- colnames(BetasCisFix$highCI) <- colnames(CosmicSigs)
rownames(ThetaCisFix$mean) <- colnames(CosmicSigs)

# Plot the betas
p_Beats_sigs_fixed <- plot_betas_capped(BetasCisFix$mean, BetasCisFix$lowCI, BetasCisFix$highCI) + theme(legend.position = "none")

# Plot the fraction of the mutations attributed to the signatures
bigProd <- SigsCisFix$mean[data$gr_Mutations$channel, ] * t(ThetaCisFix$mean[, data$gr_Mutations$sample]) *
  exp(as.matrix(mcols(data$gr_Mutations)[, -c(1:3)]) %*% BetasCisFix$mean)
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
  left_join(data.frame(best_sig = names(MuFix),
                       mu = MuFix), by = "best_sig")

p_mu_sigs_fixed <- plot_vector_facets_x(df_AssignFix, levs = colnames(BetasCisFix$mean))
p_mu_sigs_fixed + p_Beats_sigs_fixed + patchwork::plot_layout(widths = c(0.5,2))
p_Beats_sigs_fixed + p_mu_sigs_fixed + patchwork::plot_layout(widths = c(2, 0.5))
# ggsave("~/SigPoisProcess/figures/FixedBetas.pdf", , width = 9.40, height = 6.42)
ggsave("~/SigPoisProcess/figures/FixedBetas.pdf", , width = 8.42, height = 5.56)

#----------------------------------- Plot now the signal as it varies


# Calculate the BaseNMF
mutMatrix <- getTotalMutations(data$gr_Mutations)
BaseNMF <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                              K = 0,
                                              S = 1e6 * CosmicSigs + 1,
                                              alpha = 1.01, a = 1.01)
dimnames(BaseNMF$mapOutput$R) <- dimnames(CosmicSigs)
colnames(BaseNMF$mapOutput$Theta) <- colnames(mutMatrix)
rownames(BaseNMF$mapOutput$Theta) <- colnames(CosmicSigs)

# Make the plot now.
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

Rchain <- outMCMCFixed$chain$SIGSchain[-c(1:7000), ,]
THETAchain <- outMCMCFixed$chain$THETAchain[-c(1:7000), ,]
BETASchain <- outMCMCFixed$chain$BETASchain[-c(1:7000), ,]

calculate_ChannelProbs_regionCI <- function(chr,
                                            region_start, region_end,
                                            Rchain,
                                            THETAchain,
                                            BETASchain,
                                            i = 1,
                                            j = 1,
                                            data) {

  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, data$gr_SignalTrack)
  gr_subset <- data$gr_SignalTrack[subjectHits(overlaps)]
  # Estimate the loadings in each subset
  SignalTrack_subset <- as.matrix(mcols(gr_subset)[, -1])
  CopyTrack_subset <- data$CopyTrack[subjectHits(overlaps), j]
  # Calculate the probability for the sequences
  nsamples <- dim(Rchain)[1]
  res_all <- array(NA, dim = c(nsamples, nrow(SignalTrack_subset), dim(Rchain)[3]))
  pb <- txtProgressBar(min = 0, max = nsamples, style = 3)
  for(s in 1:nsamples){
    R <- Rchain[s, ,]
    Theta <- THETAchain[s, ,]
    Betas <- BETASchain[s, ,]
    ExpSolBetas <- exp(SignalTrack_subset %*% Betas)
    ExpSolBetas_copy <- apply(ExpSolBetas, 2, function(x) x * CopyTrack_subset)
    SigTheta <- t(R[i, ] * Theta[, j])
    res_all[s, ,] <- SigTheta[rep(1, nrow(ExpSolBetas_copy)), ] * ExpSolBetas_copy
    setTxtProgressBar(pb, s)
  }
  close(pb)

  lowCi <- apply(res_all, c(2,3), function(x) quantile(x, 0.025))
  highCi <- apply(res_all, c(2,3), function(x) quantile(x, 0.975))
  MeanAll <- apply(res_all, c(2,3), mean)

  return(list("MeanAll" = MeanAll, "lowCi" = lowCi, "highCi" = highCi))
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


i <- "A[C>A]C"
j <- "DO220823"#"DO218709" #"DO1006"
chr <- "chr1"
# region_start <- 1e6
# region_end   <- 22e6
region_start <- 150e6#140e6
region_end   <- 160e6#2100e6

ChannelProbsCI <- calculate_ChannelProbs_regionCI(chr,region_start, region_end,
                                                  Rchain, THETAchain, BETASchain,
                                                  i, j, data)

ChannelProbs <- calculate_ChannelProbs_region(chr,
                                              region_start,
                                              region_end,
                                              SigsCisFix$mean,
                                              ThetaCisFix$mean,
                                              BetasCisFix$mean,
                                              j = j,
                                              data$gr_SignalTrack,
                                              data$CopyTrack)

# Constant value for the NMF
pr_const <- 2000 *BaseNMF$mapOutput$R[i, ] * BaseNMF$mapOutput$Theta[, j]/(sum(data$gr_SignalTrack$bin_weight))
df_const <- data.frame(t(pr_const)[rep(1,nrow(ChannelProbs)), ]) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1)) * 2000 + region_start,
         method = "Standard NMF")

df_mod <- data.frame(ChannelProbs[, , i])
colnames(df_mod) <- colnames(BaseNMF$mapOutput$R)

df_pos <- find_mutations_in_region(data, i, j, chr,
                                   region_start, region_end,
                                   SigsCisFix$mean,
                                   ThetaCisFix$mean,
                                   BetasCisFix$mean)

window <- 1
p_intensity <- df_mod %>%
  mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
                                                   align = "center"))) %>%
  mutate(region = (0:(nrow(df_mod) - 1))*2000 + region_start,
         method = "PoissonProcess") %>%
  bind_rows(df_const) %>%
  gather(key = "Sig", value = "Prob", - region, -method) %>%
  filter(Sig %in% c("SBS3", "SBS5", "SBS8")) %>%#head(names(sort(pr_const, decreasing = TRUE)), 4))%>%
  ggplot() +
  geom_line(aes(x = region, y = Prob, color = Sig, linetype = method, alpha = method), linewidth = 0.6) +
  scale_color_manual(name = "Signature",
                     values = c("red", "antiquewhite3", "blue")) +
  scale_alpha_manual(values = c(0.5, 1))+
  ylab("Expected mutations in 2kb region") +
  #xlab(paste0("Genomic regions in ",  chr, "-",
  #            region_start,":", region_end, " (2kb)")) +
  #xlab(paste0("Genomic regions in ",  chr, " (2kb)")) +
  facet_wrap(~paste0("Mutation type ", i, " - Patient ", j)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  #geom_vline(xintercept = lines_to_add, color = "black", linewidth = 0.6)
  geom_rug(data = df_pos,
           aes(x = positions, color = Sig),
           sides = "b", # "b" for bottom
           inherit.aes = FALSE,
           length = grid::unit(0.035, "npc"),
           linewidth = 1.2)+
  theme(axis.title.x = element_blank())
p_intensity


#--------
# df_modMean <- data.frame(ChannelProbsCI$MeanAll)
# df_modLowCi <- data.frame(ChannelProbsCI$lowCi)
# df_modHighCi <- data.frame(ChannelProbsCI$highCi)
# colnames(df_modMean) <- colnames(df_modLowCi) <- colnames(df_modHighCi) <- colnames(BaseNMF$mapOutput$R)
#
# region <- (0:(nrow(df_modMean)-1))*2000 + region_start
# sig_keep <- c("SBS3", "SBS5", "SBS8")
#
# window <- 1
# df_tmp <- df_modMean %>%
#   mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
#                                                    align = "center"))) %>%
#   mutate(region = (0:(nrow(df_modMean) - 1))*2000 + region_start,
#          method = "PoissonProcess") %>%
#   gather(key = "Sig", value = "Prob", - region, -method) %>%
#   left_join(df_modLowCi %>%
#               mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
#                                                                align = "center"))) %>%
#               mutate(region = (0:(nrow(df_modLowCi) - 1))*2000 + region_start,
#                      method = "PoissonProcess") %>%
#               gather(key = "Sig", value = "lowCi", - region, -method),
#             by = c("region", "method", "Sig")) %>%
#   left_join(df_modHighCi %>%
#               mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
#                                                                align = "center"))) %>%
#               mutate(region = (0:(nrow(df_modHighCi) - 1))*2000 + region_start,
#                      method = "PoissonProcess") %>%
#               gather(key = "Sig", value = "highCi", - region, -method),
#             by = c("region", "method", "Sig"))
#
# # Plot with ggplot2:
# ggplot(df_tmp %>% filter(Sig %in% sig_keep), aes(x = region, color = Sig, fill = Sig)) +
#   geom_ribbon(aes(ymin = log(lowCi), ymax = log(highCi)), alpha = 0.5, color = NA) +
#   geom_line(aes(y = log(Prob)), size = 0.7) +
#   scale_color_manual(name = "Signature",
#                      values = c("SBS3"="red", "SBS5"="antiquewhite3", "SBS8"="blue")) +
#   scale_fill_manual(name = "Signature",
#                     values = c("SBS3"="red", "SBS5"="antiquewhite3", "SBS8"="blue")) +
#   ylab("Expected mutations in 2kb region") +
#   theme_minimal()

#--------
# Now, make the plot of the various histone marks below it.
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
#ggsave(plot = p_panels, filename = "~/SigPoisProcess/figures/Track_and_mutations.png", width = 9.40, height = 6.42)
ggsave(plot = p_panels, filename = "~/SigPoisProcess/figures/Track_and_mutations.png",
       width = 11.68, height = 5.91)
# # Now plot just the mutational signatures
# p_sigs <- CompressiveNMF::plot_SBS_signature(Sigs[, c("SBS1", "SBS3", "SBS5", "SBS8")])+
#   theme(axis.text.x = ggplot2::element_blank())
# colnames(outFixed$Betas) <- colnames(Sigs)
# p_betas <- plot_betas(outFixed$Betas[, c("SBS1", "SBS3", "SBS5", "SBS8")])+
#   theme(axis.text.y = ggplot2::element_blank())
# p_sigs + p_betas+ plot_layout(widths = c(3.5, 2))
#


# How many mutations change?
Compute_mutation_Probs <- function(gr_Mutations, R, Theta, Betas){
  X <- as.matrix(mcols(gr_Mutations)[, -c(1:3)])
  Probs <- R[gr_Mutations$channel, ] * t(Theta[, gr_Mutations$sample]) * exp(X %*% Betas)
  Probs <- t(apply(Probs, 1, function(x) x/sum(x)))
  colnames(Probs) <- colnames(Betas)
  return(Probs)
}

# Poisson Process
MutProbs <- Compute_mutation_Probs(data$gr_Mutations,
                                   R = SigsCisFix$mean,
                                   Theta = ThetaCisFix$mean,
                                   Betas = BetasCisFix$mean)



# BaselineNMF
Betas_zero <- matrix(0, nrow = nrow(BetasCisFix$mean), ncol = ncol(BetasCisFix$mean),
                     dimnames = dimnames(BetasCisFix$mean))
MutProbs_fixed <- Compute_mutation_Probs(data$gr_Mutations,
                                        R = BaseNMF$mapOutput$R,
                                        Theta = BaseNMF$mapOutput$Theta,
                                        Betas = Betas_zero)

best_sigs_ppf <- Sigs_to_use[apply(MutProbs, 1, which.max)]
best_sigs_ppf <- factor(best_sigs_ppf, levels = Sigs_to_use)

best_sigs_nmf <- Sigs_to_use[apply(MutProbs_fixed, 1, which.max)]
best_sigs_nmf <- factor(best_sigs_nmf, levels = Sigs_to_use)

df_probs %>%
  mutate(channel = rownames(MutProbs)) %>%
  filter(Best_nmf != Best_ppf) %>%
  group_by(channel, Best_nmf, Best_ppf) %>%
  summarize(n = n()) %>%
  filter(Best_nmf %in% c("SBS8", "SBS3", "SBS5"),
         Best_ppf %in% c("SBS8", "SBS3", "SBS5"))

plot_sigs(CosmicSigs[, c("SBS8", "SBS3", "SBS5")])

df_probs <- data.frame("Best_ppf" = best_sigs_ppf,
                       "Prob_ppf" = apply(MutProbs, 1, max),
                       "Best_nmf" = best_sigs_nmf,
                       "Prob_nmf" = apply(MutProbs_fixed, 1, max))

filter_pp <- apply(MutProbs, 1, max)  > 0.5
diff_max <- apply(MutProbs,1, function(x) max(x) - max(x[-which.max(x)]))

max(MutProbs[1, ]) - max(MutProbs[1, ][-which.max(MutProbs[1, ])])

p_confusion <-df_probs %>%
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


p_mu_sigs_fixed
ggsave("~/SigPoisProcess/figures/Breast_SemiSup_mu.pdf",
       width = 11.68, height = 5.91)

(p_Beats_sigs_fixed + theme(aspect.ratio = 1)) + p_confusion #+ patchwork::plot_layout(widths = c(0.1,2,2))
ggsave("~/SigPoisProcess/figures/Breast_SemiSup_Beta_confusion.pdf",
       width = 11.68, height = 5.05)

# ThetaCopy <- ThetaCisFix$mean * t(colSums(data$CopyTrack))[rep(1, 13), ]
# Thetas <- apply(ThetaCopy, 2, function(x) x/sum(x))
# Thetas <- Thetas[, names(sort(colSums(ThetaCopy)))]
# Wmat <- ThetaCopy[, names(sort(colSums(ThetaCopy)))]
#
# # groups <- 1*(colnames(Wmat) %in% names(sort(colSums(Wmat)))[1:109])
# # groups <- ifelse(groups == 1, "1.Normal", "2.Hypermutated")
# factoextra::fviz_nbclust(t(Thetas), cluster::pam, method = "silhouette")
# clust <- cluster::pam(t(Thetas), k = 3, stand = FALSE)$clustering
#
# col_values <- c("darkblue", "#4959C7",  "skyblue1", "#F2AB67",  "#F6C866",  "#ED7470",
#                 "brown" , "#999999",  "black", "#b5de2b", "#6ece58", "#35b779", "#1f9e89")
#
# pup <- plot_baseline_weights(Thetas, clust = clust) +
#   scale_fill_manual(values = col_values)
# pdown <- plot_baseline_weights(Wmat, clust = clust) +
#   scale_fill_manual(values = col_values)
# pup/pdown








table(best_sigs_nmf, best_sigs_ppf)

filter_pp <- apply(MutProbs, 1, max)  > 0.5
table(filter_pp)
table(best_sigs_nmf[filter_pp], best_sigs_ppf[filter_pp])

hist(apply(MutProbs_fixed, 1, max))
plot_betas()



# Suppose your matrix is 'tab'
tab <- table(best_sigs_nmf[filter_pp], best_sigs_ppf[filter_pp])

# Convert to tidy data frame, zeroes to NA
df <- as.data.frame(tab)
colnames(df) <- c("nmf", "ppf", "count")
df$count[df$count == 0] <- NA

# Keep row and col order for axes
df$nmf <- factor(df$nmf, levels = rownames(tab))
df$ppf <- factor(df$ppf, levels = colnames(tab))

# Plot
ggplot(df, aes(x = ppf, y = nmf)) +
  geom_point(aes(size = count), color = "steelblue", na.rm = TRUE) +
  scale_size_continuous(range = c(2, 15)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(df$nmf))) +
  theme_bw() +
  labs(x = "best_sigs_ppf", y = "best_sigs_nmf", size = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        aspect.ratio = 1)

print(p)

library(ggplot2)
library(dplyr)

# Suppose your matrix is 'tab'
tab <- table(best_sigs_nmf, best_sigs_ppf)

# Convert to tidy data frame, zeroes to NA
df <- as.data.frame(tab)
colnames(df) <- c("nmf", "ppf", "count")

df$count[df$count == 0] <- NA

# Keep row/col order for axes
df$nmf <- factor(df$nmf, levels = rownames(tab))
df$ppf <- factor(df$ppf, levels = colnames(tab))

# Create an alternating color column for tiles
df <- df %>%
  mutate(AltCol = ifelse((as.numeric(ppf) %% 2) == 0, "gray93", "gray97"))

ggplot(df, aes(x = ppf, y = nmf)) +
  geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
  scale_fill_manual(values = c("gray93", "gray97"), guide = "none") +
  geom_point(aes(size = count), color = "steelblue", na.rm = TRUE) +
  scale_size_binned(
    range = c(0.5, 10),
    breaks = c(1, 10, 100, 1000, 5000, 10000, 40000, 100000, 140000), # Choose sensible breaks for your data
    name = "Count"
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(df$nmf))) +
  theme_bw() +
  labs(x = "best_sigs_ppf", y = "best_sigs_nmf", size = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        panel.grid = element_blank())
print(p)


#---- OLD PLOTS for GENES, not to use in the paper


# Extract posterior means for all quantities
mu <- colMeans(outMCMC$chain$MUchain[-c(1:6500), ])
R <- apply(outMCMC$chain$SIGSchain[-c(1:6500), ,], c(2, 3), mean)[, mu > 0.01]
Theta <- apply(outMCMC$chain$THETAchain[-c(1:6500), ,], c(2, 3), mean)[mu > 0.01, ]
Betas <- apply(outMCMC$chain$BETASchain[-c(1:6500), ,], c(2, 3), mean)[, mu > 0.01]

# Find the location of a known oncogene


get_gene_gr <- function(gene_symbol) {
  entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=gene_symbol, keytype="SYMBOL",
                                  columns="ENTREZID")$ENTREZID[1]
  gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)[entrez]
  gr
}

gr_gene <- get_gene_gr("KRAS")

# find overlap with gr_signalTrack and gr_CopyTrack
overlaps_track <- findOverlaps(gr_gene, gr_SignalTrack_2kb_std)
gr_gene_track <- gr_SignalTrack_2kb_std[subjectHits(overlaps_track)]

overlaps_copy <- findOverlaps(gr_gene, gr_CopyTrack)
gr_gene_copy <- gr_CopyTrack[subjectHits(overlaps_copy)]

# plot(gr_gene_track$GC)
# plot(gr_gene_track$CTCF)
# plot(gr_gene_track$H3K27ac)
# plot(gr_gene_track$H3K4me1)
# plot(gr_gene_track$H3K4me3)
# plot(gr_gene_track$RepliTime)
# plot(gr_gene_track$NuclOccup)

# Now, calculate the average activity for each signature
copyGene <- as.matrix(mcols(gr_gene_copy)[, -1]) * 2000

TrackSigs <- exp(as.matrix(mcols(gr_gene_track)[, -1]) %*% Betas) * tcrossprod(copyGene, Theta)

get_long_tracksigs <- function(TrackSigs, bin_size) {
  n_bins <- nrow(TrackSigs)
  n_cols <- ncol(TrackSigs)
  starts <- seq(1, by=bin_size, length.out=n_bins)
  ends <- starts + bin_size - 1
  df_wide <- data.frame(position = as.vector(t(cbind(starts, ends))))
  for(col in 1:n_cols) {
    df_wide[[colnames(TrackSigs)[col]]] <- rep(TrackSigs[, col], each = 2)
  }
  return(df_wide)
}

df_TrackSig <- get_long_tracksigs(TrackSigs, bin_size)
df_TrackSig %>%
  gather(key = "signature", value = "value", -position) %>%
  ggplot() +
  geom_line(aes(x = position, y = value, colour = signature))
















