# This script merges covariates at the 1kb resolution with the patients.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(patchwork)

# Load genome
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:24]
gr_tumor <- readRDS(file = "~/SigPoisProcess/data/ICGC_snp_mutations/Liver-HCC_snp.rds.gzip")
blacklist <- rtracklayer::import("~/SigPoisProcess/data/ENCODE_pvalues/hg19-blacklist.v2.bed")
gr_tumor <- gr_tumor[-queryHits(findOverlaps(gr_tumor, blacklist))]

hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
hg19_1kb <- tileGenome(chrom_lengths, tilewidth = 1000, cut.last.tile.in.chrom = TRUE)

df_tumor_all <- as.data.frame(gr_tumor) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n())

df_corrs <- data.frame()

for(i in 1:length(hg19_Mb)){
  print(i)
  region_filter <- hg19_Mb[i]
  # Load skin melanoma
  overlaps <- findOverlaps(gr_tumor, region_filter)
  gr_tumor_filter <- gr_tumor[queryHits(overlaps)]

  if(length(gr_tumor_filter) > 0){
    df_tumor <- as.data.frame(gr_tumor_filter)
    regions_all <- subjectHits(findOverlaps(region_filter, hg19_1kb))

    overlaps_regions <- findOverlaps(gr_tumor_filter, hg19_1kb)
    df_tumor$region <- subjectHits(overlaps_regions)
    table(df_tumor$sample)

    df_aggr <- df_tumor %>%
      group_by(sample, region) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      tidyr::complete(sample, region = regions_all, fill = list(n = 0)) %>%
      group_by(region) %>%
      summarise(all = sum(n))

    # Load signaltrack
    #signalCTCF <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates_normalized/Skin-Melanoma_male_CTCF.bigWig", which = region_filter)
    signalCTCF <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_female_CTCF.bigWig", which = region_filter)
    df_aggr$ctcf <- signalCTCF$score

    #signalH3K9me3 <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates_normalized/Skin-Melanoma_male_H3K9me3.bigWig", which = region_filter)
    signalH3K9me3 <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_female_H3K9me3.bigWig", which = region_filter)
    df_aggr$H3K9me3 <- signalH3K9me3$score

    #signalH3K27me3 <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates_normalized/Skin-Melanoma_male_H3K27me3.bigWig", which = region_filter)
    #df_aggr$H3K27me3 <- signalH3K27me3$score

    signalH3K27ac <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_male_H3K27ac.bigWig", which = region_filter)
    df_aggr$H3K27ac <- signalH3K27ac$score

    #signalRepli <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates_normalized/GC_male.bigWig", which = region_filter)
    signalDNAse<- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_male_DNAse-seq.bigWig", which = region_filter)
    df_aggr$DNAse_seq <- signalDNAse$score

    signalMethyl<- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Liver-HCC_Methylation.bigWig", which = region_filter)
    df_aggr$Methyl <- signalMethyl$score

    signalRepli <- import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig", which = region_filter)
    df_aggr$Repli <- signalRepli$score

    cors <- cor(df_aggr[, -1])[1, 2:7]
    df_corrs <- rbind(df_corrs,
                      data.frame(region = i,
                                t(cors)))
  }
}

signalRepli <- import("~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_1kb/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig")
overlaps <- findOverlaps(gr_tumor, signalRepli)
sc <- rep(1e-18, length(gr_tumor))
sc[queryHits(overlaps)] <- signalRepli$score[subjectHits(overlaps)]
gr_tumor$score <- sc
hist(sc, breaks = 100)

plot(diff(gr_tumor[gr_tumor$sample == "DO45197"]$score), type = "l")

df_corrs <- read_tsv("tmp_corss_liver.tsv")

df_corrs %>%
  left_join(df_tumor_all, by = "region") %>%
  filter(n > 500, H3K9me3 > 0.4) %>%
  arrange(desc(H3K9me3)) %>%
  head(20)

i <- 2988 # <------ Very important

df_corrs %>%
  left_join(df_tumor_all, by = "region") %>%
  filter(Repli > 0.4) %>%
  arrange(desc(Repli)) %>%
  head(20)


# Mutations in a Macro-region
p1 <- ggplot() +
  theme_bw(base_size = 11) +
  geom_rect(data = as.data.frame(hg19_Mb) %>%
              mutate(region = 1:length(hg19_Mb)) %>%
              group_by(seqnames) %>%
              summarise(xmin = min(region), xmax = max(region)) %>%
              arrange(seqnames) %>%
              mutate(fill_color = factor(seq_along(seqnames) %% 2)),
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = fill_color),
            alpha = 0.2) +
  scale_fill_manual(values = c("#4682B4", "antiquewhite2" )) +
  geom_point(data = df_tumor_all,
             aes(x = region + 0.5, y = n), size = 1, shape = 21, color = "gray25") +
  geom_point(data = df_tumor_all %>%
               filter(region == 2988),
             aes(x = region + 0.5, y = n), size = 2, color = "red") +
  #ylab("Total mutations\n") +
  xlab("Genomic region (1Mb)") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  facet_wrap(.~"Total mutations in Liver-HCC (Mb scale)")



df_aggr$region2 <- 90000001 + 1000 * (0:999)
p2 <- ggplot(df_aggr %>%
               mutate(is_zero = (all == 0)), aes(x = region2 + 0.5, y = all)) +
  geom_segment(aes(x = region2 + 0.5, xend = region2 + 0.5, y = 0, yend = all),
               linewidth = 0.2, color = "#4682B4", alpha = 0.7) +
  geom_point(aes(color = is_zero,shape = is_zero), size = 1.2)+
  scale_color_manual(values = c("#4682B4", "gray25")) +
  theme_bw(base_size = 11) +
  scale_shape_manual(values = c(1, NA)) +
  theme(legend.position = "none", axis.title = element_blank()) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  #facet_wrap(.~"Total mutations in chrX:90000001-91000000, (Kb scale)")+
  #ylab("Total mutations\n")+
  xlab("Genomic region (1kb)")

# Let's now make the histome modification H3K9me3 example
# p3 <- ggplot(df_aggr) +
#   geom_area(aes(x = region2, y = H3K9me3), alpha = 0.9, fill = "#4682B4") +
#   theme_minimal(base_size = 11) +
#   theme(axis.title = element_blank())+
#   scale_x_continuous(expand = c(0.01, 0.02)) +
#   scale_y_continuous(expand = c(0.03, 0.02))

p3 <- df_aggr %>%
  mutate(region2_next = lead(region2)) %>%
  filter(!is.na(region2_next)) %>%
  rowwise() %>%
  do(tibble(region2 = c(.$region2, .$region2_next), H3K9me3 = rep(.$H3K9me3, 2))) %>%
  ungroup() %>%
  ggplot()+
  geom_area(aes(x = region2, y = H3K9me3), alpha = 0.9, fill = "#4682B4") +
  theme_bw(base_size = 11) +
  theme(axis.title = element_blank())+
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02))

p4 <- df_aggr %>%
  mutate(region2_next = lead(region2)) %>%
  filter(!is.na(region2_next)) %>%
  rowwise() %>%
  do(tibble(region2 = c(.$region2, .$region2_next), ctcf = rep(.$ctcf, 2))) %>%
  ungroup() %>%
  ggplot() +
  geom_area(aes(x = region2, y = ctcf), alpha = 0.9, fill = "red") +
  theme_minimal(base_size = 11) +
  theme(axis.title= element_blank())+
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02))

p1
ggsave(plot = p1, filename = "~/SigPoisProcess/figures/Total_mutations_LiverHcc.pdf", width = 9.16, height = 2.27)
(p2/p3/p4) + plot_layout(heights = c(1.5, 1.5, 1))

(p2/p3) + plot_layout(heights = c(1, 1))
ggsave(filename = "~/SigPoisProcess/figures/kb_mutations_H3K9me3_LiverHcc.pdf", width = 8.43, height = 2.31)

#### Use now replication timing
i <- 340
ggplot(df_aggr %>%
         mutate(is_zero = (all == 0)), aes(x = region + 0.5, y = all)) +
  geom_segment(aes(x = region + 0.5, xend = region + 0.5, y = 0, yend = all),
               linewidth = 0.2, color = "#4682B4", alpha = 0.6) +
  geom_point(aes(color = is_zero,shape = is_zero), size = 1.2)+
  scale_color_manual(values = c("#4682B4", "gray25")) +
  theme_bw(base_size = 11) +
  scale_shape_manual(values = c(1, NA)) +
  theme(legend.position = "none", axis.title = element_blank()) +
  scale_x_continuous(expand = c(0.03,0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  facet_wrap(.~"Total mutations in chrX:90000001-91000000, (Kb scale)")+
  #ylab("Total mutations\n")+
  xlab("Genomic region (1kb)")

ggplot(df_aggr) +
  geom_area(aes(x = region, y = Repli), alpha = 0.8, fill = "darkblue") +
  theme_bw(base_size = 11) +
  theme(axis.title= element_blank())+
  scale_x_continuous(expand = c(0.03, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02))

###### Plot for the output
TumorData <- readRDS("~/SigPoisProcess/data/ICGC_snp_merged/Liver-HCC_merged_1kb.rds.gzip")
out <- readRDS("~/SigPoisProcess/output/ICGC/Liver-HCC_MLE_init_v3.rds.gzip")
probs <- compute_SignatureCovariate_probs(out, gr_Mutations = TumorData$gr_tumor)

#plot(out)

genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:24]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)

i <- 2988
region_filter <- hg19_Mb[i]
overlaps <- findOverlaps(TumorData$gr_tumor, region_filter)
gr_tumor_filter <- TumorData$gr_tumor[queryHits(overlaps)]
probs_filter <- probs[queryHits(overlaps), ,]

df_sel <- as.data.frame(gr_tumor_filter) %>%
  group_by(sample) %>%
  summarize(n = n()) %>%
  as.data.frame()

patients <- df_sel$sample

pat <- 235

pp <- function(pat){
  probs_patient <- probs_filter[grepl(patients[pat], rownames(probs_filter)), , ]

  df_plot <- data.frame()
  for(j in 1:nrow(probs_patient)){
    df_temp <- reshape2::melt(probs_patient[j, ,])
    nm <- c(str_split(rownames(probs_patient)[j], "-", simplify = TRUE))
    df_temp$sample <- nm[1]
    df_temp$pos <- nm[2]
    df_temp$mut <- nm[3]
    df_temp$pos2 <- paste0(df_temp$pos, "\n", df_temp$mut)
    df_plot <- rbind(df_plot, df_temp)

  }

  files_all <- (TumorData$df_files %>%
                  filter(sample == patients[pat]) %>%
                  c() %>%
                  unlist())[-c(1,2)]

  df_tracks <- data.frame()
  for(s in 1:length(files_all)) {
    gr_track <- import(files_all[s], which = region_filter)
    df_temp <- as.data.frame(gr_track)
    df_temp$covariate <- names(files_all)[s]
    df_tracks <- rbind(df_tracks, df_temp)
  }
  df_tracks$covariate[df_tracks$covariate == "DNAse-seq"] <- "DNAse.seq"
  df_tracks$covariate <- factor(df_tracks$covariate, levels = levels(df_plot$Var2))


  df_alphabet <- data.frame("Var1" = levels(df_plot$Var1),
                            "Var3" = LETTERS[1:length(levels(df_plot$Var1))])

  p_probs <- ggplot(df_plot %>%
                      left_join(df_alphabet, by = "Var1")) +
    geom_tile(aes(x = Var3, y = Var2, fill = value), colour = "gray") +
    facet_wrap(.~ pos2, nrow = 1, strip.position = "bottom") +
    scale_fill_gradient(name = "Signature-\ncovariate\nprobability", low = "white", high = "#1F6097") +
    theme_bw(base_size = 9)  +
    scale_x_discrete(position = "top", guide = guide_axis(n.dodge = 2.5)) +
    theme(axis.title = element_blank())


  pos <- sort(unique(as.numeric(gsub(".*:", "", df_plot$pos))))
  p_tracks <- ggplot(df_tracks) +
    geom_raster(aes(x = start + 500, y = covariate, fill = score)) +
    #facet_wrap(.~ pos2, nrow = 1) +
    scale_fill_gradientn(name = "Signal track\n(normalized)",
                         colours =c("antiquewhite", "#F6C866", "#F2AB67", "#EF8F6B",
                                    "#ED7470", "#BF6E97", "#926AC2",
                                    "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
    theme_minimal(base_size = 9) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(xintercept = pos, col = "red", alpha = 0.5) +
    theme(axis.title = element_blank())


  p_probs/p_tracks
}

pp(235)
# ggsave("~/SigPoisProcess/figures/LiverHcc_track_235_map.pdf",
#        width = 10.91, height= 5.94)
ggsave("~/SigPoisProcess/figures/LiverHcc_track_235_mle.pdf",
       width = 10.91, height= 5.94)


probs_patient <- probs_filter[grepl(patients[pat], rownames(probs_filter)), , ]

df_plot <- data.frame()
for(j in 1:nrow(probs_patient)){
  df_temp <- reshape2::melt(probs_patient[j, ,])
  nm <- c(str_split(rownames(probs_patient)[j], "-", simplify = TRUE))
  df_temp$sample <- nm[1]
  df_temp$pos <- nm[2]
  df_temp$mut <- nm[3]
  df_temp$pos2 <- paste0(df_temp$pos, "\n", df_temp$mut)
  df_plot <- rbind(df_plot, df_temp)

}

files_all <- (TumorData$df_files %>%
  filter(sample == patients[pat]) %>%
  c() %>%
  unlist())[-c(1,2)]

df_tracks <- data.frame()
for(s in 1:length(files_all)) {
  gr_track <- import(files_all[s], which = region_filter)
  df_temp <- as.data.frame(gr_track)
  df_temp$covariate <- names(files_all)[s]
  df_tracks <- rbind(df_tracks, df_temp)
}
df_tracks$covariate[df_tracks$covariate == "DNAse-seq"] <- "DNAse.seq"
df_tracks$covariate <- factor(df_tracks$covariate, levels = levels(df_plot$Var2))


df_alphabet <- data.frame("Var1" = levels(df_plot$Var1),
                          "Var3" = LETTERS[1:length(levels(df_plot$Var1))])

p_probs <- ggplot(df_plot %>%
                    left_join(df_alphabet, by = "Var1")) +
  geom_tile(aes(x = Var3, y = Var2, fill = value), colour = "gray") +
  facet_wrap(.~ pos2, nrow = 1, strip.position = "bottom") +
  scale_fill_gradient(name = "Signature-\ncovariate\nprobability", low = "gray95", high = "#1F6097") +
  theme_bw(base_size = 10)  +
  scale_x_discrete(position = "top", guide = guide_axis(n.dodge = 2.5)) +
  theme(axis.title = element_blank())


pos <- sort(unique(as.numeric(gsub(".*:", "", df_plot$pos))))
p_tracks <- ggplot(df_tracks) +
  geom_raster(aes(x = start + 500, y = covariate, fill = score)) +
  #facet_wrap(.~ pos2, nrow = 1) +
  scale_fill_gradientn(name = "Signal track\n(normalized)",
                       colours =c("antiquewhite", "#F6C866", "#F2AB67", "#EF8F6B",
                                  "#ED7470", "#BF6E97", "#926AC2",
                                  "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
  theme_minimal(base_size = 10) +
  scale_x_continuous(expand = c(0,0)) +
  geom_vline(xintercept = pos, col = "red", alpha = 0.5) +
  theme(axis.title = element_blank())


p_probs/p_tracks


#
#
#
#
#
# table(df_aggr$sample)
#
# table(df_aggr$n)
#
# ggplot(df_aggr %>% filter(sample == "DO220906")) +
#   geom_line(aes(x = region, y = n))
#
# plot(df_aggr$region, df_aggr$mn)
# lines(df_aggr$region, df_aggr$low, type = "l", col = "red")
# lines(df_aggr$region, df_aggr$up, type = "l", col = "red")
#
#
# plot(df_aggr[df_aggr$sample == "DO220882", ]$region, df_aggr[df_aggr$sample == "DO220882", ]$n)
#
# # Load the Breast cancer data
# df_breast <- read_tsv("data/Breast_AdenoCa_icgc_snv.tsv")
# df_breast$Chromosome <- factor(paste0("chr", df_breast$Chromosome), levels = names(chrom_lengths))
# # Make breast cancer data into a genome
# grICGC_breast <- GRanges(
#   seqnames = df_breast$Chromosome,
#   IRanges(start = df_breast$Start_Position, end = df_breast$Start_Position),
#   sample = df_breast$Tumor_Sample_Barcode,
#   channel = df_breast$channel
# )
#
# hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
# # For visualization, count the number of mutations across genomic regions
# overlaps_1Mb <- findOverlaps(grICGC_breast, hg19_Mb)
# df_breast$region <- subjectHits(overlaps_1Mb)
#
# ggplot() +
#   geom_rect(data = df_breast %>%
#               group_by(Chromosome) %>%
#               summarise(xmin = min(region), xmax = max(region)) %>%
#               arrange(Chromosome) %>%
#               mutate(fill_color = factor(seq_along(Chromosome) %% 2)),
#             aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = fill_color),
#             alpha = 0.2) +
#   geom_point(data = df_breast %>%
#                group_by(Chromosome, region, type) %>%
#                summarise(nMuts = n()) %>%
#                filter(nMuts > 0),
#              aes(x = region, y = nMuts, color = type), size = 0.8, alpha = 0.5) +
#   theme_minimal() +
#   scale_fill_manual(values = c("gray", "white")) +
#   scale_color_manual(values = c("blue", "#020202", "red", "gray", "green", "pink")) +
#   theme(axis.title = element_blank()) +
#   facet_grid("Aggregate muts. (all patients)"~.)
#
#
# region_filter <- hg19_Mb[279]
#
#
#
# # Load signaltrack
# signalCTCF <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Skin-Melanoma_female_CTCF.bigWig", which = region_filter)
# plot(signalCTCF$score)
#
# signalH3K9me3 <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Skin-Melanoma_male_H3K9me3.bigWig", which = region_filter)
# plot(signalH3K9me3$score, type = "l")
#
# signalH3K27me3 <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Skin-Melanoma_female_H3K27me3.bigWig", which = region_filter)
# plot(signalH3K9me3$score)
#
# dnase <- import("~/SigPoisProcess/data/ENCODE_pvalues/Cancer_covariates/Skin-Melanoma_male_DNAse-seq.bigWig", which = region_filter)
# plot(dnase$score)
#
#
# cor(df_aggr$mn, signalH3K9me3$score)
#
#
#
#


# Is there a way to differential ADAR edited mutations and APOBEC mutations.
# High Apex1 or Apex2.
# slit high vs low, check signatuere exposure
#


tmp <- read_tsv("../data/E066-H3K9me3.imputed.narrowPeak.bed.nPk.gz", col_names = FALSE)

gImport <- GRanges(seqnames = tmp$X1,
                   IRanges(start = tmp$X2, end = tmp$X3))


tmp2 <- import("../data/E066-H3K9me3.imputed.pval.signal.bigwig")
tmp3 <- import("../data/E066-H3K9me3.imputed.pval.signal.bigwig", which = gImport[1])

table(10^-tmp3$score < 0.05)

width(gImport)






