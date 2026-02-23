################################################################################
# Reproduce Figure 1 in the paper - Breast cancer
################################################################################

#---- Load useful libraries
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(patchwork)

#---- Load the dataset and exclude the blacklisted regions
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
gr_tumor <- readRDS(file = "~/SigPoisProcess/data/data_for_figure1/Breast-AdenoCa_snp.rds.gzip")
blacklist <- rtracklayer::import("~/SigPoisProcess/data/data_for_figure1/hg19-blacklist.v2.bed")
gr_tumor <- gr_tumor[-queryHits(findOverlaps(gr_tumor, blacklist))]

#---- Partition the genome into Megabases, an into kilobases
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
hg19_2kb <- tileGenome(chrom_lengths, tilewidth = 2000, cut.last.tile.in.chrom = TRUE)

#---- Count the number of mutations in the dataset per megabase
df_tumor_all <- as.data.frame(gr_tumor) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n())

region_highlight <- 2986 # <--- we highlight this specific region

#-------------------------------- Figure 1 - Panel a
pTrack_genomeMb <- ggplot() +
  theme_bw(base_size = 11) +
  geom_rect(data = as.data.frame(hg19_Mb) %>%
              mutate(region = 1:length(hg19_Mb)) %>%
              group_by(seqnames) %>%
              summarise(xmin = min(region), xmax = max(region) + 1) %>%
              arrange(seqnames) %>%
              mutate(fill_color = factor(seq_along(seqnames) %% 2)),
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = fill_color),
            alpha = 0.2) +
  scale_fill_manual(values = c("#4682B4", "antiquewhite2" )) +
  geom_point(data = df_tumor_all,
             aes(x = region, y = n), size = 1, shape = 21, color = "gray25") +
  geom_point(data = df_tumor_all %>%
               filter(region == region_highlight),
             aes(x = region + 0.5, y = n), size = 2, color = "red") +
  xlab("Genomic region (1Mb)") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  facet_wrap(.~"Total mutations in Breast-AdenoCA (Mb scale)")

ggsave(plot = pTrack_genomeMb, filename = "~/SigPoisProcess/figures/Figure1_PanelA_Total_mutations_Breast.pdf", width = 9.16, height = 2.27)


#-------------------------------- Figure 1 - Panel b
region_filter <- hg19_Mb[region_highlight]
# Load the signal for the specific track
signalH3K9me3 <- import("~/SigPoisProcess/data/data_for_figure1/Breast-Cancer_tissue_H3K9me3_2kb.bigWig", which = region_filter)

# Overlaps between the tumor and the megabase region selected
overlaps <- findOverlaps(gr_tumor, region_filter)
gr_tumor_filter <-  gr_tumor_filter <- gr_tumor[queryHits(overlaps)]
# Overlaps between the region selected and the 2kb region
regions_all <- subjectHits(findOverlaps(region_filter, hg19_2kb))
overlaps_regions <- findOverlaps(gr_tumor_filter, hg19_2kb)

# Count the number of points per 2kb
df_aggr <- as.data.frame(gr_tumor_filter) %>%
  mutate(region = subjectHits(overlaps_regions)) %>%
  group_by(sample, region) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  tidyr::complete(sample, region = regions_all, fill = list(n = 0)) %>%
  group_by(region) %>%
  summarise(all = sum(n)) %>%
  mutate(region2 = start(region_filter) + 2000 * (0:499),
         H3K9me3 = signalH3K9me3$score)

#--- Top panel
panelB_top <- ggplot(df_aggr %>%
               mutate(is_zero = (all == 0)), aes(x = region2 + 0.5, y = all)) +
  geom_segment(aes(x = region2 + 0.5, xend = region2 + 0.5, y = 0, yend = all),
               linewidth = 0.2, color = "#4682B4", alpha = 0.7) +
  geom_point(aes(color = is_zero,shape = is_zero), size = 1.2)+
  scale_color_manual(values = c("#4682B4", "gray25")) +
  theme_bw(base_size = 11) +
  scale_shape_manual(values = c(1, NA)) +
  theme(legend.position = "none", axis.title = element_blank()) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  xlab("Genomic region (1kb)")

#--- Bottom panel
panelB_bottom <- df_aggr %>%
  ggplot()+
  geom_area(aes(x = region2, y = H3K9me3), alpha = 0.9, fill = "#4682B4") +
  theme_bw(base_size = 11) +
  theme(axis.title = element_blank())+
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02))

(panelB_top/panelB_bottom) + plot_layout(heights = c(1, 1))
ggsave(filename = "~/SigPoisProcess/figures/Figure1_PanelB_kb_mutations_H3K9me3_Breast.pdf", width = 8.43, height = 2.31)



#--------- NOT RUN --------- #
# Script to extract the correlation with between the total mutations
# in the megabase, and the signaltrack.

df_corrs <- data.frame()

for(i in 1:length(hg19_Mb)){
  print(i)
  region_filter <- hg19_Mb[i]
  overlaps <- findOverlaps(gr_tumor, region_filter)
  gr_tumor_filter <- gr_tumor[queryHits(overlaps)]

  if(length(gr_tumor_filter) > 0){
    df_tumor <- as.data.frame(gr_tumor_filter)
    regions_all <- subjectHits(findOverlaps(region_filter, hg19_2kb))

    overlaps_regions <- findOverlaps(gr_tumor_filter, hg19_2kb)
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
    signalH3K9me3 <- import("~/SigPoisProcess/data/data_for_figure1/Breast-Cancer_tissue_H3K9me3_2kb.bigWig", which = region_filter)
    df_aggr$H3K9me3 <- signalH3K9me3$score

    cors <- cor(df_aggr[, -1])[1, -1]
    df_corrs <- rbind(df_corrs,
                      data.frame(region = i,
                                 t(cors)))
  }
}

df_corrs %>%
  left_join(df_tumor_all, by = "region") %>%
  filter(n >= 200) %>%
  arrange(desc(t.cors.)) %>%
  head(20)

i <- 2986 # <------ We end up selecting this region as illustration


