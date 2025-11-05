library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)


#----
# Let's download all files for a selected group of cancer types
df_otlu <- read_tsv("~/SigPoisProcess/data/ENCODE_pvalues/Otlu_merged_with_experiments.tsv") %>%
  filter(`Cancer Type` == "Stomach-AdenoCA") %>%
  dplyr::select(`Topography Feature`, `File Name(s)`, target, biosample_summary) %>%
  distinct() %>%
  as.data.frame()

# Now, import the files
accessions <- gsub(".bed.gz", "", df_otlu$`File Name(s)`)

urls <- sprintf(
  "https://www.encodeproject.org/files/%s/@@download/%s.bed.gz",
  accessions, accessions
)
# Write all URLs into a text file
writeLines(urls, "~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bedfiles/files_Stomach.txt")
# On the terminal:
# nohup xargs -n 1 curl -O -L < files_Stomach.txt > download.log 2>&1 &

load_bedfile <- function(file){
  df <- read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  colnames(df) <- c("chrom", "chromStart", "chromEnd", "name", "score",
                    "strand", "signalValue", "pValue", "qValue", "peak")

  df <- df %>%
    mutate(chrom = factor(chrom,
                          levels = paste0("chr", c(1:22, "X", "Y")))) %>%
    arrange(chrom, chromStart, chromEnd) %>%
    drop_na()

  gr <- GRanges(
    seqnames = df$chrom,
    ranges = IRanges(
      start = df$chromStart + 1,   # convert 0-based BED start to 1-based GRanges start
      end = df$chromEnd
    ),
    strand = "*",
    name = df$name,
    score = df$score,
    signalValue = df$signalValue,
    pValue = df$pValue,
    qValue = df$qValue,
    peak = df$peak
  )
  return(gr)
}
#----

# Load gr_tumor
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")

# Remove mutations from blacklist
blacklist <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/hg19-blacklist.v2.bed")
overlaps <- findOverlaps(gr_tumor, blacklist)
gr_tumor <- gr_tumor[-queryHits(overlaps)]

# Extract the MAP solution
MutMatrix <- getTotalMutations(gr_tumor)
set.seed(10)
sigs_to_select <- paste0("SBS", c(1, 2, 3, 5, 13, 15, "17a", "17b", 18, 20, 21, 26, 28, "40a", 44))
Sigs <- 1e5 * CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37[, sigs_to_select] + 1
resMAP <- CompressiveNMF::CompressiveNMF_map(MutMatrix, S = Sigs, K = 0)


# We will take one file at a time, one covariate at a time, and understand
# how the mutations look like.
#----------------------------------------------------------------- H3K9me3
dir <- "~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bedfiles/"
files <- df_otlu[df_otlu$target == "H3K27ac", ]$`File Name(s)`

gr <- load_bedfile(paste0(dir, files[3]))

overlaps <- findOverlaps(gr_tumor, gr)
score_all <- (gr$score - min(gr$score))/(max(gr$score) - min(gr$score))
scores <- rep(0, length(gr_tumor))
scores[queryHits(overlaps)] <- (gr[subjectHits(overlaps)]$score - min(gr$score))/(max(gr$score) - min(gr$score))
gr_tumor$score <- scores

tmp <- data.frame(gr_tumor) %>%
  group_by(channel, sample) %>%
  filter(score > 0) %>%
  summarise(mn = sum(score)) %>%
  ungroup() %>%
  group_by(channel) %>%
  summarise(mn = mean(mn)) %>%
  column_to_rownames("channel") %>%
  as.matrix()
plot_96SBSsig(tmp)
table(scores>0)


mut_all <-as.matrix(rowSums(MutMatrix)); colnames(mut_all) <- "all"

plot_96SBSsig(mut_all)


#----------------------------------------------------------------- H3K27ac


#----------------------------------------------------------------- H3K27me3


#----------------------------------------------------------------- H3K36me3

# Load all files, then record the number of times a mutation falls in a peak
# for the files we have
gr_tumor_merged <- gr_tumor
dir <- "~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bedfiles/"
files <- df_otlu$`File Name(s)`
for(i in 1:length(files)){
  f <- files[i]
  gr <- load_bedfile(paste0(dir, f))
  overlaps <- findOverlaps(gr_tumor_merged, gr)
  score_all <- (gr$score - min(gr$score))/(max(gr$score) - min(gr$score))
  scores <- rep(0, length(gr_tumor_merged))
  scores[queryHits(overlaps)] <- (gr[subjectHits(overlaps)]$score - min(gr$score))/(max(gr$score) - min(gr$score))
  gr_tumor_merged$score <- scores
  colnames(mcols(gr_tumor_merged))[ncol(mcols(gr_tumor_merged))] <- paste0( df_otlu$target[i], "_", f)
}

counts <- colSums(as.matrix(mcols(gr_tumor_merged)[, -c(1:3)]) > 0)
# Extract type (before underscore)
types <- sub("_.*", "", names(counts))

# Calculate mean per type
sort(tapply(counts, types, mean))

plot_96SBSsig(resMAP$Signatures)
plot(MutMatrix, resMAP$Signatures %*% resMAP$Theta)

# Now, give a score to every mutation

Probs <- resMAP$Signatures[gr_tumor_merged$channel, ] * t(resMAP$Theta[, gr_tumor_merged$sample])
Probs <- t(apply(Probs, 1, function(x) x/sum(x)))
maxMut <- apply(Probs, 1, function(x) names(which.max(x)))
#mcols(gr_tumor_merged) <- cbind(mcols(gr_tumor_merged), Probs)

#plot(density(gr_tumor_merged$SBS15[gr_tumor_merged$H3K36me3_ENCFF569XKN.bed.gz > 0]))
#lines(density(gr_tumor_merged$SBS15[gr_tumor_merged$H3K36me3_ENCFF569XKN.bed.gz == 0]), col = "red")

# SBS15 all
ggplot(as.data.frame(gr_tumor_merged[Probs[, "SBS15"] > 0.01, ]) %>%
         mutate(H3K36me3 = (H3K36me3_ENCFF569XKN.bed.gz > 0))) +
  geom_boxplot(aes(x = H3K36me3, y = SBS15))

ggplot(as.data.frame(gr_tumor_merged[Probs[, "SBS15"] > 0.01, ]) %>%
         mutate(H3K36me3 = (H3K36me3_ENCFF569XKN.bed.gz > 0))) +
  geom_boxplot(aes(x = H3K36me3, y = SBS15))

ggplot(as.data.frame(gr_tumor_merged[Probs[, "SBS15"] > 0.01, ]) %>%
         mutate(H3K36me3 = (H3K36me3_ENCFF569XKN.bed.gz > 0))) +
  geom_density(aes(color = H3K36me3, SBS15))

plot(gr_tumor_merged$H3K36me3_ENCFF569XKN.bed.gz, gr_tumor_merged$SBS15)

gr_tumor_merged

# SBS15 max
ggplot(as.data.frame(gr_tumor_merged[maxMut == "SBS15", ]) %>%
         mutate(H3K36me3 = (H3K36me3_ENCFF569XKN.bed.gz > 0))) +
  geom_boxplot(aes(x = H3K36me3, y = SBS15))
plot(gr_tumor_merged$H3K36me3_ENCFF569XKN.bed.gz, gr_tumor_merged$SBS15)


tmp <- as.data.frame(gr_tumor_merged[maxMut == "SBS44", ]) %>%
  mutate(H3K36me3 = (H3K36me3_ENCFF473POV.bed.gz > 0))
t.test(SBS15 ~ H3K36me3, data = tmp, alternative = "less")

table(maxMut == "SBS15", gr_tumor_merged$H3K36me3_ENCFF473POV.bed.gz > 0)
sig <-  "SBS44"
counts <- colSums(as.matrix(mcols(gr_tumor_merged[maxMut == sig, ])[, -c(1:3, 49:60)]) > 0)
sort(tapply(counts/sum(maxMut == sig), types, mean))


sig <- "SBS1"
mark <- "CTCF"
plot_peak_prob <- function(sig, mark){
  p <- as.data.frame(gr_tumor_merged) %>%
    mutate(SigProb = Probs[, sig]) %>%
    filter(SigProb > 0.01) %>%
    gather(key = "peak", value = "value", -(seqnames:channel), -SigProb) %>%
    filter(grepl(mark, peak)) %>%
    mutate(value_binary = value > 0) %>%
    ggplot(aes(x = value_binary, y = SigProb, fill = value_binary)) +
    geom_boxplot(alpha = 0.6) +
    theme_bw()+
    #geom_jitter(size=0.2, alpha=0.5) +
    facet_grid(.~peak)+
    ggtitle(paste(sig, "and", mark))
  return(p)
}
plot_peak_prob("SBS1", "CTCF")
plot_peak_prob("SBS2", "H3K36me3")
plot_peak_prob("SBS2", "H3K27ac")
plot_peak_prob("SBS13", "H3K36me3")
plot_peak_prob("SBS13", "H3K27ac")
plot_peak_prob("SBS3", "CTCF")
plot_peak_prob("SBS3", "H3K9me3")
plot_peak_prob("SBS15", "H3K9me3")
plot_peak_prob("SBS15", "H3K36me3") #<---------!
plot_peak_prob("SBS15", "H3K27ac")
plot_peak_prob("SBS26", "H3K36me3")
plot_peak_prob("SBS26", "H3K9me3")
plot_peak_prob("SBS44", "H3K36me3")
plot_peak_prob("SBS44", "H3K36me3")


#----------------------------------------------------------------- GC content
df_tumor <- as.data.frame(gr_tumor)

gr <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_GC_1kb_standardized.bigWig.rds.gzip")
gr <- gr[gr$bin_weight > 0]

overlaps <- findOverlaps(gr_tumor, gr)
scores <- rep(0, length(gr_tumor))
scores[queryHits(overlaps)] <- gr[subjectHits(overlaps)]$score
plot(density(gr$score))
lines(density(scores), col = "red")

df_tumor$score <- scores
df_tumor$mut <- sub(".*\\[([A-Z]>[A-Z])\\].*", "\\1", df_tumor$channel)

ggplot() +
  theme_bw()+
  geom_density(data = df_tumor, aes(score, color = mut), linewidth = 0.8)+
  geom_density(data = data.frame(score= gr$score), aes(score),
               fill = "blue", alpha = 0.1, color = "blue")+
  scale_color_manual(values = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))

df_tmp <- rbind(data.frame(score = gr$score, mut = "AllTrack"),
                df_tumor[, c("score", "mut")])

ggplot() +
  theme_bw()+
  geom_boxplot(data = df_tmp, aes(x = mut, y = score, color = mut, fill = mut), alpha = 0.6, outlier.size = 0.5)+
  scale_fill_manual(values = c("white", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  scale_color_manual(values = c("black", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  geom_hline(yintercept = median(gr$score), linetype = "dashed")


#----------------------------------------------------------------- Methylation
df_tumor <- as.data.frame(gr_tumor)
gr <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_Methylation_1kb_standardized.bigWig.rds.gzip")
#gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/old_data/Stomach_AdenoCA/methylation_bigWig_GRCh38/Stomach_methylation_male_hg19_1kb.bigWig")
gr <- gr[gr$bin_weight > 0]

overlaps <- findOverlaps(gr_tumor, gr)
scores <- rep(0, length(gr_tumor))
scores[queryHits(overlaps)] <- gr[subjectHits(overlaps)]$score
df_tumor$score <- scores
df_tumor$mut <- sub(".*\\[([A-Z]>[A-Z])\\].*", "\\1", df_tumor$channel)

ggplot() +
  theme_bw()+
  geom_density(data = df_tumor, aes(score, color = mut), linewidth = 0.8)+
  geom_density(data = data.frame(score= gr$score), aes(score),
               fill = "blue", alpha = 0.1, color = "blue")+
  scale_color_manual(values = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))


df_tmp <- rbind(data.frame(score = gr$score, mut = "AllTrack"),
                df_tumor[, c("score", "mut")])

ggplot() +
  theme_bw()+
  geom_boxplot(data = df_tmp, aes(x = mut, y = score, color = mut, fill = mut), alpha = 0.6, outlier.size = 0.5)+
  scale_fill_manual(values = c("white", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  scale_color_manual(values = c("black", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  ylim(c(0, 5))


#----------------------------------------------------------------- Replication
df_tumor <- as.data.frame(gr_tumor)

gr <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_RepliTiming_1kb_standardized.bigWig.rds.gzip")
gr <- gr[gr$bin_weight > 0]

overlaps <- findOverlaps(gr_tumor, gr)
scores <- rep(0, length(gr_tumor))
scores[queryHits(overlaps)] <- gr[subjectHits(overlaps)]$score
plot(density(gr$score))
lines(density(scores), col = "red")

df_tumor$score <- scores
df_tumor$mut <- sub(".*\\[([A-Z]>[A-Z])\\].*", "\\1", df_tumor$channel)

ggplot() +
  theme_bw()+
  geom_density(data = df_tumor, aes(score, color = mut), linewidth = 0.8)+
  geom_density(data = data.frame(score= gr$score), aes(score),
               fill = "blue", alpha = 0.1, color = "blue")+
  scale_color_manual(values = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))


df_tmp <- rbind(data.frame(score = gr$score, mut = "AllTrack"),
                df_tumor[, c("score", "mut")])
ggplot() +
  theme_bw()+
  geom_boxplot(data = df_tmp, aes(x = mut, y = score, color = mut, fill = mut), alpha = 0.6, outlier.size = 0.5)+
  scale_fill_manual(values = c("white", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  scale_color_manual(values = c("black", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))


#----------------------------------------------------------------- Nucleosomes
df_tumor <- as.data.frame(gr_tumor)

gr <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_NuclOccup_1kb_standardized.bigWig.rds.gzip")
gr <- gr[gr$bin_weight > 0]

overlaps <- findOverlaps(gr_tumor, gr)
scores <- rep(0, length(gr_tumor))
scores[queryHits(overlaps)] <- gr[subjectHits(overlaps)]$score
plot(density(gr$score))
lines(density(scores), col = "red")

df_tumor$score <- scores
df_tumor$mut <- sub(".*\\[([A-Z]>[A-Z])\\].*", "\\1", df_tumor$channel)

ggplot() +
  theme_bw()+
  geom_density(data = df_tumor, aes(score, color = mut), linewidth = 0.8)+
  geom_density(data = data.frame(score= gr$score), aes(score),
               fill = "blue", alpha = 0.1, color = "blue") +
  scale_color_manual(values = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  xlim(c(0, 5))

df_tmp <- rbind(data.frame(score = gr$score, mut = "AllTrack"),
                df_tumor[, c("score", "mut")])
ggplot() +
  theme_bw()+
  geom_boxplot(data = df_tmp, aes(x = mut, y = score, color = mut, fill = mut), alpha = 0.6, outlier.size = 0.5)+
  scale_fill_manual(values = c("white", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  scale_color_manual(values = c("black", "#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5"))+
  ylim(c(0, 2.5))

palette <- c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")
names(palette) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
df_pal <- data.frame(channel = df_tumor$channel,
                     palette = palette[df_tumor$mut]) %>%
  distinct() %>% arrange(channel)
ggplot(df_tumor) +
  theme_bw()+
  geom_density(aes(score, color = channel), linewidth = 0.5)+
  geom_density(data = data.frame(score= gr$score), aes(score),
               fill = "blue", alpha = 0.1, color = "blue")+
  xlim(c(0, 5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = df_pal$palette)



#--------------------------------------------------------
# Let's look at total number of mutations in 1kb regions
# Import blacklist and regions with gaps
genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)
# Load gr_tumor
gr_tumor <- readRDS("~/SigPoisProcess/data/ICGC_snp_mutations/Stomach-AdenoCA_snp.rds.gzip")
# Remove mutations from blacklist
blacklist <- rtracklayer::import("~/SigPoisProcess/data/regions_to_exclude/hg19-blacklist.v2.bed")
overlaps <- findOverlaps(gr_tumor, blacklist)
gr_tumor <- gr_tumor[-queryHits(overlaps)]

tmp <- readRDS("~/SigPoisProcess/data/Covariates_standardized/Stomach-AdenoCA_male_H3K27me3_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(tmp, gr_tumor)
gr_tumor$region <- queryHits(hits)
df_regions <- as.data.frame(gr_tumor) %>%
  group_by(region) %>%
  summarize(n = n())

df_cov <- as.data.frame(tmp)
df_cov$region <- 1:nrow(df_cov)
df_cov <- df_cov %>%
  left_join(df_regions, by = "region") %>%
  mutate(n = case_when(is.na(n)~ 0,
                       TRUE ~ n)) %>%
  filter(bin_weight > 0)

cor(df_cov$n, df_cov$score)
plot(density((df_cov$score[df_cov$n > 0])))
lines(density(df_cov$score[df_cov$n ==0]), col = "red")


f <- "ENCFF194BNH.bed.gz"
gr <- load_bedfile(paste0(dir, f))
seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24]
gr_cov <- as(coverage(gr), "GRanges")
hits <- findOverlaps(gr_cov, gr_tumor)
gr_tumor$region <- queryHits(hits)
df_regions <- as.data.frame(gr_tumor) %>%
  group_by(region) %>%
  summarize(n = n()) %>%
  left_join(as.data.frame(gr_cov) %>%
              mutate(width = width(gr_cov),
                     region = 1:length(gr_cov)), by = "region")
df_regions$n_adj <- df_regions$n/df_regions$width

hist(df_regions$n/df_regions$width,
     breaks = 100)
mean(df_regions$n_adj[df_regions$score == 1])
mean(df_regions$n_adj[df_regions$score == 0])

#-------- We do the following
# For replication timing, we split it into cell type
# For histones and CTCF, simply tell if mutation falls into peak or not
# for GC content, split into medium, low and high
# Temporary pause on nucleosomes

df_areas <- data.frame(sample = colnames(MutMatrix))
#---------------------------------------------------------------------- Merge with histone marks and CTCF
gr_tumor_merged <- gr_tumor
dir <- "~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bedfiles/"
files <- df_otlu$`File Name(s)`
for(i in 1:length(files)){
  f <- files[i]
  gr <- load_bedfile(paste0(dir, f))
  overlaps <- findOverlaps(gr_tumor_merged, gr)
  score_all <- (gr$score - min(gr$score))/(max(gr$score) - min(gr$score))
  scores <- rep(0, length(gr_tumor_merged))
  scores[queryHits(overlaps)] <- (gr[subjectHits(overlaps)]$score - min(gr$score))/(max(gr$score) - min(gr$score))
  gr_tumor_merged$score <- scores
  colnames(mcols(gr_tumor_merged))[ncol(mcols(gr_tumor_merged))] <- paste0(df_otlu$target[i], "_", f)

  # Areas
  df_areas$area <- sum((gr$score - min(gr$score))/(max(gr$score) - min(gr$score)) * width(gr))
  colnames(df_areas)[ncol(df_areas)] <- paste0(df_otlu$target[i], "_", f)
}

counts <- colSums(as.matrix(mcols(gr_tumor_merged)[, -c(1:3)]) > 0)
types <- sub("_.*", "", names(counts))
df_tp <- data.frame(names = names(counts), counts, types) %>%
  group_by(types) %>%
  filter(counts == max(counts)) %>%
  ungroup()

mcols(gr_tumor_merged) <- mcols(gr_tumor_merged)[c("tumor", "sample", "channel", df_tp$names)]
df_areas <- df_areas[, c("sample", df_tp$names)]

#---------------------------------------------------------------------- Merge with Replitime
# Import all tracks
data_path <- "~/SigPoisProcess/data/ENCODE_pvalues/ReplTiming/Repliseq_WaveSignal_bigWig/"
files <- list.files(path = data_path)
files <- files[grepl("Nhek", files)]
files_full <- file.path(data_path, files)
import_bw <- function(x) rtracklayer::import(x, format="bigWig")
# Import all files and store as a named list
track_list <- setNames(lapply(files_full, import_bw), files)

migrate_to_1kb <- function(gr, std_chr = paste0("chr", c(1:22, "X"))){
  gr <- keepSeqlevels(gr, std_chr, pruning.mode = "coarse")
  seqlengths(gr) <- seqlengths(genome_1kb)[-24]
  numvar <- coverage(gr, weight = "score")
  gr <- binnedAverage(bins = keepSeqlevels(genome_1kb, std_chr, pruning.mode = "coarse"),
                                numvar = numvar,
                                varname = "score")
  gr
  #----- Remove the problematic ranges
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
  gr$bin_weight <- BinWeight
  gr[gr$bin_weight > 0]
}

# Tranform to 1kb
G1b <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekG1bPctSignalRep1.bigWig)
S1  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS1PctSignalRep1.bigWig)
S2  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS2PctSignalRep1.bigWig)
S3  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS3PctSignalRep1.bigWig)
S4  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekS4PctSignalRep1.bigWig)
G2  <- migrate_to_1kb(track_list$wgEncodeUwRepliSeqNhekG2PctSignalRep1.bigWig)

# Filter every region with the blacklist and regions to exclude
gr_repli <- G1b
gr_repli$G1b <- G1b$score
gr_repli$S1 <- S1$score
gr_repli$S2 <- S2$score
gr_repli$S3 <- S3$score
gr_repli$S4 <- S4$score
mcols(gr_repli)[, -c(1,2)] <- as.matrix(mcols(gr_repli)[, -c(1,2)])/100

overlaps <- findOverlaps(gr_tumor_merged, gr_repli)
MatScoresRepli <- matrix(0, nrow = length(gr_tumor_merged), ncol = 6)
MatScoresRepli[queryHits(overlaps), ] <- as.matrix(mcols(gr_repli[subjectHits(overlaps)][, -1]))
colnames(MatScoresRepli) <- c("width", "G1b", "S1", "S2", "S3", "S4")
mcols(gr_repli)[, -c(1,2)] <- as.matrix(mcols(gr_repli)[, -c(1,2)])/100
mcols(gr_tumor_merged) <- cbind(mcols(gr_tumor_merged), MatScoresRepli[, -1])
df_areas <- cbind(df_areas, t(crossprod(MatScoresRepli[, -1], MatScoresRepli[, 1])))

#---------------------------------------------------------------------- Merge with Methylation
gr <- import("~/SigPoisProcess/data/ENCODE_pvalues/old_data/Stomach_AdenoCA/methylation_bigWig_GRCh38/Stomach_methylation_male_hg19_1kb.bigWig")
filter_blacklist <- function(gr){
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
  gr$bin_weight <- BinWeight
  gr
}
gr <- filter_blacklist(gr)
# Normalize it between 0 and 1
sc <- gr[gr$bin_weight > 0]$score
sc2 <- (pmin(sc, quantile(sc, 0.999)) - min(sc))/(quantile(sc, 0.999) -  min(sc))
gr[gr$bin_weight > 0]$score <- sc2
gr[gr$bin_weight == 0]$score <- 0
overlaps <- findOverlaps(gr_tumor_merged, gr)
gr_tumor_merged$Methyl <- gr[subjectHits(overlaps)]$score
df_areas$Methyl <- sum(gr$score * gr$bin_weight)

#---------------------------------------------------------------------- Merge with GC content
gr <- readRDS("~/SigPoisProcess/data/Covariates_standardized/male_GC_1kb_standardized.bigWig.rds.gzip")
overlaps <- findOverlaps(gr_tumor_merged, gr)
gr_tumor_merged$GC <- gr[subjectHits(overlaps)]$score
df_areas$GC <- sum(gr$score * gr$bin_weight)


K <- ncol(resMAP$Signatures)
J <- ncol(resMAP$Theta)
p <- ncol(mcols(gr_tumor_merged)[-c(1,2,3)])
R_start <- as.matrix(resMAP$Signatures)
Xtotals <- as.matrix(df_areas[, -c(1)])
Betas_start <- array(NA, dim = c(K, J, p))
for(i in 1:p){
  Betas_start[, , i] <- t(apply(resMAP$Theta, 1, function(x) x/Xtotals[, i]))
}
Mu_start <- matrix(rowMeans(resMAP$Theta))[, rep(1, p)]


# MLE solution
set.seed(10)
out <- SigPoisProcess(gr_Mutations = gr_tumor_merged,
                      df_areas = df_areas,
                      K = K,
                      method = "mle",
                      init =  SigPoisProcess.init(R_start = R_start,
                                                  Betas_start = Betas_start,
                                                  Mu_start = Mu_start),
                      controls = SigPoisProcess.controls(maxiter = 1000,
                                                         merge_move = FALSE))
plot(out)

#---------------------------------------------------------------------- Add a baseline
gr_tumor_merged$baseline <- 1
df_areas$baseline <- sum(1 * gr$bin_weight)

K <- ncol(resMAP$Signatures)
J <- ncol(resMAP$Theta)
p <- ncol(mcols(gr_tumor_merged)[-c(1,2,3)])
R_start <- as.matrix(resMAP$Signatures)
Xtotals <- as.matrix(df_areas[, -c(1)])
Betas_start <- array(NA, dim = c(K, J, p))
for(i in 1:p){
  Betas_start[, , i] <- t(apply(resMAP$Theta, 1, function(x) x/Xtotals[, i]))
}
Mu_start <- matrix(rowMeans(resMAP$Theta))[, rep(1, p)]


# MLE solution
set.seed(10)
out_base <- SigPoisProcess(gr_Mutations = gr_tumor_merged,
                      df_areas = df_areas,
                      K = 0,
                      method = "map",
                      #init =  SigPoisProcess.init(R_start = R_start,
                      #                            Betas_start = Betas_start,
                      #                            Mu_start = Mu_start),
                      controls = SigPoisProcess.controls(maxiter = 300,
                                                         merge_move = FALSE),
                      prior_params = SigPoisProcess_mult.PriorParams(SigPrior = 1e6 * R_start + 1))
plot(out_base)













