# This file extracts the gene expression for every patient
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

# Load the Gene expression matrix for each patient, as well as the locations of
# the genes

#Mexpr <- read_tsv("https://pcawg-hub.s3.us-east-1.amazonaws.com/download/tophat_star_fpkm_uq.v2_aliquot_gl.donor.log")
#gene_coords <- read_tsv("https://pcawg-hub.s3.us-east-1.amazonaws.com/download/gencode.v19.annotation.gene.probemap")
#write_tsv(Mexpr, "data/ICGC_rawdata/gene_expression_icgc.tsv")
#write_tsv(gene_coords, "data/ICGC_rawdata/gene_names_positions_icgc.tsv")

Mexpr <- read_tsv("data/ICGC_rawdata/gene_expression_icgc.tsv") %>% as.data.frame()
gene_coords <- read_tsv("data/ICGC_rawdata/gene_names_positions_icgc.tsv")

# Genome split
genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Genes split
gr_genes <- GRanges(seqnames = gene_coords$chrom,
                    ranges   = IRanges(gene_coords$chromStart, gene_coords$chromEnd),
                    strand = gene_coords$strand)
gr_genes$gene <- gene_coords$`#id`
gr_genes$strand_orig <- gene_coords$strand
gr_genes <- keepSeqlevels(gr_genes, std_chrs, pruning.mode = "coarse")
strand(gr_genes) <- "*"

# Gene expression  by patient
patients <- colnames(Mexpr)[-1]

# Let's do it only for our tumors
df_clinical <- read_tsv("data/ICGC_rawdata/pcawg_specimen_histology_August2016_v9_donor.tsv")
pats_interested <- df_clinical$icgc_specimen_id[df_clinical$histology_abbreviation %in%
                            c("Stomach-AdenoCA", "Skin-Melanoma",
                              "Prost-AdenoCA", "Breast-AdenoCA", "Liver-HCC")]
patients <- patients[patients %in% pats_interested]

for(j in 1:length(patients)){
  print(j)
  # Append score
  scores <- rep(0, length(gr_genes))
  names(scores) <- gr_genes$gene
  scores_tmp <- c(Mexpr[, patients[j]])
  names(scores_tmp) <- Mexpr$sample
  scores[names(scores_tmp)] <- 2^unname(scores_tmp) + 0.001

  # Merge score with genes in gr
  gr_tmp <- gr_genes
  gr_tmp$score <- scores[gr_tmp$gene]

  # Calculate average score per bin
  numvar <- coverage(gr_tmp, weight = "score")
  binned_score <- binnedAverage(bins = genome_1kb,
                                numvar = numvar,
                                varname = "score")

  binned_score$score[binned_score$score < 0] <- 0
  binned_score$sample <- patients[j]
  export.bw(binned_score,
            con = paste0("~/SigPoisProcess/data/ICGC_GeneExpr_bigWig/", patients[j], ".bigWig"))
}















