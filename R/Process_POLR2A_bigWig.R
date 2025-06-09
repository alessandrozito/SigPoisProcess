# This file dowloads the data from POLR2A and aggregates them under the usual
# 1kb window.
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

# Load the data
df_polr2a <- read_tsv("data/ENCODE_pvalues/ENCODE_rawdata/POLR2A_bigWig/POLR2A_accession_files.tsv", skip = 1) %>%
  as.data.frame()

df_polr2a <- df_polr2a %>%
  filter(`Preferred Default`,
         !is.na(`Simple biosample summary`))

# Write all URLs into a text file
urls <- sprintf(
  "https://www.encodeproject.org/files/%s/@@download/%s.bigWig",
  df_polr2a$Accession, df_polr2a$Accession
)
writeLines(urls, "data/ENCODE_pvalues/ENCODE_rawdata/POLR2A_bigWig/files_polr2a_Stomach_Skin_Liver_Breast_Prost.txt")

# On the terminal
# nohup xargs -n 1 curl -O -L < files_polr2a_Stomach_Skin_Liver_Breast_Prost.txt > download.log 2>&1 &

df_polr2a <- df_polr2a %>% dplyr::select(Organ, `Simple biosample summary`)%>%
  # Separate the bio sample summary into parts
  mutate(
    organ_type = case_when(grepl("stomach", Organ)  ~ "stomach",
                           grepl("prostate", Organ)  ~ "prostate",
                           grepl("breast", Organ)  ~ "breast",
                           grepl("skin", Organ)  ~ "skin",
                           grepl("liver", Organ)  ~ "liver"),
    gender = if_else(str_detect(`Simple biosample summary`, "female"), "female", "male"))

table(df_polr2a$organ_type, df_polr2a$gender)

# Group into 1kb
# Load the genome and split it
genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Load the files
names_bigWig <- list.files("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/POLR2A_bigWig/")
names_bigWig <- names_bigWig[grepl("[.]bigWig", names_bigWig)]

# Files specific to cancers
for(i in 1:length(names_bigWig)){
  print(i)
  gr <- rtracklayer::import(paste0("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/POLR2A_bigWig/", names_bigWig[i]))
  gr <- keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  #numvar <- coverage(gr, weight = "score")
  binned_score <- binnedAverage(bins = genome_1kb,
                                numvar = coverage(gr, weight = "score"),
                                varname = "score")
  #rm(numvar)
  rm(gr)
  export.bw(binned_score,
            con = paste0("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/POLR2A_1kb/", names_bigWig[i]))
}






