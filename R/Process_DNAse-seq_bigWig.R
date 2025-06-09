# Read all DNAse files
library(tidyverse)
library(httr)
library(jsonlite)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

#--------------------------------------------------------------
# Chain file URL
chain_url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz"

# Download and decompress
chain_file <- tempfile(fileext = ".chain.gz")
download.file(chain_url, chain_file)
chain_unzipped <- R.utils::gunzip(chain_file, remove = FALSE, temporary = TRUE)

# Import the chain
chain <- import.chain(chain_unzipped)
#--------------------------------------------------------------
load_DNAse_files <- FALSE
if(load_DNAse_files){


  df_DNAse <- data.frame()

  #------- Stomach
  df_stomach <- read_tsv("data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/stomach_DNAse.tsv", skip = 1) %>%
  filter(`Preferred Default`,
         !is.na(`Simple biosample summary`),
         grepl("adult", `Simple biosample summary`)) %>%
  mutate(gender = if_else(str_detect(`Simple biosample summary`, "female"), "female", "male")) %>%
  as.data.frame()

  df_DNAse <- rbind(df_DNAse,
                  data.frame(accession = df_stomach$Accession,
                             assembly = df_stomach$`Genome assembly`,
                             assay = df_stomach$`Assay term name`,
                             organ = df_stomach$`Biosample name`,
                             biosample_summary = df_stomach$`Simple biosample summary`))

  #------- Liver
  df_liver <- read_tsv("data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/liver_DNAse.tsv")
  df_liver <- df_liver %>%
  filter(`File format` == "bigWig") %>%
  dplyr::select(`File accession`, `File assembly`, `Experiment accession`, Assay,
                `Biosample term name`, Project) %>%
  as.data.frame()
  # Find biosample summary by experiment
  url <- paste0("https://www.encodeproject.org/experiments/", df_liver$`Experiment accession`, "/?format=json")
  df_liver$`Simple biosample summary` <- sapply(url, function(x) {
  res <- GET(x, user_agent("encode-metadata-query"))
  meta <- content(res, as = "parsed", simplifyVector = TRUE)
  meta$biosample_summary})

  df_DNAse <- rbind(df_DNAse,
                  data.frame(accession = df_liver$`File accession`,
                             assembly = df_liver$`File assembly`,
                             assay = df_liver$Assay,
                             organ = df_liver$`Biosample term name`,
                             biosample_summary = df_liver$`Simple biosample summary`))

  #------- Breast
  df_breast <- read_tsv("data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/breast_DNAse.tsv")
  df_breast <- df_breast %>%
  filter(`File format` == "bigWig") %>%
  dplyr::select(`File accession`, `File assembly`, `Experiment accession`, Assay,
                `Biosample term name`, Project) %>%
  as.data.frame()
  # Find biosample summary by experiment
  url <- paste0("https://www.encodeproject.org/experiments/", df_breast$`Experiment accession`, "/?format=json")
  df_breast$`Simple biosample summary` <- sapply(url, function(x) {
  res <- GET(x, user_agent("encode-metadata-query"))
  meta <- content(res, as = "parsed", simplifyVector = TRUE)
  meta$biosample_summary})

  df_DNAse <- rbind(df_DNAse,
                  data.frame(accession = df_breast$`File accession`,
                             assembly = df_breast$`File assembly`,
                             assay = df_breast$Assay,
                             organ = df_breast$`Biosample term name`,
                             biosample_summary = df_breast$`Simple biosample summary`))


  #------- Prostate
  df_prost <- read_tsv("data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/prostate_DNAse.tsv")
  df_prost <- df_prost %>%
  filter(`File format` == "bigWig") %>%
  dplyr::select(`File accession`, `File assembly`, `Experiment accession`, Assay,
                `Biosample term name`, Project) %>%
  as.data.frame()
  # Find biosample summary by experiment
  url <- paste0("https://www.encodeproject.org/experiments/", df_prost$`Experiment accession`, "/?format=json")
  df_prost$`Simple biosample summary` <- sapply(url, function(x) {
  res <- GET(x, user_agent("encode-metadata-query"))
  meta <- content(res, as = "parsed", simplifyVector = TRUE)
  meta$biosample_summary})

  df_DNAse <- rbind(df_DNAse,
                  data.frame(accession = df_prost$`File accession`,
                             assembly = df_prost$`File assembly`,
                             assay = df_prost$Assay,
                             organ = df_prost$`Biosample term name`,
                             biosample_summary = df_prost$`Simple biosample summary`))



  #------- skin
  df_skin <- read_tsv("data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/skin_DNAse.tsv")
  df_skin <- df_skin %>%
  filter(`File format` == "bigWig") %>%
  dplyr::select(`File accession`, `File assembly`, `Experiment accession`, Assay,
                `Biosample term name`, Project) %>%
  as.data.frame()
  # Find biosample summary by experiment
  url <- paste0("https://www.encodeproject.org/experiments/", df_skin$`Experiment accession`, "/?format=json")
  df_skin$`Simple biosample summary` <- sapply(url, function(x) {
  res <- GET(x, user_agent("encode-metadata-query"))
  meta <- content(res, as = "parsed", simplifyVector = TRUE)
  meta$biosample_summary})

  df_DNAse <- rbind(df_DNAse,
                  data.frame(accession = df_skin$`File accession`,
                             assembly = df_skin$`File assembly`,
                             assay = df_skin$Assay,
                             organ = df_skin$`Biosample term name`,
                             biosample_summary = df_skin$`Simple biosample summary`))

  write_tsv(x = df_DNAse, file = "data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/DNAse-seq_BigWig.tsv")

  # File to download the
  urls <- sprintf(
  "https://www.encodeproject.org/files/%s/@@download/%s.bigWig",
  df_DNAse$accession, df_DNAse$accession
  )
  writeLines(urls, "data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/files_DNAse_seq_Stomach_Skin_Liver_Breast_Prost.txt")

  # On the terminal
  # nohup xargs -n 1 curl -O -L < files_DNAse_seq_Stomach_Skin_Liver_Breast_Prost.txt > download.log 2>&1 &

}
# Migraete the files to hg19 and shrink them to 1kb
df_DNAse <- read_tsv("data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/DNAse-seq_BigWig.tsv")
names_bigWig <- list.files("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/")
names_bigWig <- names_bigWig[grepl("[.]bigWig", names_bigWig)]
#df_DNAse$accession %in% gsub(".bigWig","", names_bigWig)

# Group into 1kb
# Load the genome and split it
genome_1kb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24],
                         tilewidth = 1000,
                         cut.last.tile.in.chrom = TRUE)
std_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Files specific to cancers
for(i in 1:length(names_bigWig)){
  print(i)
  gr <- rtracklayer::import(paste0("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_bigWig/", names_bigWig[i]))
  gr <- keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  nm <- gsub(".bigWig","", names_bigWig[i])
  assembly <- df_DNAse[df_DNAse$accession == nm, ]$assembly
  if(assembly == "hg19"){
    # Simply split into regions
    gr <- coverage(gr, weight = "score")
  } else {
    # Migrate it into hg19 with liftOver
    genome(gr) <- "hg38"
    gr <- keepSeqlevels(unlist(liftOver(gr, chain)), std_chrs, pruning.mode = "coarse")
    genome(gr) <- "hg19"; seqlengths(gr) <- seqlengths(genome_1kb)
    gr <- coverage(gr, weight = "score")
  }
  binned_score <- binnedAverage(bins = genome_1kb,
                                numvar = gr,
                                varname = "score")
  rm(gr)
  export.bw(binned_score,
            con = paste0("~/SigPoisProcess/data/ENCODE_pvalues/ENCODE_rawdata/DNAse-seq_1kb/", names_bigWig[i]))
}



