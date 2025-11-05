library(tidyverse)
df_all_encode <- read_tsv("~/SigPoisProcess/data/ENCODE_pvalues/Otlu_merged_with_experiments.tsv")

# Accession files for breast cancer (epithelial cell and breast epithelium)
accessions <- c("ENCFF143AFG", "ENCFF663JVM","ENCFF216DKX",
                "ENCFF984VUL", "ENCFF529OVO", "ENCFF764JHZ",
                "ENCFF467DBD", "ENCFF353WUO", "ENCFF353WUO",
                "ENCFF101RHP", "ENCFF709GNN", "ENCFF322JKF",
                "ENCFF981WTU", "ENCFF531OEA", "ENCFF094MQS",
                "ENCFF240TPI", "ENCFF925FUX", "ENCFF175PNA")

df_files_to_use <- df_all_encode %>%
  filter(accession %in% accessions)
write_tsv(df_files_to_use, "~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/Files_to_use_breast.tsv")

urls <- sprintf(
  "https://www.encodeproject.org/files/%s/@@download/%s.bigWig",
  df_files_to_use$accession, df_files_to_use$accession
)
writeLines(urls, "~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/files_breast.txt")
# On the terminal:
# nohup xargs -n 1 curl -O -L < files_breast.txt > download.log 2>&1 &


# table(df_files_to_use$target, df_files_to_use$term_name)
# df_files_to_use %>%
#   filter(target == "H3K9me3") %>%
#   as.data.frame()
#
# region <- GRanges("chr1", IRanges(start = 1000000, end = 2000000))
# track1 <- import("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/ENCFF764JHZ.bigWig", which = region)
# track2 <- import("~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/ENCFF709GNN.bigWig", which = region)
#
# plot(start(track1), track1$score, type = "l")
# lines(start(track2), track2$score, col = "red")
#
#
