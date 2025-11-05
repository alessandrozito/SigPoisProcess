library(tidyverse)
df_all_encode <- read_tsv("~/SigPoisProcess/data/ENCODE_pvalues/Otlu_merged_with_experiments.tsv")

df_files_to_use <- df_all_encode %>%
  filter(file_type == "bigWig",
         assembly == "hg19",
         output_type == "fold change over control",
         `Cancer Type` %in% c("Stomach-AdenoCA")) %>%
         #target %in% c("CTCF", "H3K27ac", "H3K27me3", "H3K36me3",
        #               "H3K4me1", "H3K4me3", "H3K9me3"))
  filter(grepl("male adult", biosample_summary) | grepl("female adult", biosample_summary)) %>%
  filter(grepl("stomach", term_name) | target == "H3K27ac") %>%
  group_by(`File Name(s)`) %>%
  mutate(n = row_number()) %>%
  filter(n == 1)
write_tsv(df_files_to_use, "~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/Files_to_use.tsv")

urls <- sprintf(
  "https://www.encodeproject.org/files/%s/@@download/%s.bigWig",
  df_files_to_use$accession, df_files_to_use$accession
)
writeLines(urls, "~/SigPoisProcess/data/ENCODE_foldchange/ENCODE_rawdata/CTCF_Histones_raw/files_Stomach.txt")
# On the terminal:
# nohup xargs -n 1 curl -O -L < files_Stomach.txt > download.log 2>&1 &

