# This file takes the table in the Otlu paper and extracts the name
# of all the bigWig files associated to the bw
library(tidyverse)
library(httr)
library(jsonlite)

#----
get_encode_metadata <- function(accession_id) {
  url <- paste0("https://www.encodeproject.org/files/", accession_id, "/?format=json")
  res <- GET(url, user_agent("encode-metadata-query"))
  if (status_code(res) == 200) {
    content(res, as = "parsed", simplifyVector = TRUE)
  } else {
    warning(paste("Failed to fetch:", accession_id))
    NULL
  }
}

get_experiment_files_metadata <- function(experiment_id) {
  # Query the ENCODE API
  url <- paste0("https://www.encodeproject.org/experiments/", experiment_id, "/?format=json")
  res <- GET(url, user_agent("encode-metadata-query"))

  if (status_code(res) != 200) {
    warning(paste("Failed to fetch metadata for:", experiment_id))
    return(NULL)
  }

  # Parse JSON response
  meta <- content(res, as = "parsed", simplifyVector = TRUE)

  # Extract experiment-level metadata
  exp_description <- ifelse(is.null(meta$description), NA, meta$description)
  target_label <- ifelse(is.null(meta$target$label), NA, meta$target$label)
  classification <- ifelse(is.null(meta$biosample_ontology$classification), NA, meta$biosample_ontology$classification)
  term_name <- ifelse(is.null(meta$biosample_ontology$term_name), NA, meta$biosample_ontology$term_name)
  biosample_summary <- ifelse(is.null(meta$biosample_summary), NA, meta$biosample_summary)

  # Loop over each file associated with the experiment
  # If meta$files is a dataframe
  if (!is.data.frame(meta$files)) {
    warning("meta$files is not a data frame for: ", experiment_id)
    return(NULL)
  }

  # Minimal columns to extract from meta$files
  file_meta <- meta$files

  # Keep only useful columns, and add experiment-level annotation
  df <- data.frame(
    accession_experiment = experiment_id,
    description = exp_description,
    target = target_label,
    classification = classification,
    term_name = term_name,
    biosample_summary = biosample_summary,
    accession = file_meta$accession,
    file_format = file_meta$file_format,
    file_type = file_meta$file_type,
    preferred_default = file_meta$preferred_default,
    output_type = file_meta$output_type,
    assay_title = file_meta$assay_title,
    assembly = file_meta$assembly,
    stringsAsFactors = FALSE
  )

  # Combine into one dataframe
  df
}

#----

# Let's download all files for a selected group of cancer types
df_files_Cell2023 <- readxl::read_xlsx(path = "data/ENCODE_pvalues/Otlu_2013_listfiles.xlsx", skip = 3)

df_files <- df_files_Cell2023 %>%
  group_by(`Cancer Type`, `Topography Feature`, `Data Source`, `File Format(s)`) %>%
  separate_rows(`File Name(s)`, sep = ", ") %>%
  mutate(`File Name(s)` = str_trim(`File Name(s)`)) %>%
  filter(`Data Source` == "ENCODE",
         `File Format(s)` == "bed",
         `Topography Feature` %in% c("CTCF Binding", "Histone Modifications")) %>%
  as.data.frame()

accession_ids <- gsub(".bed.gz", "", df_files$`File Name(s)`)


# Fetch metadata for all accessions
df_all_encode <- data.frame()
for(i in 1:length(accession_ids)){
  print(i)
  x <- get_encode_metadata(accession_ids[i])
  experiment <- gsub("/", "", gsub("/experiments/", "", x$dataset))
  experiment_metadata <- get_experiment_files_metadata(experiment)
  experiment_metadata <- cbind(df_files[rep(i, nrow(experiment_metadata)), ], experiment_metadata)
  rownames(experiment_metadata) <- NULL
  df_all_encode <- rbind(df_all_encode, experiment_metadata)
}

write_tsv(df_all_encode, "data/ENCODE_pvalues/Otlu_merged_with_experiments.tsv")

#--------------------------------------------------

# Filter dataset based on bigWig, preferred default and hg19 file.
# We only focus on 5 main tumors
df_files_to_use <- df_all_encode %>%
  filter(file_type == "bigWig",
         assembly == "hg19",
         preferred_default,
         `Cancer Type` %in% c("Skin-Melanoma", "Liver-HCC", "Breast-Cancer",
                              "Stomach-AdenoCA", "Prost-AdenoCA"),
         target %in% c("CTCF", "H3K27ac", "H3K27me3", "H3K36me3",
                       "H3K4me1", "H3K4me3", "H3K9me3"))

write_tsv(df_files_to_use, "data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bigWigfiles_to_use.tsv")

# Create full ENCODE download URLs
urls <- sprintf(
  "https://www.encodeproject.org/files/%s/@@download/%s.bigWig",
  df_files_to_use$accession, df_files_to_use$accession
)
# Write all URLs into a text file
writeLines(urls, "data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bigWig/files_Stomach_Skin_Liver_Breast_Prost.txt")
# On the terminal:
# nohup xargs -n 1 curl -O -L < files_Stomach_Skin_Liver_Breast_Prost.txt > download.log 2>&1 &







