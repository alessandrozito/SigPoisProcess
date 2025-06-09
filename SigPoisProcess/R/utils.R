get_df_trinucl <- function(maf, genome = BSgenome.Hsapiens.UCSC.hg19){

  # Baseline mutational channels
  nucleotides <- c("A", "C", "G", "T")
  substitutions <- c("C>A", "C>T", "C>G",
                     "T>A", "T>G", "T>C")
  channels_96 <- apply(expand.grid(nucleotides, substitutions[1:6], nucleotides),
                       1, function(x) paste0(x[1], "[", x[2], "]", x[3]))
  channels_96 <- sort(channels_96)

  # Obtain mutational channels from data
  mutation_data <- maf@data %>%
    dplyr::filter(Variant_Type == "SNP")
  # Create the trinucleotide contect
  data <- GRanges(sample = mutation_data$Tumor_Sample_Barcode,
                  mutation_data$Chromosome,
                  pos = mutation_data$Start_Position,
                  IRanges(start = mutation_data$Start_Position - 1, width = 3),
                  ref = DNAStringSet(mutation_data$Reference_Allele),
                  alt = DNAStringSet(mutation_data$Tumor_Seq_Allele2))

  data$context <- Biostrings::getSeq(genome, data)

  # Obtain the context and the reverse context
  data$cref <- complement(data$ref)
  data$calt <- complement(data$alt)
  data$rccontext <- reverseComplement(data$context)

  # Extract the final mutation category (mutational channel)
  data$channel <- ifelse(as.character(data$ref) %in% c("C", "T"),
                         paste0(subseq(data$context, 1, 1), "[", data$ref, ">", data$alt, "]", subseq(data$context, 3, 3)),
                         paste0(subseq(data$rccontext,1, 1), "[", data$cref, ">", data$calt, "]", subseq(data$rccontext, 3, 3)))
  data$channel <- factor(data$channel, levels = channels_96)

  df_trinucl <- as.data.frame(data) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Chromosome = seqnames,
                  sample = as.factor(sample)) %>%
    dplyr::select(sample, Chromosome, start, end, channel) %>%
    dplyr::arrange(sample, Chromosome, start)

  return(df_trinucl)
}

#' @export
match_to_RefSigs <- function(sigs, ref){
  all <- sigminer::cosine(sigs, ref)
  ids <- cbind(1:nrow(all), apply(all, 1, which.max))
  data.frame("sigs" = colnames(sigs), "best" = colnames(ref)[ids[, 2]], "cosine" = round(all[ids], 3))
}

#' @export
estimateTotalCounts.SigPoisProcess <- function(object, ...){
  R <- object$Signatures
  Betas <- object$Betas
  X_total <- object$Xbar[1, , ]
  MatOut <-sapply(1:ncol(Betas), function(j) R %*% Betas[, j, ] %*% X_total[j, ], simplify = TRUE)
  colnames(MatOut) <- colnames(Betas)
  rownames(MatOut) <- rownames(R)
  return(MatOut)
}

#' @export
getTotalMutations <- function(gr_Mutations){
  Mutations <- as.data.frame(mcols(gr_Mutations)) %>%
    dplur::select(sample, channel)
  Mutations$sample <- factor(Mutations$sample, levels = rownames(X_totals))
  Mutations$channel <- as.factor(Mutations$channel)
  X_all <- table(Mutations$channel, Mutations$sample)
  matrix(X_all, nrow = nrow(X_all), dimnames = dimnames(X_all))
}







