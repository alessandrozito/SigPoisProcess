get_df_trinucl <- function(maf, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){

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

#' Function to match the signatures with a reference set. Defaults to cosmic
#' @param sigs Matrix of estimated signatures
#' @param ref  Matrix of reference signatures
#' @export
match_to_RefSigs <- function(sigs, ref = COSMIC_v3.4_SBS96_GRCh37){
  all <- sigminer::cosine(sigs, ref)
  ids <- cbind(1:nrow(all), apply(all, 1, which.max))
  data.frame("sigs" = colnames(sigs), "best" = colnames(ref)[ids[, 2]], "cosine" = round(all[ids], 3))
}


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Posterior credible intervals for the quantities
#' @param chain  Array containg the MCMC samples.
#' @export
get_PosteriorCI <- function(chain){
  chainMean <- apply(chain, c(2,3), mean)
  chainHighCI <- resMean <- apply(chain, c(2,3), function(x) quantile(x, 0.975))
  chainLowCI <- apply(chain, c(2,3), function(x) quantile(x, 0.025))
  return(list("mean" = chainMean, "lowCI" = chainLowCI, "highCI"= chainHighCI))
}



#' SBS96 informative prior parameters for GRCh37 v3.4
#'
#' Beta parameters for 96 single base substitution (SBS) mutational signatures
#' from COSMIC v3.4 using the GRCh37 genome reference.
#'
#' @format A vector of length 96
"Betah_SBS96_GRCh37_v3.4"

#' COSMIC v3.4 SBS96 Mutational Signatures for GRCh37
#'
#' SBS96 signatures from the COSMIC database, v3.4, for GRCh37.
#'
#' @format A data frame with 96 rows and 86 columns.
#' @source \url{https://cancer.sanger.ac.uk/signatures/downloads/}
"COSMIC_v3.4_SBS96_GRCh37"
