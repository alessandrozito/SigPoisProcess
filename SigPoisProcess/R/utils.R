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
estimateTotalCounts.SigPoisProcess_mult <- function(object, SignalTrack, bin_weight, ...){
  R <- object$Signatures
  Betas <- object$Betas
  Theta <- object$Thetas
  #CumSums <- apply(exp(SignalTrack %*% Betas), 2, function(x) cumsum(bin_weight * x))
  sums <- bin_weight %*% exp(SignalTrack %*% Betas)
  MatOut <- R %*% (Theta * c(sums))
  colnames(MatOut) <- colnames(Theta)
  rownames(MatOut) <- rownames(R)
  return(MatOut)
}

#' @export
estimateMutationRate.SigPoisProcess_mult <- function(i, j, R, Betas, Theta, SignalTrack, bin_weight, ...){
  CumSums <- apply(exp(SignalTrack %*% Betas), 2, function(x) cumsum(bin_weight * x))
  track <- c(crossprod(R[i, ], Theta[, j] * t(CumSums)))
  return(track)
}


#' @export
getTotalMutations <- function(gr_Mutations, df_areas = NULL){
  Mutations <- as.data.frame(GenomicRanges::mcols(gr_Mutations)) %>%
    dplyr::select(sample, channel)
  if(!is.null(df_areas)){
    Mutations$sample <- factor(Mutations$sample, levels = df_areas$sample)
  }
  #Mutations$channel <- as.factor(Mutations$channel)
  X_all <- table(Mutations$channel, Mutations$sample)
  matrix(X_all, nrow = nrow(X_all), dimnames = dimnames(X_all))
}

#' @export
compute_SignatureCovariate_probs <- function(object, gr_Mutations) {

  Mutations <- as.data.frame(GenomicRanges::mcols(gr_Mutations)) %>%
    dplyr::select(sample, channel)
  Mutations$sample <- factor(Mutations$sample, levels = colnames(object$Betas))
  Mutations$channel <- as.factor(Mutations$channel)
  channel_id <- as.numeric(Mutations$channel)
  sample_id <- as.numeric(Mutations$sample)

  Sigs <- object$Signatures
  Loads <- object$Betas
  X <- as.data.frame(GenomicRanges::mcols(gr_Mutations)) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()
  Sig_Cov_prob <- calculate_Sig_covariate_prob(X, Sigs, Loads, channel_id - 1, sample_id - 1)
  Probs <- Sig_Cov_prob$Probs

  dimnames(Probs)[[1]] <- paste0(gr_Mutations$sample, "-",  as.data.frame(granges(gr_Mutations))$seqnames,
                                 ":", start(gr_Mutations), "-", gr_Mutations$channel)
  dimnames(Probs)[[2]] <- colnames(Sigs)
  dimnames(Probs)[[3]] <- dimnames(Loads)[[3]]

  return(Probs)
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


#' @export
plot_betas <- function(Betas_sol){
  df_betas <- as.data.frame(Betas_sol) %>%
    rownames_to_column("Covariate") %>%
    gather(key = "Signature", value = "Beta", -Covariate) %>%
    dplyr::mutate(
      Covariate = factor(Covariate, levels=unique(Covariate)),
      Signature = as.factor(Signature),
      AltCol = ifelse((as.numeric(Covariate) %% 2) == 0, "gray93", "gray97"))

  plot_out <- ggplot2::ggplot(df_betas, aes(x = Covariate, y = Signature)) +
    ggplot2::geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
    ggplot2::scale_fill_manual(
      values = c("gray93" = "gray93", "gray97" = "gray97"),
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(aes(fill = Beta, size = abs(Beta)),
                        shape = 21, color = "grey30", stroke = 0.2,
                        na.rm = TRUE) +
    #ggplot2::scale_fill_gradientn(name = "Beta",
    #                              colours =c("blue","lightblue", "white", "red")) +
    scale_fill_gradient2(mid = "white", low = "blue", high = "red", midpoint = 0)+
    ggplot2::scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      axis.title = ggplot2::element_blank(),
      #axis.text.y = element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      axis.ticks.x = ggplot2::element_blank())+
    guides(size = "none")
  return(plot_out)
}

#' @export
plot_betasCI <- function(Betas_sol, lowCI, highCI){
  df_betas <- as.data.frame(Betas_sol) %>%
    rownames_to_column("Covariate") %>%
    gather(key = "Signature", value = "Beta", -Covariate)

  # Add CI columns
  df_betas$lowCI <- as.numeric(as.data.frame(lowCI) %>%
                                 gather(key = "Signature", value = "lowCI") %>%
                                 pull(lowCI))
  df_betas$highCI <- as.numeric(as.data.frame(highCI) %>%
                                  gather(key = "Signature", value = "highCI") %>%
                                  pull(highCI))


  # Add alternating color
  df_betas <- df_betas %>%
    mutate(
      Covariate = factor(Covariate, levels=unique(Covariate)),
      Signature = as.factor(Signature),
      AltCol = ifelse((as.numeric(Covariate) %% 2) == 0, "gray93", "gray97"),
      MarkZero = (lowCI <= 0 & highCI >= 0)
    )


  plot_out <- ggplot2::ggplot(df_betas, aes(x = Covariate, y = Signature)) +
    ggplot2::geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
    ggplot2::scale_fill_manual(
      values = c("gray93" = "gray93", "gray97" = "gray97"),
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(aes(fill = Beta, size = abs(Beta)),
                        shape = 21, color = "grey30", stroke = 0.2,
                        na.rm = TRUE) +
    geom_point(data = df_betas %>% filter(MarkZero),
               aes(x = Covariate, y = Signature), shape = 4, size = 2, color = "gray30") +
    #ggplot2::scale_fill_gradientn(name = "Beta",
    #                              colours =c("blue","lightblue", "white", "red")) +
    scale_fill_gradient2(mid = "white", low = "blue", high = "red", midpoint = 0)+
    ggplot2::scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0),
      axis.title = ggplot2::element_blank(),
      #axis.text.y = element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      axis.ticks.x = ggplot2::element_blank())+
    guides(size = "none")
  return(plot_out)
}

get_PosteriorCI <- function(chain){
  chainMean <- apply(chain, c(2,3), mean)
  chainHighCI <- resMean <- apply(chain, c(2,3), function(x) quantile(x, 0.975))
  chainLowCI <- apply(chain, c(2,3), function(x) quantile(x, 0.025))
  return(list("mean" = chainMean, "lowCI" = chainLowCI, "highCI"= chainHighCI))
}


