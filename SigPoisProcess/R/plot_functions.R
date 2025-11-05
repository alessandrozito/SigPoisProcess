#' Plot the regression coefficient matrix from the Poisson process factorization
#'
#' @param Betas Betas matrix
#' @param lowCI Credible interval for the Betas matrix
#' @param highCI Credible interval for the Betas matrix
#'
#' @export
#' @import patchwork
#' @import ggplot2
plot_Betas <- function(Betas, lowCI = NULL, highCI = NULL) {

  df_betas <- as.data.frame(Betas) %>%
    tibble::rownames_to_column("Covariate") %>%
    tidyr::gather(key = "Signature", value = "Beta", -Covariate)

  if(!(is.null(lowCI) & is.null(highCI))) {
    df_low <- as.data.frame(lowCI) %>%
      tibble::rownames_to_column("Covariate") %>%
      tidyr::gather(key = "Signature", value = "lowCI", -Covariate)

    df_high <- as.data.frame(highCI) %>%
      tibble::rownames_to_column("Covariate") %>%
      tidyr::gather(key = "Signature", value = "highCI", -Covariate)

    df_betas <- df_betas %>%
      dplyr::left_join(df_low, by = c("Covariate", "Signature")) %>%
      dplyr::left_join(df_high, by = c("Covariate", "Signature")) %>%
      dplyr::mutate(Covariate = factor(Covariate, levels = unique(Covariate)),
             Signature = as.factor(Signature),
             Beta = ifelse(lowCI < 0 & highCI > 0, NA, Beta))
  } else {
    df_betas <- df_betas %>%
      mutate(Covariate = factor(Covariate, levels = unique(Covariate)),
             Signature = as.factor(Signature))
  }


  max_abs_beta <- max(abs(df_betas$Beta), na.rm = TRUE)

  oldopt <- options()
  options(ggplot2.continuous.colour = "ignore",
          ggplot2.continuous.fill = "ignore")

  plot_out <- ggplot2::ggplot(df_betas, aes(x = Covariate, y = Signature, fill = Beta)) +
    ggplot2::geom_tile(color = "black", width = 1, height = 1, na.rm = TRUE) +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"),
      limits = c(-max_abs_beta, max_abs_beta),
      na.value = "gray90"
    ) +
    ggplot2::geom_text(aes(label = ifelse(is.na(Beta), "", round(Beta, 2))), size = 3, na.rm = TRUE) +
    ggplot2::scale_y_discrete(limits = rev(levels(df_betas$Signature))) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = -1, size = 10,
                                 colour = "black"),
      axis.title = ggplot2::element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank())

  options(oldopt)

  return(plot_out)
}

#' Plot the matrix of mutational signatures
#'
#' @param signatures Matrix of signature to plot
#' @param lowCI Matrix containing the higher end of the posterior credible interval
#' @param highCI Matrix containing the higher end of the posterior credible interval
#' @param palette Palette for the substitutions. Defaults to the values in COSMIC
#'
#' @export
#' @import patchwork
#' @import ggplot2
plot_Signatures <- function(signatures,
                            lowCI = NULL,
                            highCI = NULL,
                            palette = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) {

  signatures <- as.matrix(signatures)
  names_sig <- rownames(signatures)
  df_plot <- data.frame(signatures) %>%
    dplyr::mutate(
      Channel = names_sig,
      Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                      function(x) paste0(x[c(1,3,7)], collapse = "")),
      Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                       function(x) paste0(x[c(3,4,5)], collapse = "")),
      Mutation = as.factor(Mutation)
    ) %>%
    tidyr::gather(key = "Sig", value = "Prob", -Channel, -Triplet, -Mutation)

  if(!is.null(lowCI) & !is.null(highCI)){
    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         tidyr::gather(key = "Sig", value = "lowCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         tidyr::gather(key = "Sig", value = "highCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig"))
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Triplet, y = Prob, fill = Mutation))+
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::facet_grid(Sig~Mutation, scales = "free", switch = "y") +
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_manual(values = palette)+
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, color = "gray35",
                                          vjust = .5, size = 6.5, margin = ggplot2::margin(t = -4)),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(0, "lines"),
      panel.spacing.y = ggplot2::unit(0, "lines"),
      # The next line makes y facet labels horizontal
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1)
    )

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey65")
  }

  return(p)

}

#' Plot the matrix of baseline values
#'
#' @param Theta Matrix of baseline activities
#' @param clust Cluster profiles of the observations.
#' @param nclust Number of clusters to split the data.
#'               Default is 10 for visualization. This is valid only if \code{clust = NULL}
#'
#' @export
#' @import patchwork
#' @import ggplot2
plot_ThetaBaseline <- function(Theta, clust = NULL, nclust = 10) {
  Theta <- t(Theta)
  if(is.null(clust)){
    clust <- stats::cutree(stats::hclust(dist(Theta)), k = nclust)
  } else {
    if(length(clust) != nrow(Theta)){
      stop("Length of clust must be equal to ncol(Theta)")
    }
  }
  Theta <- Theta[order(clust), ]
  pW <- Theta %>%
    as.data.frame() %>%
    dplyr::mutate(patient = rownames(Theta),
                  group =1,
                  clust = clust[order(clust)]) %>%
    tidyr::gather(key = "Signature", value = "weight", -patient, -group, -clust) %>%
    dplyr::mutate(Signature = as.factor(Signature),
                  patient = as.factor(patient))%>%
    dplyr::mutate(patient = ordered(patient, unique(patient))) %>%
    ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::geom_bar(ggplot2::aes(x=patient, y = weight, fill = Signature), stat = "identity") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(axis.text.x =  element_blank(),
                   axis.title.x = element_blank(),
                   legend.position = "right")+
    ggplot2::ylab("Baseline activities")+
    ggforce::facet_row(~clust, scales = "free", space = "free")

  return(pW)
}


