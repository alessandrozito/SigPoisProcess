#' @export
plot_avgLoads <- function(Betas, Xbar, AttrSummary, cutoff_excluded = 0.002){
  # Average values for Beta
  AvgBeta <- apply(Betas * Xbar, c(1,3), mean)
  rownames(AvgBeta) <- rownames(Betas)
  df <- reshape2::melt(as.matrix(AvgBeta))
  colnames(df) <- c("Signature", "Covariate", "Value")

  df <- df %>%
    mutate(
      Covariate = factor(Covariate, levels=unique(Covariate)),
      AltCol = ifelse((as.numeric(Covariate) %% 2) == 0, "gray93", "gray97")) %>%
    left_join(AttrSummary, by = c("Signature", "Covariate")) %>%
    mutate(Value = dplyr::case_when(Value < cutoff_excluded ~ NA_real_,
                                    TRUE ~ Value))
  df$Signature <- as.factor(df$Signature)
  df$Covariate <- as.factor(df$Covariate)

  p <- ggplot2::ggplot(df, aes(x = Covariate, y = Signature)) +
    geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
    scale_fill_manual(
      values = c("gray93" = "gray93", "gray97" = "gray97"),
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(fill = Value, size = Fraction),
               shape = 21, color = "grey30", stroke = 0.2,
               na.rm = TRUE) +
    scale_fill_gradientn(name = "Average\nloading",
                         colours =c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                                    "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
    scale_y_discrete(limits = rev(levels(df$Signature))) +
    scale_size_continuous(name = "Fraction \nof total\nmutations")+
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
      axis.title = element_blank(),
      #axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.ticks.x = element_blank())

  return(p)
}


#' @export
plot.SigPoisProcess <- function(object){
  p_sig <- plot_96SBSsig(object$Signatures, remove_strip_text_y = TRUE)
  p_Loads <- plot_avgLoads(object$Betas, object$Xbar, object$AttrSummary)
  p_sig + p_Loads
}

#' @export
plot_96SBSsig <- function(signatures,
                          lowCI = NULL,
                          highCI = NULL,
                          returnList = FALSE,
                          remove_strip_text_y = FALSE,
                          palette = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) {
  signatures <- as.matrix(signatures)
  names_sig <- rownames(signatures)
  df_plot <- data.frame(signatures) %>%
    dplyr::mutate(Channel = names_sig,
                  Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                  function(x) paste0(x[c(1,3,7)], collapse = "")),
                  Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                   function(x) paste0(x[c(3,4,5)], collapse = "")),
                  Mutation = as.factor(Mutation)) %>%
    gather(key = "Sig", value = "Prob", -Channel, -Triplet, -Mutation)

  if(!is.null(lowCI) & !is.null(highCI)){
    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         gather(key = "Sig", value = "lowCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         gather(key = "Sig", value = "highCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig"))
  }

  p_list <- vector(mode = "list", length = ncol(signatures))
  for(j in 1:length(p_list)){
    p <-   ggplot(df_plot %>% filter(Sig == colnames(signatures)[j]),
                  aes(x = Triplet, y = Prob, fill = Mutation))+
      geom_bar(stat = "identity", width = 0.5) +
      #geom_hline(yintercept = -0.001)+
      #geom_segment(x = 1, xend = 1, y = -0.01, yend = Inf, linewidth = 0.3)+
      facet_grid(Sig~Mutation, scales = "free")+
      theme_minimal()+
      scale_fill_manual(values = palette)+
      theme(legend.position = "none",aspect.ratio = 1.2,
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            plot.margin = margin(0, 0, 0, 0),
            #axis.text.x = element_text(angle = 90, color = "gray35",
            #                          vjust = .5, size = 6.5, margin = margin(t = -4)),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            panel.spacing.x=unit(0.1, "lines"),
            panel.spacing.y=unit(0,"lines"),
            strip.background=element_blank(),
            axis.line.x.bottom = element_line(linewidth = 0.2),
            axis.line.y.left = element_line(linewidth = 0.2),
            strip.text.x = element_blank())
    if(!is.null(lowCI) & !is.null(highCI)){
      p <- p +
        geom_linerange(aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey65")
    }
    if(remove_strip_text_y){
      p <- p + theme(strip.text.y = element_blank())
    }
    p_list[[j]]  <- p
  }
  if(returnList){
    return(p_list)
  }
  p_all <- ggpubr::ggarrange(plotlist = p_list, ncol = 1, common.legend = TRUE, legend = "none")
  return(p_all)
}






