
library(tidyverse)
library(SigPoisProcess)
# Make the following plots.

# 0 - check that loadings and signatures are matched, and calculate
# the RMSE
data <- readRDS("data/temp_simulation_corr.rds.gzip")
set.seed(42)
object <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
                         df_areas = data$df_areas_norm,
                         baseline_total = 2e6,
                         K = 10,
                         prior_params = list(a = 1.1, alpha = 1.1),
                         controls = SigPoisProcess.controls(maxiter = 5000))
plot(object)

TrueSigs <- data$R


match_MutSign <- function(R_true, R_hat) {
  # Need to make sure that R_hat and R_true are the same
  k_true <- ncol(R_true)
  k_hat <- ncol(R_hat)
  k_tot <- max(c(k_hat, k_true))
  I <- nrow(R_true)
  mat0 <- matrix(100, nrow = I, ncol = abs(k_hat - k_true)) #
  if (k_hat > k_true) {
    colnames(mat0) <- paste0("new_extra", 1:ncol(mat0))
    R_true <- cbind(R_true, mat0)
  } else if (k_hat < k_true) {
    R_hat <- cbind(R_hat, mat0)
  }

  # Match mutational signatures using the Hungarian algorithm
  CosMat <- matrix(1, k_tot, k_tot)
  for (i in 1:k_tot) {
    for (j in 1:k_tot) {
      CosMat[i, j] <- 1 - lsa::cosine(R_true[, i], R_hat[, j])
    }
  }
  match <- RcppHungarian::HungarianSolver(CosMat)$pairs[, 2]
  R_hat_matched <- R_hat[, match]

  # Change 100 with 0s
  R_hat_matched[R_hat_matched == 100] <- 0
  R_true[R_true == 100] <- 0

  return(list("R_hat" = R_hat_matched, "R_true" = R_true, "match" = match))
}

X_bar <- array(NA, dim = dim(data$Betas))
for (k in 1:nrow(data$Betas)) X_bar[k, , ] <- as.matrix(data$df_areas[, 2:6])
Betas_true <- data$Betas * X_bar
Sigs_true <- data$R

postProcess_model <- function(object, Sigs_true, Betas_true){
  Mu_hat <- object$Mu
  ids <- rowSums(Mu_hat > 1.5 * 0.001) > 0
  Mu_hat <- Mu_hat[ids, ]
  R_hat <- object$Signatures[, ids]
  Betas_hat <- object$Betas * object$Xbar
  Betas_hat <- Betas_hat[ids, ,]
  match <- match_MutSign(R_hat = R_hat, R_true = Sigs_true)$match
  R_hat <- R_hat[, match]
  Betas_hat <- Betas_hat[match, , -1]
  Mu_hat <- Mu_hat[match, ]

  i <- 7
  plot(R_hat[, i], Sigs_true[, i] )

  i <- 2
  plot(Betas_hat[, , i], Betas_true[, , i])
  abline(a = 0, b = 1)
  points(Betas_hatMAP[, , 1], Betas_true[, , 1], col = "red")
  abline(a = 0, b = 1)
}


grX <- data$gr_SignalTrack
grX_norm <- grX
mcols(grX_norm) <- apply(mcols(grX_norm), 2, function(x) x/mean(x))

J <- 50
CumSums <- apply(cbind("baseline" = 1, as.matrix(mcols(grX_norm))), 2, cumsum) * 100
CumSums_all <- array(NA, dim = c(nrow(CumSums), J, ncol(CumSums)))
for(j in 1:J) CumSums_all[, j, ] <- CumSums


estimateVaryingMutations <- function(object, CumSums_all){
  R <- object$Signatures
  Betas <- object$Betas
  EstTotals <- array(NA, dim = c(nrow(CumSums_all), ncol(CumSums_all), nrow(R)))
  for(j in 1:J) EstTotals[, j, ] <- t(R %*% Betas[, j, ] %*% t(CumSums_all[, j, ]))

  plot(EstTotals[, 4, 5], type = "l")

  tmp <- apply(MutTotals, c(1,2), sum)

  plot(tmp[, 3], type = "l")
}

# The real one
df_muts <- data.frame(data$gr_Mutations)
df_muts$region <- subjectHits(findOverlaps(data$gr_Mutations, data$gr_SignalTrack))
df_muts$region <- factor(df_muts$region, levels = 1:20000)
df_muts <- df_muts %>%
  group_by(region, sample, channel) %>%
  summarise(n = n())
df_muts$region <- as.numeric(df_muts$region)

MutTotals <- array(0, dim = c(nrow(CumSums_all), ncol(CumSums_all), nrow(R)),
                   dimnames = list(rownames(CumSums_all), colnames(Betas), rownames(R)))
MutTotals[df_muts$region, df_muts$sample, df_muts$channel] <- df_muts$n

sp <- df_muts$sample[i]
for(i in 1:nrow(df_muts)) {
  MutTotals[df_muts$region[i], df_muts$sample[i], df_muts$channel[i]] <- df_muts$n[i]
}




df_muts$sample
data$gr_SignalTrack
data$gr_Mutations







#' @export
estimateTotalCounts <- function(object, ...){
  R <- object$Signatures
  Betas <- object$Betas
  X_total <- object$Xbar[1, , ]
  MatOut <-sapply(1:ncol(Betas), function(j) R %*% Betas[, j, ] %*% X_total[j, ], simplify = TRUE)
  colnames(MatOut) <- colnames(Betas)
  rownames(MatOut) <- rownames(R)
  return(MatOut)
}





# 1 - variation along the genome




# 2 - Point plots with signatures

# 3 - mutation-specific plots



