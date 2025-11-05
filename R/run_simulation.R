# This file runs the first batch of simulations
library(tidyverse)
library(SigPoisProcess)
library(corrplot)

# Source relevant files
source("R/simulate_data.R")

# Generate covariates with no correlation
set.seed(10)
data_nocorr <- generate_MutationData()
saveRDS(data_nocorr, file = "data/temp_simulation.rds.gzip", compress = "gzip")

# Generate covariates with some correlation
set.seed(10)
data_corr <- generate_MutationData(corr = "onion")
saveRDS(data_corr, file = "data/temp_simulation_corr.rds.gzip", compress = "gzip")


#----------------------------------------------
# Part 1 - Uncorrelated covariates
#----------------------------------------------
data_nocorr <- readRDS("~/SigPoisProcess/data/temp_simulation.rds.gzip")
corrplot(cor(as.matrix(mcols(data_nocorr$gr_SignalTrack))))

#--- Model is correctly specified with covariates
set.seed(10)
res_correct <- SigPoisProcess(gr_Mutations = data_nocorr$gr_Mutations,
                      df_areas = data_nocorr$df_areas_norm,
                      K = 10,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 10000))
plot(res_correct)

#--- Exclude the first two important covariates, keep the others
gr_Mutations <- data_nocorr$gr_Mutations
mcols(gr_Mutations) <- mcols(gr_Mutations)[, -(c(1:2) + 2)]
df_areas <- data_nocorr$df_areas_norm[, -(c(1:2) + 1)]
set.seed(10)
res_misp1 <- SigPoisProcess(gr_Mutations = gr_Mutations,
                      df_areas = df_areas,
                      baseline_total = 2e6,
                      K = 10,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 10000))
plot(res_misp1)

#--- Exclude all covariates
gr_Mutations <- data_nocorr$gr_Mutations
mcols(gr_Mutations) <- mcols(gr_Mutations)[, -(c(1:5) + 2)]
df_areas <- data_nocorr$df_areas_norm[, -(c(1:5) + 1)]
set.seed(10)
res_misp_all <- SigPoisProcess(gr_Mutations = gr_Mutations,
                            df_areas = df_areas,
                            baseline_total = 2e6,
                            K = 10,
                            prior_params = list(a = 1.1, alpha = 1.1),
                            controls = SigPoisProcess.controls(maxiter = 10000))
plot(res_misp_all)
match_to_RefSigs(res_misp_all$Signatures, ref = data_nocorr$R)
#----------------------------------------------
# Part 2 - Correlated covariates
#----------------------------------------------
data_corr <- readRDS("~/SigPoisProcess/data/temp_simulation_corr.rds.gzip")
corrplot(cor(as.matrix(mcols(data_corr$gr_SignalTrack))))

#--- Model is correctly specified with covariates
set.seed(10)
res_correct_corr <- SigPoisProcess(gr_Mutations = data_corr$gr_Mutations,
                              df_areas = data_corr$df_areas_norm,
                              K = 10,
                              prior_params = list(a = 1.1, alpha = 1.1),
                              controls = SigPoisProcess.controls(maxiter = 1000))
plot(res_correct_corr)

#--- Exclude the first two important covariates, keep the others
gr_Mutations <- data_corr$gr_Mutations
mcols(gr_Mutations) <- mcols(gr_Mutations)[, -(c(1:2) + 2)]
df_areas <- data_corr$df_areas_norm[, -(c(1:2) + 1)]
set.seed(10)
res_misp1_corr <- SigPoisProcess(gr_Mutations = gr_Mutations,
                            df_areas = df_areas,
                            baseline_total = 2e6,
                            K = 10,
                            prior_params = list(a = 1.1, alpha = 1.1),
                            controls = SigPoisProcess.controls(maxiter = 10000))
plot(res_misp1_corr)

#--- Exclude all covariates
gr_Mutations <- data_corr$gr_Mutations
mcols(gr_Mutations) <- mcols(gr_Mutations)[, -(c(1:5) + 2)]
df_areas <- data_corr$df_areas_norm[, -(c(1:5) + 1)]
set.seed(10)
res_misp_all_corr <- SigPoisProcess(gr_Mutations = gr_Mutations,
                               df_areas = df_areas,
                               baseline_total = 2e6,
                               K = 10,
                               prior_params = list(a = 1.1, alpha = 1.1),
                               controls = SigPoisProcess.controls(maxiter = 10000))
plot(res_misp_all_corr)

data <- readRDS("../data/temp_simulation.rds.gzip")
res <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
                      df_areas = data$df_areas_norm,
                      K = 10,
                      prior_params = list(a = 1.1, alpha = 1.1),
                      controls = SigPoisProcess.controls(maxiter = 1000))


res2 <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
                      df_areas = data$df_areas_norm,
                      K = 12, compressType = "sig",
                      prior_params = list(a = 1.1, alpha = 1.1, a0 = 1.1 * 50 + 1, b0 = 1.1 * 0.001 * 50),
                      controls = SigPoisProcess.controls(maxiter = 5000))
plot(res2)
round(res2$Mu, 3)

res3 <- SigPoisProcess(gr_Mutations = data$gr_Mutations,
                       df_areas = data$df_areas_norm,
                       K = 12, compressType = "sig_cov",
                       prior_params = list(a = 1.1, alpha = 1.1, a0 = 1.1 * 50 + 1, b0 = 1.1 * 0.001 * 50),
                       controls = SigPoisProcess.controls(maxiter = 5000))

plot(res3)

round(res2$Betas, 3)
round((res2$Betas * res2$Xbar)[, , "Encode06"], 5)
dim(res2$Betas)

apply(res2$Betas * res2$Xbar, 1, mean)

tmp <- res2$Betas * res2$Xbar
j <- 1
round(colSums(tmp[, j, ])/sum(tmp[,j,]), 5)

Xnorm <- tmp
Attr <- data.frame()
for(j in 1:50){
  Attr <- rbind(Attr,
                data.frame("patient" = colnames(tmp)[j],
                           t(colSums(tmp[, j, ])/sum(tmp[,j,]))))
  Xnorm[, j, ] <- Xnorm[, j, ]/sum(Xnorm[, j, ])
}

apply(Xnorm > 0.005, c(1,3), sum)/50

Attr %>%
  tidyr::gather(key = "Covariate", value = "value", -patient)


sum(apply(tmp, 2, function(x) x/sum(x))[,1,])


round(apply(res2$Betas * res2$Xbar, 1, mean), 3)
round(res2$Mu[,1], 3)

