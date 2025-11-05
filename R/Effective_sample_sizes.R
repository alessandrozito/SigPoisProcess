library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(foreach)
library(sigminer)
library(patchwork)
library(doParallel)
library(GenomicRanges)
library(AnnotationDbi)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(org.Hs.eg.db); #BiocManager::install("org.Hs.eg.db")
devtools::document()

get_PosteriorEffectiveSize <- function(chain, burnin = 1500){
  if(is.null(dim(chain))){
    coda::effectiveSize(chain[-c(1:burnin)])
  } else if(length(dim(chain)) == 2) {
    apply(chain[-c(1:burnin), ], 2, function(x) coda::effectiveSize(x))
  } else if (length(dim(chain)) == 3){
    apply(chain[-c(1:burnin), ,], c(2,3), function(x) coda::effectiveSize(x))
  }
}



#----- The tumor data
file_data <- "~/SigPoisProcess/data/ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip"
data <- readRDS(file_data)


# ------ Part 1: de novo analysis
outMCMC <- readRDS("~/SigPoisProcess/output/Application_Breast/MapSolutions/Replicate03/MCMCSolution.rds.gzip")
#### Calculate the effective sample sizes
mu <- colMeans(outMCMC$chain$MUchain[-c(1:7600), ])

mat_allCopy <- t(colSums(data$CopyTrack))[rep(1, sum(mu > 0.01)), ]
ThetaChain_adj <- outMCMC$chain$THETAchain[-c(1:7600), mu > 0.01, ]
array_AllCopy <- array(NA, dim = c(3000, 9, 113))
for(s in 1:3000) ThetaChain_adj[s,,] <- mat_allCopy * ThetaChain_adj[s, , ]

REffects <- get_PosteriorEffectiveSize(outMCMC$chain$SIGSchain[, , mu > 0.01], burnin = 7600)
ThetaEffects <- get_PosteriorEffectiveSize(ThetaChain_adj, burnin = 0)
BetasEffects <- get_PosteriorEffectiveSize(outMCMC$chain$BETASchain[, , mu > 0.01], burnin = 7600)
Sigma2Effects <- get_PosteriorEffectiveSize(outMCMC$chain$SIGMA2chain[, mu > 0.01], burnin = 7600)
MuEffects <- get_PosteriorEffectiveSize(outMCMC$chain$MUchain[, mu > 0.01], burnin = 7600)
lpEffects <- get_PosteriorEffectiveSize(outMCMC$chain$logPostchain, burnin = 7600)


effects_df <- tibble(
  value = c(REffects, ThetaEffects, BetasEffects, MuEffects),
  quantity = factor(rep(
    c("R", "Theta([0, T))", "B", "mu"),
    times = c(length(REffects), length(ThetaEffects), length(BetasEffects),
              length(MuEffects))
  ),
  levels = c("R", "Theta([0, T))", "B", "mu"))
)

p_effective <- ggplot(effects_df, aes(x = "", y = value)) +
  geom_boxplot() +
  facet_wrap(~ quantity, nrow = 1, scales = "free") +
  labs(x = "", y = "Effective Sample Size") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_logPost <- ggplot(data = data.frame(iter = 1:3000,
                                      logPosterior = outMCMC$chain$logPostchain[-c(1:7600)]))+
  geom_line(aes(x = iter, y = logPosterior)) +
  theme_bw()+
  facet_wrap(~ "logposterior trace", nrow = 1) +
  xlab("post-burnin iteration")

p_ess <- p_effective + p_logPost + plot_layout(nrow = 1, widths =  c(2, 1))

# ------ Part 2: fixed signatures
outMCMCFixed <- readRDS("~/SigPoisProcess/output/Application_Breast/SemiSupervised/MCMCSolution_SemiSupervised_noMap.rds.gzip")
burnin <- 2000

# ---- recalculate mu for this solution (post-burnin)
mu <- colMeans(outMCMCFixed$chain$MUchain[-c(1:burnin), ])

# ---- Compute mat_allCopy and adjusted ThetaChain as before
mat_allCopy <- t(colSums(data$CopyTrack))[rep(1, sum(mu > 0.01)), ]
ThetaChain_adj <- outMCMCFixed$chain$THETAchain[-c(1:burnin), mu > 0.01, ]
for (s in 1:dim(ThetaChain_adj)[1]) {
  ThetaChain_adj[s,,] <- mat_allCopy * ThetaChain_adj[s,,]
}

# ---- Calculate effective sample sizes
ThetaEffects <- get_PosteriorEffectiveSize(ThetaChain_adj, burnin = 0)
BetasEffects <- get_PosteriorEffectiveSize(outMCMCFixed$chain$BETASchain[,, mu > 0.01], burnin = burnin)
Sigma2Effects <- get_PosteriorEffectiveSize(outMCMCFixed$chain$SIGMA2chain[, mu > 0.01], burnin = burnin)
MuEffects <- get_PosteriorEffectiveSize(outMCMCFixed$chain$MUchain[, mu > 0.01], burnin = burnin)
lpEffects <- get_PosteriorEffectiveSize(outMCMCFixed$chain$logPostchain, burnin = burnin)

# ---- Make tidy data frame
effects_df <- tibble(
  value = c(ThetaEffects, BetasEffects, MuEffects),
  quantity = factor(rep(
    c("Theta([0, T))", "B", "mu"),
    times = c(length(ThetaEffects), length(BetasEffects), length(MuEffects))
  ),
  levels = c("Theta([0, T))", "B", "mu"))
)

# ---- Effective sample size boxplot
p_effective <- ggplot(effects_df, aes(x = "", y = value)) +
  geom_boxplot() +
  facet_wrap(~ quantity, nrow = 1, scales = "free") +
  labs(x = "", y = "Effective Sample Size") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# ---- Log posterior trace plot
p_logPost <- ggplot(
  data = data.frame(
    iter = 1:(length(outMCMCFixed$chain$logPostchain) - burnin),
    logPosterior = outMCMCFixed$chain$logPostchain[-c(1:burnin)]
  )
) +
  geom_line(aes(x = iter, y = logPosterior)) +
  theme_bw() +
  facet_wrap(~ "logposterior trace", nrow = 1) +
  xlab("post-burnin iteration")


p_ess_fixed <- p_effective + p_logPost + plot_layout(nrow = 1, widths = c(2, 1))



p_ess

p_ess_fixed

