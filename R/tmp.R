


cosmic <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37
SigPoisProcess::plot_96SBSsig(cosmic[, c("SBS1", "SBS2", "SBS3")])

library(tidyverse)
library(rtracklayer)

bed <- as.data.frame(read.table("data/Breast/ENCODE_HistoneModification/ENCFF908LTV.bed.gz",
                                header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                  "signalValue", "pValue", "qValue", "peak")

big <- import("data/ENCODE_pvalues/ENCODE_rawdata/Otlu_bigWig/ENCFF007QRU.bigWig")

big

bed <- bed %>%
  arrange(chrom, chromStart)

head(bed)

file <- import("data/ENCODE_pvalues/Cancer_covariates/Breast-Cancer_male_POLR2A.bigWig")

sc <- file$score
signal <- (sc - quantile(sc, 0.001))/(quantile(sc, 0.99) - quantile(sc, 0.001)) * 10
table(signal > 10)
signal[signal > 10] <- 10
signal[signal < 0] <- 0
hist(signal, breaks = 50)
area = sum(width(file) *signal)


sc <- signal[1e6:2e6]
plot(sc, type = "l")

plot(bed$score, bed$signalValue)
plot(bed$score, bed$qValue)
plot(bed$score, bed$qValue)
fit <- lm(score ~ qValue, data = bed)

coef(fit)

- coef(fit)[1]/coef(fit)[2]
(1 - coef(fit)[1])/coef(fit)[2]

min(bed$qValue)
max(bed$qValue)

hist(bed$qValue)
quantile(bed$qValue)
quantile(bed$score)


bed[10^(-bed$qValue) > 0.05, ]
coef()
hist(bed$score)
