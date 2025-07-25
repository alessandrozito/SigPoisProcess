

reshape2::melt(Loads)

MaxProbIds <- reshape2::melt(apply(Probs, c(2,3), mean))

Xbar <- object$Xbar

res <- reshape2::melt(Loads * Xbar) %>%
  dplyr::mutate(in_present = value > 0.002) %>%
  dplyr::group_by(Var1, Var3) %>%
  dplyr::summarise(n = mean(in_present))

library(tidyverse)
plot(object)
