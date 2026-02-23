# Poisson process factorization

Code to reproduce figures, simulation results, and real data analysis results from the paper 

> Zito, A., Parmigiani, G. and Miller, J. W. (2025) - Poisson process factorization for mutational signature analysis with genomic covariates, arxiv: https://arxiv.org/abs/2510.26090 

## As a first step, install the SigPoisProcess package

```
install.packages("SigPoisProcess_0.1.0.tar.gz")
library(SigPoisProcess)
```

## Code to reproduce the results of the paper

### Figure 1 - Patterns of mutations along the genome

* Figure 1 panel a and - Aggregate mutations at the Mb scale, and relationship between mutations at the 2kb scale and signal from H3K9me3
  - `R/reproduce_Figure_1_panels_a_b.R`

### Simulation
* Figure 2, and Figure S6.1 to S6.6 in the Supplementary material - Main simulation
  - `R/main_Simulation_analysis.R`
* Figure S6.7 - Sparse indel simulation
  - `R/Simulation_sparsity_indels_suppl.R`
* Figure S7.2, S7.3 and S7.4 - Fixed strength vs Compressive hyperprior
  - `R/Simulation_fixed_vs_compressive.R`
* Figure S7.5 - Sensitivity to epsilon and K
  - `R/Simulation_sensitivity_epsilon_K.R`
* Figure S7.6 - Sensitivity to a and alpha
  - `R/Simulation_sensitivity_alpha_a_suppl.R`










