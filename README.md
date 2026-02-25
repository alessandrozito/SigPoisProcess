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

### Simulation study in the paper

* Figure 2, and Figure S1, Table S1 and S2 in the Supplementary material - Main simulation
  - `R/Simulaton.R` 
  - `R/Simulation_functions.R`

### Application 1 - de novo signature extraction

* Files to assemble the data and run the model
  - `R/Load_ICGC_BreastAdenoCA_data.R`
  - `R/Application_denovo.R`

* Figure 3 and Figure 4 - prediction of mutations at Mb rate and posterior estimates
  - `R/Reproduce_figures_Application_denovo.R`

### Application 2 - Fixed signatures analysiss

* Files to assemble the data and run the model
  - `R/Application_refit.R`

* Figure 5 - output of the posterior quantities
  - `R/Reproduce_figures_Application_refit.R`


