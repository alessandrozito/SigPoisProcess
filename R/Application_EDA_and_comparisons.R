# This file prepares the dataset for for the exploratory data analysis

# Load packages
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(foreach)
library(patchwork)
library(doParallel)
library(sigminer)

# Load SigPoisProcess
library(SigPoisProcess)

#
