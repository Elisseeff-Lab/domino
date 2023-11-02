## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.cap="",
  tidy = "styler"
)

load("../data/pbmc_dom.rda")
dom <- pbmc_dom

## ----setup, include = FALSE---------------------------------------------------
library(domino2)

## ----include = FALSE----------------------------------------------------------
# More of Jacob's stuff that I'm cutting:
# Active transcription factors in each cluster are determined by conducting wilcoxon rank sum tests for each transcription factor where the trascription factor activity scores amongst all cells in the cluster are tested against the activity scores of all cells outside of the cluster. The p-value for the one-sided test for greater activity within the cluster compared to outside is stored in pbmc_dom\@clust_de. Linkage between receptors and transcription factors is assessed by Spearman correlation between transcription factor activity scores and scaled expression of receptor-encoding genes across all cells in the data set. Spearman coefficients are stored in pbmc_dom\@cor.

## -----------------------------------------------------------------------------
Sys.Date()
sessionInfo()

