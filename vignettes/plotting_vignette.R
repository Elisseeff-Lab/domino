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

## ----setup, include = FALSE, eval = TRUE--------------------------------------
library(domino2)
library(patchwork)

## ----corheatmap---------------------------------------------------------------
cor_heatmap(dom, title = "PBMC R-TF Correlations")

## ----corheatmap-options, fig.show="hold", out.width = "50%"-------------------
cor_heatmap(dom, bool = TRUE, bool_thresh = 0.25)
cor_heatmap(dom, bool = FALSE, mark_connections = TRUE)

## ----corheatmap-subset--------------------------------------------------------
receptors <- c("CSF1R", "CSF3R", "CCR7", "FCER2")
tfs <- c("PAX5", "JUNB", "FOXJ3", "FOSB")
cor_heatmap(dom, feats = tfs, recs = receptors)

## ----corheatmap-nmf-args------------------------------------------------------
cor_heatmap(dom, title = FALSE, Rowv = NA, Colv = NA,
    main = "Heatmap Without Clustering")

## ----featheatmap--------------------------------------------------------------
feat_heatmap(dom)

## ----featheatmap-options, fig.show="hold", out.width = "50%"------------------
feat_heatmap(dom, min_thresh = 0.1, max_thresh = 0.6, 
    norm = TRUE, bool = FALSE)
feat_heatmap(dom, bool = TRUE)

