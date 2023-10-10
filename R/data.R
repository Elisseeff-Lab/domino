#' PBMC3K single cell RNA sequencing data
#'
#' A small example single cell data set from 10x genomcs
#'
#' @format ## `pbmc`
#' A Seurat object that has been through preliminary analysis, including normalization, scaling, PCA, UMAP, clustering, and cluster annotation
#' @keywords internal
#' 
#' @source <https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>
"pbmc"

#' Domino object created from PBMC3K single cell RNA sequencing data
#'
#' A small example domino object created with the `pbmc` data
#'
#' @format ## `pbmc_dom`
#' A domino object that has been built using SCENIC, CellphoneDB v4, and the domino2 package
#' @keywords internal
"pbmc_dom"