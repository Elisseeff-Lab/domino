#' SCENIC AUC subset
#'
#' A subset of SCENIC AUCs as applied to PBMC data.
#'
#' @format A list of:
#' \describe{
#'  \item{auc_tiny}{A subset of SCENIC AUCs}
#'  \item{regulons_tiny}{A subset of SCENIC regulons}
#' }
#'
#' @source <https://zenodo.org/records/10951634/files>
#' @usage data("SCENIC")
"SCENIC"


#' PBMC RNAseq data subset
#'
#' A subset of the results of PBMC RNA-seq data.
#'
#' @format A list of::
#' \describe{
#'  \item{RNA_count_tiny}{A subset of PBMC RNA-seq data: counts assay}
#'  \item{RNA_zscore_tiny}{A subset of PBMC RNA-seq data: zscore assay}
#'  \item{clusters_tiny}{A subset of PBMC RNA-seq data: clusters as defined by cell_type}
#' }
#'
#' @source <https://zenodo.org/records/10951634/files/pbmc3k_sce.rds>
#' @usage data("PBMC")
"PBMC"


#' CellPhoneDB subset
#'
#' A list of four subsets of CellPhoneDB data.
#'
#'
#' @format A list of:
#' \describe{
#'  \item{genes_tiny}{A subet of CellPhoneDB gene_input.csv}
#'  \item{proteins_tiny}{A subset of CellPhoneDB protein_input.csv}
#'  \item{complexes_tiny}{A subset of CellPhoneDB complex_input.csv}
#'  \item{interactions_tiny}{A subset of CellPhoneDB interaction_input.csv}
#' }
#' 
#' @source <https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz>
#' @usage data("CellPhoneDB")
"CellPhoneDB"
