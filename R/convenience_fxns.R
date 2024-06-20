#' @import plyr
#' @import methods
#'
NULL

#' Renames clusters in a domino object
#'
#' This function reads in a receptor ligand signaling database, cell level
#' features of some kind (ie. output from pySCENIC), z-scored single cell data,
#' and cluster id for single cell data, calculates a correlation matrix between
#' receptors and other features (this is transcription factor module scores if
#' using pySCENIC), and finds features enriched by cluster. It will return a
#' domino object prepared for [build_domino()], which will calculate a signaling
#' network.
#'
#' @param dom Domino object to rename clusters in
#' @param clust_conv Named vector of conversions from old to new clusters. Values are taken as new clusters IDs and names as old cluster IDs.
#' @param warning Logical. If TRUE, will warn if a cluster is not found in the conversion table. Default is FALSE.
#' @return A domino object with clusters renamed in all applicable slots.
#' @export
#' @examples 
#' data(pbmc_dom_built_tiny)
#' new_clust <- c("CD8_T_cell" = "CD8+ T Cells",
#'  "CD14_monocyte" = "CD14+ Monocytes", "B_cell" = "B Cells")
#' pbmc_dom_built_tiny <- rename_clusters(pbmc_dom_built_tiny, new_clust)
#'
rename_clusters <- function(dom, clust_conv, warning = FALSE) {
    if (is.null(dom@clusters)) {
        stop("There are no clusters in this domino object")
    }
    if (dom@misc$create) {
        dom@clusters <- plyr::revalue(dom@clusters, clust_conv, warn_missing = warning)
        colnames(dom@clust_de) <- plyr::revalue(colnames(dom@clust_de), clust_conv, warn_missing = warning)
        names(colnames(dom@clust_de)) <- c()
        colnames(dom@misc$cl_rec_percent) <- plyr::revalue(colnames(dom@misc$cl_rec_percent),
            clust_conv,
            warn_missing = warning
        )
    }
    if (dom@misc$build) {
        names(dom@linkages$clust_tf) <- plyr::revalue(names(dom@linkages$clust_tf), clust_conv, warn_missing = warning)
        names(dom@linkages$clust_rec) <- plyr::revalue(names(dom@linkages$clust_rec), clust_conv, warn_missing = warning)
        names(dom@linkages$clust_incoming_lig) <- plyr::revalue(names(dom@linkages$clust_incoming_lig),
            clust_conv,
            warn_missing = warning
        )
        names(dom@linkages$clust_tf_rec) <- plyr::revalue(names(dom@linkages$clust_tf_rec),
            clust_conv,
            warn_missing = warning
        )
        sig_ligands <- colnames(dom@signaling)
        sig_rec <- rownames(dom@signaling)
        # Remove L_ prefix
        sig_ligand_clust <- gsub("^L_", "", sig_ligands)
        # Remove R_ prefix
        sig_rec_clust <- gsub("^R_", "", sig_rec)
        new_lig_clust <- plyr::revalue(sig_ligand_clust, clust_conv, warn_missing = warning)
        new_rec_clust <- plyr::revalue(sig_rec_clust, clust_conv, warn_missing = warning)
        colnames(dom@signaling) <- paste0("L_", new_lig_clust)
        rownames(dom@signaling) <- paste0("R_", new_rec_clust)
        names(dom@cl_signaling_matrices) <- plyr::revalue(names(dom@cl_signaling_matrices),
            clust_conv,
            warn_missing = warning
        )
        for (cl in names(dom@cl_signaling_matrices)) {
            cl_sig_ligands <- colnames(dom@cl_signaling_matrices[[cl]])
            # Remove L_ prefix
            cl_sig_lig_clust <- gsub("^L_", "", cl_sig_ligands)
            cl_sig_lig_new <- plyr::revalue(cl_sig_lig_clust, clust_conv, warn_missing = warning)
            colnames(dom@cl_signaling_matrices[[cl]]) <- paste0("L_", cl_sig_lig_new)
        }
    }
    return(dom)
}


#' Convert Genes Using Table
#'
#' Takes a vector of gene inputs and a conversion table  and returns a 
#' converted gene table
#'
#' @param genes The genes to convert.
#' @param from  Gene symbol type of the input (ENSG, ENSMUSG, HGNC, MGI)
#' @param to    Desired gene symbol type for the output (HGNC, MGI)
#' @param conversion_table A data.frame with column names corresponding to gene symbol types (mm.ens, hs.ens, mgi, hgnc)
#' and rows corresponding to the gene symbols themselves
#' @return Data frame of genes with original and corresponding converted symbols
#' @keywords internal
#'
table_convert_genes <- function(genes, from, to, conversion_table) {
  # Check inputs:
  stopifnot(`Genes must be a vector of characters` = (is(genes, "character") & is(genes, "vector")))
  stopifnot(`From must be one of ENSMUSG, ENSG, MGI, or HGNC` = from %in% c(
    "ENSMUSG", "ENSG", "MGI",
    "HGNC"
  ))
  stopifnot(`To must be one of MGI or HGNC` = to %in% c("MGI", "HGNC"))
  stopifnot(`Conversion table must be provided with at least two of column names mm.ens, hs.ens, mgi and/or hgnc` = (is(
    conversion_table,
    "data.frame"
  ) & length(which(colnames(conversion_table) %in% c(
    "mm.ens", "hs.ens", "mgi",
    "hgnc"
  ))) > 1))
  if (from == "ENSMUSG") {
    col1 <- conversion_table$mm.ens
  }
  if (from == "ENSG") {
    col1 <- conversion_table$hs.ens
  }
  if (from == "MGI") {
    col1 <- conversion_table$mgi
  }
  if (from == "HGNC") {
    col1 <- conversion_table$hgnc
  }
  if (to == "MGI") {
    col2 <- conversion_table$mgi
  }
  if (to == "HGNC") {
    col2 <- conversion_table$hgnc
  }
  genesV2 <- cbind(col1[which(col1 %in% genes)], col2[which(col1 %in% genes)])
  return(genesV2)
}
