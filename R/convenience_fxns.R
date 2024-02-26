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
#' @return A domino object with clusters renamed in all applicable slots.
#' @keywords internal
#' @export
#' @examples 
#' new_clust <- c("CD8_T_cell" = "CD8+ T Cells",
#'  "CD14_monocyte" = "CD14+ Monocytes", "B_cell" = "B Cells")
#' pbmc_dom_built_tiny <- rename_clusters(domino2:::pbmc_dom_built_tiny, new_clust)
#'
rename_clusters <- function(dom, clust_conv) {
  if (is.null(dom@clusters)) {
    stop("There are no clusters in this domino object")
  }
  if (dom@misc$create) {
    dom@clusters <- revalue(dom@clusters, clust_conv)
    colnames(dom@clust_de) <- clust_conv
    names(colnames(dom@clust_de)) <- c()
    colnames(dom@misc$cl_rec_percent) <- clust_conv
  }
  if (dom@misc$build) {
    names(dom@linkages$clust_tf) <- clust_conv
    names(dom@linkages$clust_rec) <- clust_conv
    names(dom@linkages$clust_incoming_lig) <- clust_conv
    names(dom@linkages$clust_tf_rec) <- clust_conv
    colnames(dom@signaling) <- paste0("L_", clust_conv)
    rownames(dom@signaling) <- paste0("R_", clust_conv)
    names(dom@cl_signaling_matrices) <- clust_conv
    for (cl in clust_conv) {
      colnames(dom@cl_signaling_matrices[[cl]]) <- paste0("L_", clust_conv)
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

make_rl_reading <- function(rl_map) {
  
  rl_reading <- NULL
  for (i in 1:nrow(rl_map)) {
    rl <- list()
    inter <- rl_map[i, ]
    p <- ifelse(inter[["type_A"]] == "R", "A", "B")
    q <- ifelse(p == "A", "B", "A")
    R.gene <- inter[[paste0("gene_", p)]]
    L.gene <- inter[[paste0("gene_", q)]]
    rl[["R.gene"]] <- R.gene
    rl[["L.gene"]] <- L.gene
    if (paste0("uniprot_", p) %in% names(inter)) {
      rl[["R.uniprot"]] <- inter[[paste0("uniprot_", p)]]
    }
    if (paste0("uniprot_", q) %in% names(inter)) {
      rl[["L.uniprot"]] <- inter[[paste0("uniprot_", q)]]
    }
    if (paste0("name_", p) %in% names(inter)) {
      rl[["R.name"]] <- inter[[paste0("name_", p)]]
    }
    if (paste0("name_", q) %in% names(inter)) {
      rl[["L.name"]] <- inter[[paste0("name_", q)]]
    }
    rl <- as.data.frame(rl)
    rl_reading <- rbind(rl_reading, rl)
  }
  return(rl_reading)
} 

