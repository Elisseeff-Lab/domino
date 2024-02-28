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

#' Resolve Names
#' @description
#' Return gene names of ligands with non-standard names
#' @param dom domino object
#' @param genes vector of gene names to resolve
#' @details Ligand names (which are stored in dom@linkages$rec_lig) do not always match the gene name
#' Search all provided names and output the gene name if the ligand name is non-standard
#' @return vector of length(genes) with applicable values replaced using dom@misc$rl_map information
#' @export
#' 
resolve_names <- function(dom, genes, rec_lig = "lig") {
  rl_map = dom@misc[["rl_map"]]
  if (rec_lig == "lig") {
    gene <- "L.gene"
    name <- "L.name"
  } else if (rec_lig == "rec") {
    gene <- "R.gene"
    name <- "R.name"
  } else {
    stop("rec_lig must be one of 'lig' or 'rec'.\n")
  }
  
  genes_resolved <- sapply(genes, function(l){
    int <- rl_map[rl_map[[name]] == l, ][1,] 
    if((int[[name]] != int[[gene]]) & !grepl("\\,", int[[gene]])){
      int[[gene]]
    } else { 
      int[[name]]
    }
  })
  return(genes_resolved)
} # resolve_names

#' Resolve Complexes
#' @description
#' Expand any complex names into their component gene names
#' @param dom domino object
#' @param genes vector of complex names to resolve
#' @details Ligand names (which are stored in dom@linkages$rec_lig) can refer to complexes that are
#' detailed in dom@linkages$complexes. Search all provided names and if the name is a complex, output component genes, otherwise output the name
#' @return list of length(genes), with names(list) == genes. List values are the same as the names if not in a complex and are the complex genes if they are in a complex.
#' @export
#' 
resolve_complexes <- function(dom, genes) {
  genes_list <- lapply(genes, function(l){
    if(l %in% names(dom@linkages$complexes)){
      return(dom@linkages$complexes[[l]])
    } else {
      return(l)
    }
  })
  names(genes_list) <- genes
  return(genes_list)
}

#' Get All Ligands
#' @description
#' Get all unique ligands present in dom@linkages$rec_lig
#' @param dom domino object
#' @param expressed_only logical indicating whether to subset ligands based on expression in dom@z_scores
#' @details Get all unique ligands in a domino object, expanding all ligand complexes as well. Optionally subset by expression.
#' @return vector of ligands if expressed_only = T, list if F
get_all_reclig <- function(dom, expressed_only = T, rec_lig = "rec") {
  
  ### Vector of all ligands expressed
  if (rec_lig == "rec") {
    all <- unlist(dom@linkages$rec_lig)
    resolve_rec_lig <- "lig"
  } else if (rec_lig == "lig") {
    if (!"lig_rec" %in% names(dom@linkages)) stop("Must run invert_rec_lig_linkages if rec_lig is set to 'lig'.\n")
    all <- unlist(dom@linkages$lig_rec)
    resolve_rec_lig <- "rec"
  } else {
    stop("lig_rec must be one of 'lig' or 'rec")
  }
  
  all <- unique(all)
  all <- all[!all == ""]
  
  ### Resolve non-standard ligand names
  all_names_resolved <- resolve_names(dom, all, rec_lig = resolve_rec_lig)
  all_names_resolved <- unique(all_names_resolved)
  
  ### Resolve complexes
  if(length(dom@linkages$complexes) > 0){
    all_complexes_resolved_list <- resolve_complexes(dom, all_names_resolved)
    all_names_resolved <- unlist(all_complexes_resolved_list)
  }
  
  if (expressed_only) {
    ### Subset for ligands expressed in the data
    genes <- intersect(all_names_resolved, rownames(dom@z_scores))
  } else {
    genes <- all_names_resolved
    #genes <- all_complexes_resolved_list
  } 
  
  out_ls <- list("genes" = genes, "complex" = all_complexes_resolved_list)
  return(out_ls)
  
} # get_all_reclig
