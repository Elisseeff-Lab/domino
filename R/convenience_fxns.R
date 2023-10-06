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
#' @export 
#' 
rename_clusters <- function(dom, clust_conv) {
    if (is.null(dom@clusters)) {
        stop("There are no clusters in this domino object")
    }
    if (dom@misc$create) {
        dom@clusters <- plyr::revalue(dom@clusters, clust_conv)
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
#' Extracts all features, receptors, or ligands present in a signaling network.
#' 
#' This function collates all of the features, receptors, or ligands found in a
#' signaling network anywhere in a list of clusters. This can be useful for
#' comparing signaling networks across two separate conditions. In order to run
#' this [build_domino()] must be run on the object previously.
#' 
#' @param dom Domino object containing a signaling network (i.e. [build_domino()] run)
#' @param return String indicating where to collate 'features', 'receptors', or 'ligands'. If 'all' then a list of all three will be returned.
#' @param clusters Vector indicating clusters to collate network items from. If left as NULL then all clusters will be included.
#' @return A vector containing all features, receptors, or ligands in the data set or a list containing all three.
#' @export 
#' 
collate_network_items <- function(dom, clusters = NULL, return = NULL) {
    if (!dom@misc[["build"]]) {
        stop("Please run domino_build prior to generate signaling network.")
    }
    if (is.null(clusters) & is.null(dom@clusters)) {
        stop("There are no clusters in this domino object. Please provide clusters.")
    }
    if (is.null(clusters)) {
        clusters <- levels(dom@clusters)
    }
    # Get all enriched TFs and correlated + expressed receptors for specified clusters
    all_recs <- c()
    all_tfs <- c()
    all_ligs <- c()
    for (cl in clusters) {
        all_recs <- c(all_recs, unlist(dom@linkages$clust_tf_rec[[cl]]))
        tfs <- names(dom@linkages$clust_tf_rec[[cl]])
        tf_wo_rec <- which(sapply(dom@linkages$clust_tf_rec[[cl]], length) == 0)
        if (length(tf_wo_rec > 0)) {
            tfs <- tfs[-tf_wo_rec]
        }
        all_tfs <- c(all_tfs, tfs)
        all_ligs <- c(all_ligs, rownames(dom@cl_signaling_matrices[[cl]]))
    }
    all_recs <- unique(all_recs)
    all_tfs <- unique(all_tfs)
    all_ligs <- unique(all_ligs)
    # Make list and return whats asked for
    list_out <- list(features = all_tfs, receptors = all_recs, ligands = all_ligs)
    if (is.null(return)) {
        return(list_out)
    } else {
        return(list_out[[return]])
    }
}
#' Convert Genes Using Table
#' 
#' Takes a vector of gene inputs and returns converted gene table
#' 
#' @param genes The genes to convert.
#' @param from  Gene symbol type of the input (ENSG, ENSMUSG, HGNC, MGI)
#' @param to    Desired gene symbol type for the output (HGNC, MGI)
#' @param conversion_table A data.frame with column names corresponding to gene symbol types (mm.ens, hs.ens, mgi, hgnc)
#' and rows corresponding to the gene symbols themselves
#' @return Data frame of genes with original and corresponding converted symbols
#' @export
table_convert_genes <- function(genes, from, to, conversion_table) {
    # Check inputs:
    stopifnot(`Genes must be a vector of characters` = (is(test, "character") & is(test, "vector")))
    stopifnot(`From must be one of ENSMUSG, ENSG, MGI, or HGNC` = from %in% c("ENSMUSG", "ENSG", "MGI",
        "HGNC"))
    stopifnot(`To must be one of MGI or HGNC` = to %in% c("MGI", "HGNC"))
    stopifnot(`Conversion table must be provided with at least two of column names mm.ens, hs.ens, mgi and/or hgnc` = (is(conversion_table,
        "data.frame") & length(which(colnames(conversion_table) %in% c("mm.ens", "hs.ens", "mgi",
        "hgnc"))) > 1))
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
## Back up gene table for back up function: mouse_human_genes =
## read.csv('http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt',sep='\t')
