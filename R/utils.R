#' Access database
#'
#' A function to pull database information from a domino object
#'
#' @param dom A domino object that has been created
#' @param name_only A boolean for whether to return only the name of the database used
#'                  or the entire database that is stored. Default TRUE.
#' @return  A vector of unique databases used in building the domino object OR
#'          a data frame that includes the database information used in the domino object creation
#' @export
#' @examples
#' load("R/sysdata.rda")
#' database_name <- dom_database(pbmc_dom_built_tiny)
#' full_database <- dom_database(pbmc_dom_built_tiny, name_only = FALSE)
#' 
dom_database <- function(dom, name_only = TRUE) {
    db <- slot(dom, "db_info")
    if (name_only) {
        return(unique(db$database_name))
    } else {
        return(db)
    }
}

#' Access z-scores
#'
#' A function to pull z-scored expression from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @return  A matrix containing the z-scored gene expression values for each gene (row) by cell (column)
#' @export
#' @examples
#' load("R/sysdata.rda")
#' zscores <- dom_zscores(pbmc_dom_built_tiny)
#'
dom_zscores <- function(dom) {
    slot(dom, "z_scores")
}

#' Access counts
#'
#' A function to pull gene expression from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @return  A matrix containing the gene expression values for each gene (row) by cell (column)
#' @export
#' @examples
#' load("R/sysdata.rda")
#' counts <- dom_counts(pbmc_dom_built_tiny)
#' 
dom_counts <- function(dom) {
    as.matrix(slot(dom, "counts"))
}

#' Access clusters
#'
#' A function to pull cluster information from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @param labels A boolean for whether to return the cluster labels for each cell or the clusters used for inferring communication
#' @return  A vector containing either the names of the clusters used OR factors of the cluster label for each individual cell
#' @export
#' @examples
#' load("R/sysdata.rda")
#' cluster_names <- dom_clusters(pbmc_dom_built_tiny)
#' cell_cluster_label <- dom_clusters(pbmc_dom_built_tiny, labels = TRUE)
#' 
dom_clusters <- function(dom, labels = FALSE) {
    clust <- slot(dom, "clusters")
    if (labels) {
        return(clust)
    } else {
        return(levels(clust))
    }
}

#' Access transcription factor activation
#'
#' A function to pull transcription factor activation scores from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @return  A matrix containing the transcription factor activation scores for each feature (row) by cell (column)
#' @export
#' @examples
#' load("R/sysdata.rda")
#' tf_activation <- dom_tf_activation(pbmc_dom_built_tiny)
#' 
dom_tf_activation <- function(dom) {
    slot(dom, "features")
}

#' Access correlations
#'
#' A function to pull receptor-transcription factor correlations from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @return  A matrix containing the correlation values for each receptor (row) by transcription factor (column)
#' @export
#' @examples
#' load("R/sysdata.rda")
#' cor_matrix <- dom_correlations(pbmc_dom_built_tiny)
#' 
dom_correlations <- function(dom) {
    slot(dom, "cor")
}

#' Access linkages
#'
#' A function to pull linkages from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @param link_type One value (out of "complexes", "receptor-ligand",
#'                  "tf-target", "tf-receptor", "receptor", "incoming-ligand") used
#'                  to select the desired type of linkage
#' @param by_cluster A boolean to indicate whether the linkages should be returned overall or by cluster
#' @return  A list containing linkages between some combination of receptors, ligands, transcription factors, and clusters
#' @export
#' @examples
#' load("R/sysdata.rda")
#' complexes <- dom_linkages(pbmc_dom_built_tiny, "complexes")
#' tf_rec_by_cluster <- dom_linkages(pbmc_dom_built_tiny, "tf-receptor", TRUE)
dom_linkages <- function(dom, link_type = c(
                            "complexes", "receptor-ligand",
                            "tf-target", "tf-receptor", "receptor", "incoming-ligand"
                        ), by_cluster = FALSE) {
    links <- slot(dom, "linkages")
    if (by_cluster) {
        if (link_type == "tf-receptor") {
            return(links$clust_tf)
        } else if(link_type == "receptor") {
            return(links$clust_rec) 
        } else if(link_type == "incoming-ligand") {
            return(links$clust_incoming_lig)
        } else {
            stop("This linkage type is not available")
        }
    } else {
        if (link_type == "complexes") {
            return(links$complexes)
        } else if (link_type == "receptor-ligand") {
            return(links$rec_lig)
        } else if (link_type == "tf-target") {
            return(links$tf_targets)
        } else if (link_type == "tf-receptor") {
            return(links$tf_rec)
        } else {
            stop("This linkage type is not available.")
        }
    }
}

#' Access signaling
#'
#' A function to pull signaling matrices from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @param cluster Either NULL to indicate global signaling or a specific cluster for which a
#' @return  A data.frame containing the signaling score through each ligand (row) by each cluster (column) OR
#'          a data.frame containing the global summed signaling scores between receptors (rows) and ligands (columns) of each cluster
#' @export
#' @examples
#' load("R/sysdata.rda")
#' monocyte_signaling <- dom_signaling(pbmc_dom_built_tiny, cluster = "CD14_monocyte")
#' 
dom_signaling <- function(dom, cluster = NULL) {
    if (is.null(cluster)) {
        as.data.frame(slot(dom, "signaling"))
    } else {
        as.data.frame(slot(dom, "cl_signaling_matrices")[cluster])
    }
}

#' Access differential expression
#'
#' A function to pull differential expression p values from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @return  A matrix containing the p values for differential expression of transcription factors (rows) in each cluster (columns)
#' @export
#' @examples
#' load("R/sysdata.rda")
#' de_mat <- dom_de(pbmc_dom)
#' 
dom_de <- function(dom) {
    slot(dom, "clust_de")
}

#' Access build information
#'
#' A function to pull the parameters used when running [build_domino()] from a domino object
#'
#' @param dom A domino object that has been created with [create_domino()]
#' @return  A list containing booleans for whether the object has been created and build, and a list of the
#'          build parameters that were used in [build_domino()] to infer the signaling network
#' @export
#' @examples
#' load("R/sysdata.rda")
#' build_details <- dom_info(pbmc_dom_built_tiny)
#' 
dom_info <- function(dom) {
    info <- slot(dom, "misc")
    return(list(
        "create" = info$create, "build" = info$build,
        "build_variables" = info$build_vars
    ))
}

#' Access all features, receptors, or ligands present in a signaling network.
#'
#' This function collates all of the features, receptors, or ligands found in a
#' signaling network anywhere in a list of clusters. This can be useful for
#' comparing signaling networks across two separate conditions. In order to run
#' this [build_domino()] must be run on the object previously.
#'
#' @param dom Domino object containing a signaling network (i.e. [build_domino()] was run)
#' @param return String indicating where to collate "features", "receptors", or "ligands". If "all" then a list of all three will be returned.
#' @param clusters Vector indicating clusters to collate network items from. If left as NULL then all clusters will be included.
#' @return A vector containing all features, receptors, or ligands in the data set or a list containing all three.
#' @export
#' @examples
#' load("R/sysdata.rda")
#' monocyte_receptors <- collate_network_items(pbmc_dom_built_tiny, "CD14_monocyte", "receptors")
#' all_tfs <- collate_network_items(pbmc_dom_built_tiny, return = "features")
#' 
dom_network_items <- function(dom, clusters = NULL, return = NULL) {
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
