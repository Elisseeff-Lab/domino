#' Access database
#'
#' A function to pull database information from a domino object
#'
#' @param dom a domino object that has been created
#' @param name_only a boolean for whether to return only the name of the database used
#'                  or the entire database that is stored. Default TRUE.
#' @return  A vector of unique databases used in building the domino object OR
#'          a data frame that includes the database information used in the domino object creation
#' @export
#' @examples
#' example(build_domino)
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
#' @param dom a domino object that has been created with [create_domino()]
#' @return  A matrix containing the z-scored gene expression values for each gene (row) by cell (column)
#' @export
#' @examples
#' example(build_domino)
#' zscores <- dom_zscores(pbmc_dom_built_tiny)
#'
dom_zscores <- function(dom) {
    slot(dom, "z_scores")
}

#' Access counts
#'
#' A function to pull gene expression from a domino object
#'
#' @param dom a domino object that has been created with [create_domino()]
#' @return  A matrix containing the gene expression values for each gene (row) by cell (column)
#' @export
#' @examples
#' example(build_domino)
#' counts <- dom_counts(pbmc_dom_built_tiny)
#' 
dom_counts <- function(dom) {
    as.matrix(slot(dom, "counts"))
}

#' Access clusters
#'
#' A function to pull cluster information from a domino object
#'
#' @param dom a domino object that has been created with [create_domino()]
#' @param labels a boolean for whether to return the cluster labels for each cell or the clusters used for inferring communication
#' @return  A vector containing either the names of the clusters used OR factors of the cluster label for each individual cell
#' @export
#' @examples
#' example(build_domino)
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
#' @param dom a domino object that has been created with [create_domino()]
#' @return  A matrix containing the transcription factor activation scores for each TF (row) by cell (column)
#' @export
#' @examples
#' example(build_domino)
#' tf_activation <- dom_tf_activation(pbmc_dom_built_tiny)
#' 
dom_tf_activation <- function(dom) {
    slot(dom, "features")
}

#' Access correlations
#'
#' A function to pull receptor-transcription factor correlations from a domino object
#'
#' @param dom a domino object that has been created with [create_domino()]
#' @param type either "rl" or "complex", to select between the receptor-ligand or complex correlation matrix
#' @return  A matrix containing the correlation values for each receptor (row) by transcription factor (column)
#' @export
#' @examples
#' example(build_domino)
#' cor_matrix <- dom_correlations(pbmc_dom_built_tiny, "rl")
#' 
dom_correlations <- function(dom, type = "rl") {
    if (type == "complex") {
        corrs <- slot(dom, "cor")
    } else if (type == "rl") {
        misc <- slot(dom, "misc")
        corrs <- misc$rec_cor
    } else {
        stop("Type must be either 'rl' or 'complex'")
    }
    return(corrs)
}

#' Access linkages
#'
#' A function to pull linkages from a domino object
#'
#' @param dom a domino object that has been created with [create_domino()]
#' @param link_type one value (out of "complexes", "receptor-ligand",
#'                  "tf-target", "tf-receptor", "receptor", "incoming-ligand") used
#'                  to select the desired type of linkage
#' @param by_cluster a boolean to indicate whether the linkages should be returned overall or by cluster
#' @return  A list containing linkages between some combination of receptors, ligands, transcription factors, and clusters
#' @export
#' @examples
#' example(build_domino)
#' complexes <- dom_linkages(pbmc_dom_built_tiny, "complexes")
#' tf_rec_by_cluster <- dom_linkages(pbmc_dom_built_tiny, "tf-receptor", TRUE)
#' 
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
#' @param dom a domino object that has been created with [create_domino()]
#' @param cluster either NULL to indicate global signaling or a specific cluster for which a signaling matrix is desired
#' @return  A data frame containing the signaling score through each ligand (row) by each cluster (column) OR
#'          a data frame containing the global summed signaling scores between receptors (rows) and ligands (columns) of each cluster
#' @export
#' @examples
#' example(build_domino)
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
#' A function to pull differential expression p-values from a domino object
#'
#' @param dom a domino object that has been created with [create_domino()]
#' @return  A matrix containing the p-values for differential expression of transcription factors (rows) in each cluster (columns)
#' @export
#' @examples
#' example(build_domino)
#' de_mat <- dom_de(pbmc_dom_built_tiny)
#' 
dom_de <- function(dom) {
    slot(dom, "clust_de")
}

#' Access build information
#'
#' A function to pull the parameters used when running [build_domino()] from a domino object
#'
#' @param dom a domino object that has been created with [create_domino()]
#' @return  A list containing booleans for whether the object has been created and built and a list of the
#'          build parameters that were used in [build_domino()] to infer the signaling network
#' @export
#' @examples
#' example(build_domino)
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
#' @param dom a domino object containing a signaling network (i.e. [build_domino()] was run)
#' @param return string indicating whether to collate "features", "receptors", or "ligands". If "all" then a list of all three will be returned.
#' @param clusters vector indicating clusters to collate network items from. If left as NULL then all clusters will be included.
#' @return A vector containing all features, receptors, or ligands in the data set or a list containing all three.
#' @export
#' @examples
#' example(build_domino)
#' monocyte_receptors <- dom_network_items(pbmc_dom_built_tiny, "CD14_monocyte", "receptors")
#' all_tfs <- dom_network_items(pbmc_dom_built_tiny, return = "features")
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
        tf_wo_rec <- which(
          vapply(dom@linkages$clust_tf_rec[[cl]], FUN.VALUE = numeric(1), FUN = length) == 0
        )
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

#' Check input arguments
#'
#' Accepts an object and rules to check against; stops if requirements are not met
#'
#' @param arg the argument to check
#' @param allow_class vector of allowed classes
#' @param allow_len vector of allowed lengths
#' @param allow_range range of minimum and maximum values i.e. c(1, 5)
#' @param allow_values vector of allowed values
#' @param need_vars vector of required variables
#' @param need_colnames vogical for whether colnames are required
#' @param need_rownames logical for whether rownames are required
#' @param need_names logical for whether names are required
#' @return Logical indicating whether the argument meets the requirements
#' @keywords internal
#' 
check_arg <- function(arg, allow_class = NULL, allow_len = NULL,
                      allow_range = NULL, allow_values = NULL,
                      need_vars = c(NULL), need_colnames = FALSE,
                      need_rownames = FALSE, need_names = FALSE) {
  argname <- deparse(substitute(arg))
  classes <- paste(allow_class, collapse = ",")
  lengths <- paste(allow_len, collapse = ",")

  if (!is.null(allow_class)) {
    if (!any((class(arg) %in% allow_class))) {
      stop(sprintf("Class of %s must be one of: %s", argname, classes))
    }
  }

  if (!is.null(allow_len)) {
    if (!(length(arg) %in% allow_len)) {
      stop(sprintf("Length of %s must be one of: %s", argname, lengths))
    }
  }

  if (!is.null(need_vars)) {
    if (!all(need_vars %in% names(arg))) {
      stop(sprintf("Required variables %s not found in %s",
                   paste0(need_vars, collapse = ", "), argname))
    }
  }

  if (need_rownames) {
    if (is.null(rownames(arg))) {
      stop(sprintf("No rownames found in %s", argname))
    }
  }

  if (need_colnames) {
    if (is.null(colnames(arg))) {
      stop(sprintf("No colnames found in %s", argname))
    }
  }

  if (need_names) {
    if (is.null(names(arg))) {
      stop(sprintf("No names found in %s", argname))
    }
  }

  if (!is.null(allow_range)) {
    if (all(arg < allow_range[1]) || all(arg > allow_range[2])) {
      stop(sprintf("All values in %s must be between %s and %s",
                   argname, allow_range[1], allow_range[2]))
    }
  }

  if (!is.null(allow_values)) {
    if (!all(arg %in% allow_values)) {
      stop(sprintf("All values in %s must be one of: %s",
                   argname, paste(allow_values, collapse = ", ")))
    }
  }

}

#' Read in data if an object looks like path to it
#'
#' @param obj object to read if not already object
#' @return Object itself or data read in from path
#' @keywords internal
read_if_char <- function(obj) {
  if (is(obj, "character")) {
    check_arg(obj, allow_class = "character", allow_len = 1)
    obj <- read.csv(obj, stringsAsFactors = FALSE)
  }
  return(obj)
}

#' Change cases of True/False syntax from Python to TRUE/FALSE R syntax
#'
#' @param obj object that will be converted
#' @return The converted object
#' @keywords internal
conv_py_bools <- function(obj) {
  for (x in colnames(obj)) {
    bools <- sort(unique(obj[[x]]))
    if (identical(bools, c("False", "True"))) {
      obj[[x]] <- ifelse(obj[[x]] == "True", TRUE, FALSE)
    }
  }
  return(obj)
}


#' Create a mock linkage summary object
#'
#' @return obj a linkage summary object
#' @export
mock_linkage_summary <- function() {
  linkage_sum_tiny <- new("linkage_summary",
  subject_meta = data.frame(
    "subject_names" = paste0("P",1:6),
    "group" = c(rep("G1", 3), rep("G2", 3))
  ), 
  subject_names = factor(
    paste0("P",1:6), levels = paste0("P",1:6)
  ), 
  subject_linkages = list(
    "P1" = list(
      "C1" = list(
        "tfs" = c("TF1", "TF2", "TF3", "TF4"),
        "rec" = c("R1", "R2", "R3", "R4"),
        "incoming_lig" = c("L1", "L2", "L3", "L4"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
      ),
      "C2" = list(
        "tfs" = c("TF2", "TF3", "TF4"),
        "rec" = c("R2", "R3", "R4"),
        "incoming_lig" = c("L2", "L3", "L4"),
        "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
      )
    ),
    "P2" = list(
      "C1" = list(
        "tfs" = c("TF1", "TF2"),
        "rec" = c("R1", "R2"),
        "incoming_lig" = c("L1", "L2"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2")
      ),
      "C2" = list(
        "tfs" = c("TF3", "TF4"),
        "rec" = c("R3", "R4"),
        "incoming_lig" = c("L3", "L4"),
        "tfs_rec" = c("TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R3 <- L3", "R4 <- L4")
      )
    ),
    "P3" = list(
      "C1" = list(
        "tfs" = c("TF1", "TF2"),
        "rec" = c("R1", "R2"),
        "incoming_lig" = c("L1", "L2"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2")
      ),
      "C2" = list(
        "tfs" = c("TF3"),
        "rec" = c("R3"),
        "incoming_lig" = c("L3"),
        "tfs_rec" = c("TF3 <- R3"),
        "rec_lig" = c("R3 <- L3")
      )
    ),
    "P4" = list(
      "C1" = list(
        "tfs" = c("TF2", "TF3", "TF4"),
        "rec" = c("R2", "R3", "R4"),
        "incoming_lig" = c("L2", "L3", "L4"),
        "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
      ),
      "C2" = list(
        "tfs" = c("TF1", "TF2", "TF3", "TF4"),
        "rec" = c("R1", "R2", "R3", "R4"),
        "incoming_lig" = c("L1", "L2", "L3", "L4"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
      )
    ),
    "P5" = list(
      "C1" = list(
        "tfs" = c("TF3"),
        "rec" = c("R3"),
        "incoming_lig" = c("L3"),
        "tfs_rec" = c("TF3 <- R3"),
        "rec_lig" = c("R3 <- L3")
      ),
      "C2" = list(
        "tfs" = c("TF1", "TF2", "TF3", "TF4"),
        "rec" = c("R1", "R2", "R3", "R4"),
        "incoming_lig" = c("L1", "L2", "L3", "L4"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
      )
    ),
    "P6" = list(
      "C1" = list(
        "tfs" = c(),
        "rec" = c(),
        "incoming_lig" = c(),
        "tfs_rec" = c(),
        "rec_lig" = c()
      ),
      "C2" = list(
        "tfs" = c("TF1", "TF2", "TF3"),
        "rec" = c("R1", "R2", "R3"),
        "incoming_lig" = c("L1", "L2", "L3"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3")
      )
    )
  )
)
return(linkage_sum_tiny)
}
