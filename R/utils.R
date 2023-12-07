#' Check input arguments
#'
#' Accepts an object and rules for checking, stops if rules not met.
#'
#' @param arg The argument to check
#' @param allow_class Vector of allowed classes
#' @param allow_len Vector of allowed lengths
check_arg <- function(arg, allow_class = c("character"), allow_len = NULL) {
  argname <- deparse(substitute(arg))
  classes <- paste(allow_class, collapse = ",")
  lengths <- paste(allow_len, collapse = ",")

  if (!(class(arg) %in% allow_class)) {
    stop(sprintf("Class of %s must be one of: %s", argname, classes))
  }

  if (!is.null(allow_len)) {
    if (!(length(arg) %in% allow_len)) {
      stop(sprintf("Length of %s must be one of: %s", argname, lengths))
    }
  }
}

#' Read in data if an object looks like path to it.
#'
#' @param obj Object to read if not already object
#' @return obj Object itself in case its not a character
read_if_char <- function(obj) {
  if (is(obj, "character")) {
    obj <- read.csv(obj, stringsAsFactors = FALSE)
  }
  return(obj)
}


#' Change cases of True/False syntax from Python to TRUE/FALSE R syntax
#'
#' @param obj Object that will be converted
#' @return obj The converted object
conv_py_bools <- function(obj) {
  for (x in colnames(obj)) {
    bools <- sort(unique(obj[[x]]))
    if (identical(bools, c("False", "True"))) {
      obj[[x]] <- ifelse(obj[[x]] == "True", TRUE, FALSE)
    }
  }
  return(obj)
}

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

dom_database <- function(dom, name_only = TRUE) {
    db = slot(dom, "db_info")
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

dom_clusters <- function(dom, labels = FALSE) {
    clust = slot(dom, "clusters")
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
dom_linkages <- function(dom, link_type = c("complexes", "receptor-ligand",
                "tf-target", "tf-receptor", "receptor", "incoming-ligand"), by_cluster = FALSE) {
    links = slot(dom, "linkages")
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
        } else if (link_type == "tf_targets") {
            return(links$tf_targets)
        } else if (link_type == "tf_receptor") {
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
dom_signaling <- function(dom, cluster = NULL) {
    if(is.null(cluster)) {
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
dom_info <- function(dom) {
    info = slot(dom, "misc")
    return(list("create" = info$create, "build"= info$build, 
        "build_variables"= info$build_vars))
}

