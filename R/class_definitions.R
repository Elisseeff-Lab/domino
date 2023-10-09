#' The domino Class
#' 
#' The domino class contains all information necessary to calculate receptor
#' ligand signaling. It contains z scored expression, cluster, feature values,
#' and formatted receptor ligand databases as well as all calculated 
#' intermediates.
#' 
#' @slot db_info List of data sets from lr database.
#' @slot counts Raw count gene expression data
#' @slot z_scores Matrix of z-scored expression data with cells as columns
#' @slot clusters Named factor with cluster identity of each cell
#' @slot features Matrix of features to correlate receptor-ligand expression with. Cells are columns and features are rows.
#' @slot df List containing transcriptional targets by transcription factor.
#' @slot cor Correlation matrix of receptor expression to features.
#' @slot linkages List of lists containing info linking cluster->tf->rec->lig
#' @slot clust_de Data frame containing differential expression results for features by cluster.
#' @slot misc List of miscellaneous info pertaining to run parameters etc.
#' @slot cl_signaling_matrices Incoming signaling matrix for each cluster
#' @slot signaling Signaling matrix between all clusters.
#' 
#' @name domino-class
#' @rdname domino-class
#' @exportClass domino
#' 
domino <- setClass(Class="domino", slots=c(db_info="list", z_scores="matrix", counts="dgCMatrix",
  clusters="factor", features="matrix", cor="matrix", linkages="list", clust_de="matrix",
  misc="list", cl_signaling_matrices="list", signaling="matrix"))
#' The Domino linkage summary class
#' 
#' The linkage summary class contains linkages established in multiple domino
#' objects through gene regulatory network inference and reference to receptor-
#' ligand data bases. A data frame summarizing meta features that describe the
#' domino objects compared in the linkage summary facilitates comparisons of
#' established linkages and differential signaling interactions across categorical
#' sample covariates.
#' 
#' @slot subject_names unique names for each domino result included in the summary
#' @slot subject_meta data.frame with each row describing one subject and columns describing features of the subjects by which to draw comparisons of signaling networks
#' @slot subject_linkages nested list of linkages inferred for each subject. Lists are stored in a heirarchical structure of subject-cluster-linkage where linkages include transcription factors (tfs), linkages between transcription factors and receptors (tfs_rec), active receptors (rec), possible receptor-ligand interactions (rec_lig), and incoming ligands (incoming_lig)
#' 
#' @name linkage-summary-class
#' @rdname linkage-summary-class
#' @exportClass linkage_summary
#' 
linkage_summary <- setClass(Class="linkage_summary", slots=c(subject_names="factor", subject_meta="data.frame",
  subject_linkages="list"))
#' Print domino object
#' 
#' Prints a summary of a domino object
#' 
#' @param x Domino object
#' @keywords internal
setMethod("print", "domino", function(x, ...) {
  if (x@misc$build) {
    cat("A domino object of", length(x@clusters), "cells
                Contains signaling between",
      length(levels(x@clusters)), "clusters
                Built with a maximum of", as.integer(x@misc$build_vars["max_tf_per_clust"]),
      "TFs per cluster
                and a maximum of", as.integer(x@misc$build_vars["max_rec_per_tf"]),
      "receptors per TF\n")
  } else {
    cat(c("A domino object of", length(x@clusters), "cells\n", "A signaling network has not been built\n"),
      sep="")
  }
})
#' Show domino object information
#' 
#' Shows content overview of domino object
#' 
#' @param object Domino object
#' @keywords internal
setMethod("show", "domino", function(object) {
  if (object@misc$build) {
    cat(c("A domino object of ", length(object@clusters), " cells\n", "Built with signaling between ",
      length(levels(object@clusters)), " clusters\n"), sep = "")
  } else {
    cat(c("A domino object of", length(object@clusters), "cells\n", "A signaling network has not been built\n"),
      sep = "")
  }
})
