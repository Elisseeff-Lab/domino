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
domino <- setClass(
    Class = 'domino',
    slots = c(
        db_info = 'list',
        z_scores = 'matrix',
        counts = 'dgCMatrix',
        clusters = 'factor',
        features = 'matrix',
        cor = 'matrix',
        linkages = 'list',
        clust_de = 'matrix',
        misc = 'list',
        cl_signaling_matrices = 'list',
        signaling = 'matrix'
    )
)

#' Print domino object
#' 
#' Prints a summary of a domino object
#' 
#' @param x Domino object
#' @S4method Print domino object
print.domino = function(x, ...){
    cat('A domino object of', length(x@clusters), 'cells
    Contains signaling between', length(levels(x@clusters)), 'clusters
    Built with a maximum of', as.integer(x@misc$build_vars["max_tf_per_clust"]), 'TFs per cluster
    and a maximum of', as.integer(x@misc$build_vars["max_rec_per_tf"]), 'receptors per TF\n')
}

#' Show domino object information
#' 
#' Shows content overview of domino object
#' 
#' @param object Domino object
#' @S4method Show domino information
setMethod("show", "domino",
    function(object){
        if (object@misc$build) {
            cat(c('A domino object of ', length(object@clusters), ' cells\n',
            'Built with signaling between ',
            length(levels(object@clusters)), ' clusters\n'), sep = '')
        } else {
            cat(c('A domino object of', length(object@clusters), 'cells\n',
            'A signaling network has not been built\n'), sep = '')
        }
    }
) 