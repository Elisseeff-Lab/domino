#' Renames clusters in a domino object
#' 
#' This function reads in a receptor ligand signaling database, cell level 
#' features of some kind (ie. output from pySCENIC), z-scored single cell data, 
#' and cluster id for single cell data, calculates a correlation matrix between 
#' receptors and other features (this is transcription factor module scores if 
#' using pySCENIC), and finds features enriched by cluster. It will return a 
#' domino object prepared for build_domino, which will calculate a signaling 
#' network.
#' 
#' @param dom A domino object to rename clusters in
#' @param clust_conv A named vector of conversions from old to new clusters. Values are taken as new clusters IDs and names as old cluster IDs.
#' @return A domino object with clusters renamed in all applicable slots.
#' @export 
#' 
rename_clusters = function(dom, clust_conv){
    if(dom@misc$create){
        dom@clusters = plyr::revalue(dom@clusters, clust_conv)
        colnames(dom@clust_de) = clust_conv
        names(colnames(dom@clust_de)) = c()
    }
    if(dom@misc$build){
        names(dom@linkages$clust_tf) = clust_conv
        colnames(dom@signaling) = paste0('L_', clust_conv)
        rownames(dom@signaling) = paste0('R_', clust_conv)
        names(dom@cl_signaling_matrices) = clust_conv
        for(cl in clust_conv){
            colnames(dom@cl_signaling_matrices[[cl]]) = paste0('L_', clust_conv)
        }
    }
    return(dom)
}