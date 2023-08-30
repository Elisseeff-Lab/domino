#' Summarize linkages from multiple domino objects
#' 
#' @param domino_results list of domino result with one domino object per subject. Names from the list must match subject_names.
#' @param subject_meta dataframe that includes the subject features by which the objects could be grouped. The first column should must be subject names
#' @param subject_names vector of subject names in domino_results. If NULL, defaults to first column of subject_meta.
#' @return A linkage summary class object consisting of nested lists of the active transcription factors, active receptors, and incoming ligands for each cluster across multiple domino results.
#' @export
#' 
summarize_linkages = function(domino_results, subject_meta, subject_names = NULL){
  if(!is(domino_results, "list")){
    stop("domino_results must be provided as a named list where names correspond to subject names")
  }
  if(is.null(subject_names)){
    subject_names = subject_meta[,1]
  }
  if(sum(subject_names %in% names(domino_results)) == 0){
    stop("No provided subject names match names from the domino results list")
  }
  if(sum(!subject_names %in% names(domino_results))){
    extra_names = subject_names[!subject_names %in% names(domino_results)]
    warning(paste0("Provided subject names included names not present in domino_results: ",
                   paste(extra_names, collapse = ", ")))
    subject_names = subject_names[subject_names %in% names(domino_results)]
  }
  if(length(subject_names) < length(names(domino_results))){
    warning(paste0("Linkage summary includes results only for provided subject names: ",
                   paste(subject_names, collapse = ", ")))
    subject_meta = subject_meta[subject_meta[,1] %in% subject_names,]
  }
  
  subject_linkages = list()
  for(id in subject_names){
    dom <- domino_results[[id]]
    clusters <- levels(dom@clusters)
    
    c_features <- list()
    for(cluster in clusters){
      # list of t.factors active in cell type
      tfs <- dom@linkages$clust_tf[[cluster]]
      
      # obtain all receptors linked to active t. factors
      rec <- dom@linkages$clust_rec[[cluster]]
      
      # limit to unique entries
      tfs <- unique(tfs)
      rec <- unique(rec)
      
      # obtain all incoming ligands that interact with the receptors
      # limited to those present in data set
      lig <- rownames(dom@cl_signaling_matrices[[cluster]])
      
      # linkages of t.factors and receptors per cluster
      tfs_rec <- c()
      for(t in tfs){
        for(r in dom@linkages$clust_tf_rec[[cluster]][[t]]){
          tfs_rec <- c(tfs_rec, t, r)
        }
      }
      
      rec_lig <- c()
      for(r in rec){
        for(l in dom@linkages$rec_lig[[r]]){
          if(!l %in% lig){next}
          rec_lig <- c(rec_lig, r, l)
        }
      }
      
      # stitch linked t.factors-receptors, receptors-ligands
      int_tfs_rec <- c()
      int_rec_lig <- c()
      for(i in 0:((length(tfs_rec)/2)-1)){
        # count by twos and paste together with a <- denoting direction
        s <- i * 2
        interact <- paste(tfs_rec[1 + s], tfs_rec[2 + s], sep = " <- ")
        int_tfs_rec <- c(int_tfs_rec, interact)
      }
      for(i in 0:((length(rec_lig)/2)-1)){
        # count by twos and paste together with a "<-" denoting direction
        s <- i * 2
        interact <- paste(rec_lig[1 + s], rec_lig[2 + s], sep = " <- ")
        int_rec_lig <- c(int_rec_lig, interact)
      }
      
      # save the features of this cluster
      c_features[[cluster]] <- list(
        "tfs" = tfs,
        "rec" = rec,
        "incoming_lig" = lig,
        "tfs_rec" = int_tfs_rec,
        "rec_lig" = int_rec_lig
      )
    }
    subject_linkages[[id]] <- c_features
  }
  return(
    linkage_summary(
      subjects = factor(subject_names),
      subject_meta = subject_meta,
      subject_linkages = subject_linkages
    )
  )
}
