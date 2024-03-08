# internal scripts for the create_domino() function

test_tfs_rec_linkage <- function(
    features, z_scores, counts, 
    feature_de,
    receptors, 
    method = "spearman.correlation", verbose = TRUE) {
  tfs <- colnames(feature_de)
  names(tfs) <- tfs
  n_tfs <- length(tfs)
  
  rec_z_scores <- z_scores[rownames(z_scores) %in% receptors,]
  rec_counts <- counts[rownames(counts) %in% receptors,]
  confirmed_recs <- rownames(rec_z_scores)
  
  if(verbose) {message("Calculating feature-receptor linkages")}
  cl_link_list <- lapply(
    tfs,
    function(tf) {
      a <- which(tf == tfs)
      tf_scores <- features[tf,]
      
      if(verbose) {message(paste0(a, " of ", n_tfs))}
      linkage_score_ls <- lapply(
        confirmed_recs,
        function(r) {
          r_exp <- rec_z_scores[r,]
          if(method == "spearman.correlation") {
            test_res <- stats::cor.test(rexp, tf_scores, method = "spearman", alternative = "greater")
            return(test_res[["estimate"]])
          }
        }
      )
      linkage_score <- unlist(linkage_score_ls)
      return(linkage_score)
    }
  )
  linkage_score_mat <- do.call(cbind, cl_link_list)
  colnames(linkage_score_mat) <- tfs
  return(linkage_score_mat)
}

filter_tf_regulon_receptors <- function(cor_mat, tf_targets) {
  
}

