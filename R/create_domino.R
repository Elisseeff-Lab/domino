# internal scripts for the create_domino() function

select_cluster_tf <- function(features, clusters, verbose = TRUE, method = "one.sided.wilcox") {
  cluster_lvls <- levels(clusters)
  n_lvls <- length(cluster_lvls)
  feat_names <- rownames(features)
  names(feat_names) <- feat_names
  if(verbose) {message("Calculating feature enrichment by cluster")}
  cl_pval_list <- lapply(
    cluster_lvls,
    function(cl) {
      a <- which(cluster_lvls == cl)
      if(verbose) {message(paste0(a, " of ", n_lvls))}
      cl_cells <- which(clusters == cl)
      out_cells <- which(clusters != cl)
      p_val_ls <- lapply(
        feat_names,
        function(f) {
          cl_value <- features[f, cl_cells]
          out_value <- features[f, out_cells]
          if(method == "one.sided.wilcox") {
            test_res <- stats::wilcox.test(
              x = cl_value, y = out_value, alternative = "greater"
            )
            return(test_res[["p.value"]])
          }
        }
      )
      p_vals <- unlist(p_val_ls)
      return(p_vals)
    }
  )
  pval_mat <- do.call(cbind, cl_pval_list)
  colnames(pval_mat) <- cluster_lvls
  return(pval_mat)
}
