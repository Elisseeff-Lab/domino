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
      subject_names = factor(subject_names),
      subject_meta = subject_meta,
      subject_linkages = subject_linkages
    )
  )
}

#' Count occurrences of linkages across multiple domino results from a linkage summary
#' 
#' @param linkage_summary a linkage_summary object
#' @param cluster the name of the cell cluster being compared across multiple domino results
#' @param group.by the name of the column in linkage_summary\@subject_meta by which to group subjects for counting. If NULL, only total counts of linkages for linkages in the cluster across all subjects is given.
#' @param linkage a stored linkage from the domino object. Can compare ("tfs", "rec", "incoming_lig", "tfs_rec", "rec_lig")
#' @param subject_names a vector of subject_names from the linkage_summary to be compared. If NULL, all subject_names in the linkage summary are included in counting.
#' @return a data frame with columns for the unique linkage features and the counts of how many times the linkage occured across the compared domino results. If group.by is used, counts of the linkages are also provided as columns named by the unique values of the group.by variable.
#' @export
#' 
count_linkage <- function(linkage_summary, cluster, group.by = NULL, 
                          linkage = "rec_lig", subject_names = NULL){
  if(is.null(subject_names)){
    subject_names = linkage_summary@subject_names
  }
  all_int <- sapply(
    linkage_summary@subject_linkages, 
    FUN = function(x){return(x[[cluster]][[linkage]])}
  )
  feature <- table(unlist(all_int))
  df <- data.frame(
    feature = names(feature),
    total_count = as.numeric(feature)
  )
  
  if(!is.null(group.by)){
    if(!group.by %in% colnames(linkage_summary@subject_meta)){
      stop("group.by variable not present in subject_meta")
    }
    groups <- levels(factor(linkage_summary@subject_meta[[group.by]]))
    
    for(g in groups){
      g_index <- linkage_summary@subject_meta[[group.by]] == g
      g_subjects <- linkage_summary@subject_meta[g_index, 1]
      
      int_count <- list()
      for(f in df[["feature"]]){
        count<- sapply(
          g_subjects, 
          function(x){f %in% linkage_summary@subject_linkages[[x]][[cluster]][[linkage]]}
        )
        int_count[[f]] <- sum(count)
      }
      df[[g]] <- int_count
    }
  }
  return(df)
}

#' Statistical test for differential linkages across multiple domino results
#' 
#' @param linkage_summary a linkage_summary object
#' @param cluster the name of the cell cluster being compared across multiple domino results
#' @param group.by the name of the column in linkage_summary\@subject_meta by which to group subjects for counting.
#' @param linkage a stored linkage from the domino object. Can compare ("tfs", "rec", "incoming_lig", "tfs_rec", "rec_lig")
#' @param subject_names a vector of subject_names from the linkage_summary to be compared. If NULL, all subject_names in the linkage summary are included in counting.
#' @param test_name the statistical test used for comparison.
#' \itemize{
#'  \item{"fishers.exact"} : Fisher's exact test for the dependence of the proportion of subjects with an active linkage in the cluster on which group the subject belongs to in the group.by variable. Provides an odds ratio, p-value, and a Benjamini-Hochberg FDR-adjusted p-value (p.adj) for each linkage tested.
#' }
#' @return a data frame of results from the test of the differential linkages. Rows correspond to each linkage tested. Columns correspond to:
#' \itemize{
#'  \item{"cluster"} : the name of the cell cluster being compared
#'  \item{"linkage"} : the type of linkage being compared
#'  \item{"group.by"} : the grouping variable
#'  \item{"test_name"} : the test used for comparison
#'  \item{"feature"} : individual linkages compared
#'  \item{"test statistics"} : test statistics provided are based on test method. "fishers.exact" provides a odds ratio, p-value, and fdr-adjusted p-value.
#'  \item{"total_count"} : total number of subjects where the linkage is active
#'  \item{"*_count"} : number of subjects in each category of group.by where the linkage is active
#' }
#' @export
#' 
test_differential_linkages <- function(linkage_summary, cluster, group.by, 
                                       linkage = "rec_lig", subject_names = NULL,
                                       test_name = "fishers.exact"){
  valid_tests <- c("fishers.exact")
  if(!test_name %in% valid_tests){
    stop("test_name invalid")
  }
  if(is.null(subject_names)){
    subject_names = linkage_summary@subject_names
  }
  
  # count the number of groups
  subject_count <- as.data.frame(table(linkage_summary@subject_meta[[group.by]]))
  colnames(subject_count) <- c(group.by, "total")
  group_levels <- subject_count[[group.by]]
  level_n <- nrow(subject_count)
  count_link <- count_linkage(linkage_summary = linkage_summary,
                              cluster = cluster, linkage = linkage,
                              group.by = group.by, subject_names = subject_names)
  # initiate data frame for storing results
  n <- nrow(count_link)
  result_df <- data.frame(
    cluster = rep(cluster, n),
    linkage = rep(linkage, n),
    group.by = rep(group.by, n),
    test_name = rep(test_name, n),
    feature = count_link[["feature"]]
  )
  # empty contigency table
  test_mat <- matrix(data = NA, nrow = nrow(subject_count), ncol = 2)
  rownames(test_mat) <- subject_count[[group.by]]
  colnames(test_mat) <- c("linkage_present", "linkage_absent")
  test_template <- as.data.frame(test_mat)
  
  if(test_name == "fishers.exact"){
    test_result <- 
      as.data.frame(t(
        sapply(
          result_df[["feature"]],
          FUN = function(x){
            feat_count = count_link[count_link[["feature"]] == x, !colnames(count_link) %in% c("feature", "total_count")]
            g_names = colnames(feat_count)
            feat_count = sapply(feat_count, as.numeric)
            # fill contingency table
            test_df <- test_template
            test_df[["linkage_present"]] = feat_count
            test_df[["linkage_absent"]] = subject_count[["total"]] - feat_count
            # conduct test
            test = fisher.test(test_df)
            odds.ratio = test$estimate
            p.value = test$p.value
            res = c(odds.ratio, p.value)
            res = setNames(res, c("odds.ratio", "p.value"))
            return(res)
          }
        )
      ))
    # include fdr-adjusted p-values
    test_result[["p.adj"]] = p.adjust(
      p = test_result[["p.value"]], method = "fdr"
    )
  }
  # append test result columns to result
  result_df = cbind(result_df, test_result)
  
  # append counts of active linkages for each group
  count_append = count_link[,!colnames(count_link) == "feature"]
  colnames(count_append) <- sapply(colnames(count_append), function(x){if(!grepl("_count$", x)){return(paste0(x, "_count"))}else{return(x)}})
  result_df = cbind(result_df, count_append)
  
  return(result_df)
}
