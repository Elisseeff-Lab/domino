#' Calculate a signaling network for a domino object
#'
#' This function calculates a signaling network. It requires a domino object
#' preprocessed from create_domino and returns a domino object prepared for
#' plotting with the various plotting functions in this package.
#'
#' @param dom Domino object from create_domino.
#' @param max_tf_per_clust Maximum number of transcription factors called active in a cluster.
#' @param min_tf_pval Minimum p-value from differential feature score test to call a transcription factor active in a cluster.
#' @param max_rec_per_tf Maximum number of receptors to link to each transcription factor.
#' @param rec_tf_cor_threshold Minimum pearson correlation used to consider a receptor linked with a transcription factor. Increasing this will decrease the number of receptors linked to each transcription factor.
#' @param min_rec_percentage Minimum percentage of cells in cluster expressing a receptor for the receptor to be linked to trancription factors in that cluster.
#' @return A domino object with a signaling network built
#' @export
#' @examples
#' pbmc_dom_tiny_built <- build_domino(
#'  dom = domino2:::pbmc_dom_tiny, min_tf_pval = .001, max_tf_per_clust = 25,
#'  max_rec_per_tf = 25, rec_tf_cor_threshold = .25, min_rec_percentage = 0.1
#' )
#' 
build_domino <- function(
    dom, max_tf_per_clust = 5, min_tf_pval = 0.01, max_rec_per_tf = 5, rec_tf_cor_threshold = 0.15,
    min_rec_percentage = 0.1) {
  if (dom@misc[["create"]] == FALSE) {
    stop("Please run domino_create to create the domino object.")
  }
  dom@misc[["build"]] <- TRUE
  dom@misc[["build_vars"]] <- c(
    max_tf_per_clust = max_tf_per_clust, min_tf_pval = min_tf_pval,
    max_rec_per_tf = max_rec_per_tf, rec_tf_cor_threshold = rec_tf_cor_threshold, min_rec_percentage = min_rec_percentage
  )
  if (length(dom@clusters)) {
    # Get transcription factors for each cluster
    clust_tf <- list()
    for (clust in levels(dom@clusters)) {
      ordered <- sort(dom@clust_de[, clust], decreasing = FALSE)
      # If there are more ties than max tf per cluster then rank by logfc
      if (sum(ordered == 0) > max_tf_per_clust) {
        zeros <- names(ordered)[which(ordered == 0)]
        fcs <- c()
        for (zero in zeros) {
          fc <- mean(dom@features[zero, which(dom@clusters == clust)]) - mean(dom@features[
            zero,
            which(dom@clusters != clust)
          ])
          fcs <- c(fcs, fc)
        }
        names(fcs) <- zeros
        sorted <- sort(fcs, decreasing = TRUE)[1:max_tf_per_clust]
      } else {
        sorted <- ordered[which(ordered < min_tf_pval)]
      }
      if (length(sorted) > max_tf_per_clust) {
        sorted <- sorted[1:max_tf_per_clust]
      }
      clust_tf[[clust]] <- names(sorted)
    }
    dom@linkages[["clust_tf"]] <- clust_tf
    # Get receptors for each transcription factor
    tf_rec <- list()
    for (tf in colnames(dom@cor)) {
      ordered <- sort(dom@cor[, tf], decreasing = TRUE)
      filtered <- ordered[which(ordered > rec_tf_cor_threshold)]
      if (length(filtered) > max_rec_per_tf) {
        top_receptors <- names(filtered)[1:max_rec_per_tf]
      } else {
        top_receptors <- names(filtered)
      }
      tf_rec[[tf]] <- top_receptors
    }
    dom@linkages[["tf_rec"]] <- tf_rec
    # If cluster methods are used, provide cluster-specific tf_rec linkages
    cl_tf_rec <- list()
    for (clust in levels(dom@clusters)) {
      percent <- dom@misc$cl_rec_percent[, clust]
      pass_genes <- names(percent[percent > min_rec_percentage])
      expressed <- c()
      for (rec in names(dom@linkages$rec_lig)) {
        if (rec %in% names(dom@linkages$complexes)) {
          rec_gene <- dom@linkages$complexes[[rec]]
        } else {
          rec_gene <- rec
        }
        if (length(rec_gene) == sum(rec_gene %in% pass_genes)) {
          expressed <- c(expressed, rec)
        }
      }
      active_tf <- dom@linkages$clust_tf[[clust]]
      cl_tf_rec[[clust]] <- lapply(dom@linkages$tf_rec[active_tf], FUN = function(x) {
        return(x[x %in% expressed])
      })
    }
    dom@linkages[["clust_tf_rec"]] <- cl_tf_rec
    # Get a list of active receptors for each cluster
    clust_rec <- list()
    for (clust in levels(dom@clusters)) {
      vec <- lc(dom@linkages$clust_tf_rec[[clust]], lc(clust_tf, clust))
      vec <- unique(vec[!is.na(vec)])
      clust_rec[[clust]] <- vec
    }
    dom@linkages[["clust_rec"]] <- clust_rec
    # Get a list of incoming ligands for each cluster
    clust_ligs <- list()
    for (clust in levels(dom@clusters)) {
      vec <- lc(dom@linkages$rec_lig, lc(clust_rec, clust))
      vec <- unique(vec[!is.na(vec)])
      clust_ligs[[clust]] <- vec
    }
    dom@linkages[["clust_incoming_lig"]] <- clust_ligs
    # Build signaling matrices for each cluster
    cl_signaling_matrices <- list()
    signaling <- matrix(0, ncol = length(levels(dom@clusters)), nrow = length(levels(dom@clusters)))
    rownames(signaling) <- paste0("R_", levels(dom@clusters))
    colnames(signaling) <- paste0("L_", levels(dom@clusters))
    for (clust in levels(dom@clusters)) {
      inc_ligs <- clust_ligs[[clust]]
      rl_map <- dom@misc[["rl_map"]]
      inc_ligs <- sapply(inc_ligs, function(l) {
        int <- rl_map[rl_map$L.name == l, ][1, ]
        if ((int$L.name != int$L.gene) & !grepl("\\,", int$L.gene)) {
          int$L.gene
        } else {
          int$L.name
        }
      })
      if (length(dom@linkages$complexes) > 0) {
        # if complexes were used
        inc_ligs_list <- lapply(inc_ligs, function(l) {
          if (l %in% names(dom@linkages$complexes)) {
            return(dom@linkages$complexes[[l]])
          } else {
            return(l)
          }
        })
        names(inc_ligs_list) <- inc_ligs
        inc_ligs <- unlist(inc_ligs_list)
      }
      lig_genes <- intersect(inc_ligs, rownames(dom@z_scores))
      if (length(lig_genes) %in% c(0, 1)) {
        lig_genes <- numeric(0)
      }
      cl_sig_mat <- matrix(0, ncol = length(levels(dom@clusters)), nrow = length(lig_genes))
      colnames(cl_sig_mat) <- colnames(signaling)
      rownames(cl_sig_mat) <- lig_genes
      for (c2 in levels(dom@clusters)) {
        n_cell <- length(which(dom@clusters == c2))
        if (n_cell > 1) {
          expr <- matrix(dom@z_scores[lig_genes, which(dom@clusters == c2)], nrow = length(lig_genes))
          sig <- rowMeans(expr)
        } else if (n_cell == 1) {
          sig <- dom@z_scores[lig_genes, which(dom@clusters == c2)]
        } else {
          sig <- rep(0, length(lig_genes))
          names(sig) <- lig_genes
        }
        # mean scaled expression less than 0 is brought up to 0 as a floor
        sig[which(sig < 0)] <- 0
        cl_sig_mat[, paste0("L_", c2)] <- sig
      }
      if (length(dom@linkages$complexes) > 0) {
        # if complexes were used
        cl_sig_list <- lapply(seq_along(inc_ligs_list), function(x) {
          if (all(inc_ligs_list[[x]] %in% lig_genes)) {
            # Some of the ligands in the list object may not be present in the data
            if (length(inc_ligs_list[[x]]) > 1) {
              return(colMeans(cl_sig_mat[inc_ligs_list[[x]], ]))
            } else {
              return(cl_sig_mat[inc_ligs_list[[x]], ])
            }
          }
        })
        names(cl_sig_list) <- names(inc_ligs_list)
        if (length(cl_sig_list) > 1) {
          cl_sig_mat <- do.call(rbind, cl_sig_list)
        }
      }
      cl_signaling_matrices[[clust]] <- cl_sig_mat
      if (nrow(cl_sig_mat) > 1) {
        signaling[paste0("R_", clust), ] <- colSums(cl_sig_mat)
      } else {
        signaling[paste0("R_", clust), ] <- 0
      }
    }
    dom@cl_signaling_matrices <- cl_signaling_matrices
    dom@signaling <- signaling
  } else {
    # If clusters are not defined, take all TFs selected previously
    dom@linkages[["clust_tf"]] <- list(clust = rownames(dom@features))
    # ID receptors for transcription factors
    tf_rec <- list()
    for (tf in colnames(dom@cor)) {
      ordered <- sort(dom@cor[, tf], decreasing = TRUE)
      filtered <- ordered[which(ordered > rec_tf_cor_threshold)]
      if (length(filtered) > max_rec_per_tf) {
        top_receptors <- names(filtered)[1:max_rec_per_tf]
      } else {
        top_receptors <- names(filtered)
      }
      tf_rec[[tf]] <- top_receptors
    }
    dom@linkages[["tf_rec"]] <- tf_rec
  }
  return(dom)
}

#' Pulls all items from a list pooled into a single vector
#'
#' Helper function to convert from a nested series of lists to a single vector.
#'
#' @param list List to pull items from
#' @param list_names Names of items in list to pool
#' @return A vector contaning all items in the list by list_names
#' @keywords internal
#'
lc <- function(list, list_names) {
  vec <- c()
  for (name in list_names) {
    vec <- c(vec, list[[name]])
  }
  return(vec)
}

#' Invert Receptor Ligand Data
#' @description
#' Reformat exisiting domino object to also record expression data/receptor mapping
#' from the ligand-receptor direction rather than receptor-ligand
#' @param dom domino object built by domino build
#' @return domino object with updated `rec_signaling` slot and `cr_signaling_matrices` slot
#' @export
#' 
invertRecLig <- function(dom) {
  
  ### Check
  if (!dom@misc[["build"]]) {
    stop("Must build a signaling network with domino_build prior to reformatting")
  }
  
  ### Get map
  rl_reading <- make_rl_reading(dom@db_info)
  
  ### Add linkages, if missing
  if (!"lig_rec" %in% names(dom@linkages)) {
    
    lig_genes <- unique(unlist(strsplit(rl_reading[["L.gene"]], split = "\\,")))
    lig_names <- rl_reading[["L.name"]]
    
    lig_rec_linkage <- list()
    for (lig in lig_names) {
      inter <- rl_reading[rl_reading[["L.name"]] == lig, ]
      recs <- inter[["R.name"]]
      lig_rec_linkage[[lig]] <- recs
    } # for lig
    
    dom@linkages[["lig_rec"]] <- lig_rec_linkage
    
  } # fi
  
  ### Grab object stuff
  clust_ligs <- dom@linkages$clust_incoming_lig
  clust_recs <- dom@linkages$clust_rec
  
  # Build signaling matrices for each cluster (This is identical to lines 111-185 in processing_fxns.R/build_domino, just switch all receptor-ligand references)
  cr_signaling_matrices <- list()
  signaling <- matrix(0, ncol = length(levels(dom@clusters)), nrow = length(levels(dom@clusters)))
  rownames(signaling) <- paste0("L_", levels(dom@clusters))
  colnames(signaling) <- paste0("R_", levels(dom@clusters))
  for (clust in levels(dom@clusters)) {
    inc_recs <- clust_recs[[clust]]
    rl_map <- dom@misc[["rl_map"]]
    inc_recs <- sapply(inc_recs, function(r) {
      int <- rl_map[rl_map$R.name == r, ][1, ]
      if ((int$R.name != int$R.gene) & !grepl("\\,", int$R.gene)) {
        int$R.gene
      } else {
        int$R.name
      }
    })
    if (length(dom@linkages$complexes) > 0) {
      # if complexes were used
      inc_recs_list <- lapply(inc_recs, function(r) {
        if (r %in% names(dom@linkages$complexes)) {
          return(dom@linkages$complexes[[r]])
        } else {
          return(r)
        }
      })
      names(inc_recs_list) <- inc_recs
      inc_recs <- unlist(inc_recs_list)
    }
    rec_genes <- intersect(inc_recs, rownames(dom@z_scores))
    if (length(rec_genes) %in% c(0, 1)) {
      rec_genes <- numeric(0)
    }
    cr_sig_mat <- matrix(0, ncol = length(levels(dom@clusters)), nrow = length(rec_genes))
    colnames(cr_sig_mat) <- colnames(signaling)
    rownames(cr_sig_mat) <- rec_genes
    for (c2 in levels(dom@clusters)) {
      n_cell <- length(which(dom@clusters == c2))
      if (n_cell > 1) {
        expr <- matrix(dom@z_scores[rec_genes, which(dom@clusters == c2)], nrow = length(rec_genes))
        sig <- rowMeans(expr)
      } else if (n_cell == 1) {
        sig <- dom@z_scores[rec_genes, which(dom@clusters == c2)]
      } else {
        sig <- rep(0, length(rec_genes))
        names(sig) <- rec_genes
      }
      # mean scaled expression less than 0 is brought up to 0 as a floor
      sig[which(sig < 0)] <- 0
      cr_sig_mat[, paste0("R_", c2)] <- sig
    }
    if (length(dom@linkages$complexes) > 0) {
      # if complexes were used
      cr_sig_list <- lapply(seq_along(inc_recs_list), function(x) {
        if (all(inc_recs_list[[x]] %in% rec_genes)) {
          # Some of the ligands in the list object may not be present in the data
          if (length(inc_recs_list[[x]]) > 1) {
            return(colMeans(cr_sig_mat[inc_recs_list[[x]], ]))
          } else {
            return(cr_sig_mat[inc_recs_list[[x]], ])
          }
        }
      })
      names(cr_sig_list) <- names(inc_recs_list)
      if (length(cr_sig_list) > 1) {
        cr_sig_mat <- do.call(rbind, cr_sig_list)
      }
    }
    cr_signaling_matrices[[clust]] <- cr_sig_mat
    if (nrow(cr_sig_mat) > 1) {
      signaling[paste0("L_", clust), ] <- colSums(cr_sig_mat)
    } else {
      signaling[paste0("L_", clust), ] <- 0
    }
  }
  dom@cr_signaling_matrices <- cr_signaling_matrices
  dom@rec_signaling <- signaling
  
  return(dom)
  
}
