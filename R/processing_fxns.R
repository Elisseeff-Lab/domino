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
#'  dom = dominoSignal:::pbmc_dom_tiny, min_tf_pval = .001, max_tf_per_clust = 25,
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
        sorted <- sort(fcs, decreasing = TRUE)[seq(max_tf_per_clust)]
      } else {
        sorted <- ordered[which(ordered < min_tf_pval)]
      }
      if (length(sorted) > max_tf_per_clust) {
        sorted <- sorted[seq(max_tf_per_clust)]
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
        top_receptors <- names(filtered)[seq(max_rec_per_tf)]
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
      inc_ligs <- vapply(inc_ligs, FUN.VALUE = character(1), FUN = function(l) {
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
        top_receptors <- names(filtered)[seq(max_rec_per_tf)]
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
