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

write_rec_lig_linkages <- function(rl_parse) {
  rec_names <- unique(rl_parse[["R.name"]])
  names(rec_names) <- rec_names
  rec_lig_ls <- lapply(
    rec_names,
    function(r) {
      rl_int_row <- rl_parse[rl_parse[["R.name"]] == r,]
      ligs <- rl_int_row[["L.name"]]
      return(ligs)
    }
  )
  return(rec_lig_ls)
}

read_rl_map_complexes <- function(rl_parse, use_complexes) {
  if(use_complexes) {
    complex_lig <- rl_parse[grepl("\\,", rl_parse$L.gene), "L.name"]
    complex_rec <- rl_parse[grepl("\\,", rl_parse$R.gene), "R.name"]
    complex_names <- unique(c(complex_lig, complex_rec))
    names(complex_names) <- complex_names
    gene_complexes <- lapply(
      complex_names, FUN = function(x) {
        comp_cs <- c(
          rl_parse[rl_parse$L.name == x, "L.gene"],
          rl_parse[rl_parse$R.name == x, "R.gene"]
        )
        if(length(unique(comp_cs)) > 1) {
          comp_first <- comp_cs[1]
          warning(
            paste0(
              "protein complex ", x, 
              " has multiple, differing descriptions of component genes. \n",
              "Defaulting to first encountered component description: \n",
              comp_first
            )
          )
          comp <- unique(unlist(strsplit(comp_first, split = "\\,")))
          return(comp)
        } else {
          comp <- unique(unlist(strsplit(comp_cs, split = "\\,")))
          return(comp)
        }
      }
    )
    return(gene_complexes)
  } else {
    return(NULL)
  }
}

read_rl_map_genes <- function(rl_map) {
  rl_list <- apply(
    rl_map, MARGIN = 1,
    FUN = function(x) {
      rl <- c()
      p <- ifelse(x[["type_A"]] == "R", "A", "B")
      q <- ifelse(p == "A", "B", "A")
      R.gene <- x[[paste0("gene_", p)]]
      L.gene <- x[[paste0("gene_", q)]]
      rl[["R.gene"]] <- R.gene
      rl[["L.gene"]] <- L.gene
      if (paste0("uniprot_", p) %in% names(x)) {
        rl[["R.uniprot"]] <- x[[paste0("uniprot_", p)]]
      }
      if (paste0("uniprot_", q) %in% names(x)) {
        rl[["L.uniprot"]] <- x[[paste0("uniprot_", q)]]
      }
      if (paste0("name_", p) %in% names(x)) {
        rl[["R.name"]] <- x[[paste0("name_", p)]]
      }
      if (paste0("name_", q) %in% names(x)) {
        rl[["L.name"]] <- x[[paste0("name_", q)]]
      }
      return(as.data.frame(rl))
    }
  )
  rl_df <- as.data.frame(do.call(rbind, rl_list))
  return(rl_df)
}
