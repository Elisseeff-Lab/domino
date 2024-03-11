# internal scripts for the create_domino() function

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
