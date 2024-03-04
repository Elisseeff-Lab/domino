# internal scripts for the create_domino() function

parse_rl_map_genes <- function(rl_map) {
  rl_list <- apply(
    rl_map, MARGIN = 1,
    FUN = function(inter) {
      rl <- c()
      p <- ifelse(inter[["type_A"]] == "R", "A", "B")
      q <- ifelse(p == "A", "B", "A")
      R.gene <- inter[[paste0("gene_", p)]]
      L.gene <- inter[[paste0("gene_", q)]]
      rl[["R.gene"]] <- R.gene
      rl[["L.gene"]] <- L.gene
      if (paste0("uniprot_", p) %in% names(inter)) {
        rl[["R.uniprot"]] <- inter[[paste0("uniprot_", p)]]
      }
      if (paste0("uniprot_", q) %in% names(inter)) {
        rl[["L.uniprot"]] <- inter[[paste0("uniprot_", q)]]
      }
      if (paste0("name_", p) %in% names(inter)) {
        rl[["R.name"]] <- inter[[paste0("name_", p)]]
      }
      if (paste0("name_", q) %in% names(inter)) {
        rl[["L.name"]] <- inter[[paste0("name_", q)]]
      }
      return(rl)
    }
  )
  rl_df <- as.data.frame(do.call(rbind, rl_list))
  return(rl_df)
}
