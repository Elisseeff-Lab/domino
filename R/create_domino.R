# internal scripts for the create_domino() function

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

