# internal scripts for the create_domino() function

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

