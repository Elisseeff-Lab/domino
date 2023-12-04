#' Check input arguments
#'
#' Accepts an object and rules for checking, stops if rules not met.
#'
#' @param arg The argument to check
#' @param allow_class Vector of allowed classes
#' @param allow_len Vector of allowed lengths
check_arg <- function(arg, allow_class = c("character"), allow_len = NULL) {
  argname <- deparse(substitute(arg))
  classes <- paste(allow_class, collapse = ",")
  lengths <- paste(allow_len, collapse = ",")

  if (!(class(arg) %in% allow_class)) {
    stop(sprintf("Class of %s must be one of: %s", argname, classes))
  }

  if (!is.null(allow_len)) {
    if (!(length(arg) %in% allow_len)) {
      stop(sprintf("Length of %s must be one of: %s", argname, lengths))
    }
  }
}

#' Read in data if an object looks like path to it.
#'
#' @param obj
#' @return obj the object itself in case its not a character
read_if_char <- function(obj) {
  if (is(obj, "character")) {
    obj <- read.csv(obj, stringsAsFactors = FALSE)
  }
  return(obj)
}


#' Change cases of True/False syntax from Python to TRUE/FALSE R syntax
#'
#' @param obj Object that will be converted
#' @return obj The converted object
conv_py_bools <- function(obj) {
  for (x in colnames(obj)) {
    bools <- sort(unique(obj[[x]]))
    if (identical(bools, c("False", "True"))) {
      obj[[x]] <- ifelse(obj[[x]] == "True", TRUE, FALSE)
    }
  }
  return(obj)
}
