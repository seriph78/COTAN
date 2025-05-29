#' @importFrom zeallot zeallous
.onLoad <- function(libname, pkgname) {
  zeallous()
}

#' @importFrom rlang .data
utils::globalVariables(c(".data"))
utils::globalVariables(c(".x"))
