.onLoad <- function(libname, pkgname) {
  if (requireNamespace("zeallot") &&
      utils::packageVersion("zeallot") >= "2.0.0") {
    zeallot::zeallous()
  }
}

#' @importFrom rlang .data
utils::globalVariables(c(".data", ".x"))
