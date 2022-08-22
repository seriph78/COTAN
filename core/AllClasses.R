#' scCOTAN-class
#'
#' Define my COTAN structure
#' @slot raw ANY. To store the raw data matrix
#' @slot raw.norm ANY. To store the raw data matrix divided for the cell
#' efficiency estimated (nu)
#' @slot coex ANY. The coex matrix (sparce)
#' @slot nu vector.
#' @slot lambda vector.
#' @slot dispertion vector.
#' @slot hk vector.
#' @slot n_cells numeric.
#' @slot meta data.frame..
#' @slot clusters vector. 
#' @slot cluster_data data.frame.
#'
#' @return the object class
#' @examples
#'
#' data("ERCCraw")
#' obj <- new("scCOTAN", raw = data)
#'
setClass("COTAN",
  slots = c(
    raw        = "dgCMatrix",
    rawNorm    = "dgCMatrix",
    nu         = "vector",
    lambda     = "vector",
    hKGenes    = "vector",
    dispertion = "vector",
    nCells     = "numeric",
    meta.dataset = "data.fame",
    meta.cells = "data.frame",
    clusters.coex = "list"
  )
)

# constructor of the COTAN CLASSS
#' @export 
COTAN <- function(raw = "ANY") {
  raw <- methods::as(as.matrix(raw), "sparseMatrix")
  new("COTAN", raw = raw)
}
