setClass("COTAN",
  slots = c(
    raw        = "dgCMatrix",
    rawNorm    = "dgCMatrix",
    nu         = "vector",
    lambda     = "vector",
    hKGenes    = "vector",
    dispertion = "vector",
    nCells     = "numeric"
  )
)

# constructor of the COTAN CLASSS
#' @export 
COTAN <- function(raw = "ANY") {
  raw <- methods::as(as.matrix(raw), "sparseMatrix")
  new("COTAN", raw = raw)
}
