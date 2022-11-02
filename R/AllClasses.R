#' Definition of COTAN class
#' @slot raw raw UMI count matrix 𝑛×𝑚 (gene number × cell number)
#' @slot rawNorm raw UMI count matrix divided by UDE, 𝑛×𝑚
#' @slot coex correlation of COTAN between genes, 𝑛×𝑛
#' @slot nu vector that stores the estimated UDE, size 𝑚
#' @slot lambda vector to store the average for the gene expression, size 𝑛
#' @slot dispetrion vector to store all
#' the negative binomial dispersion factors, size 𝑛.
#' @slot hKGenes house-keeping genes. It is a vector to store the name 
#' of the genes with positive UMI count in every single cell of the sample
#' @slot nCells number of the cells in the sample (𝑚)
#' @slot metaDataset data.frame
#' @slot metaCells data.frame
#' @slot clustersCoex coex
setClass("COTAN",
         slots = c(
           raw          = "dgCMatrix",
           rawNorm      = "dgCMatrix",
           coex         = "ANY",
           cellsCoex    = "ANY",
           nu           = "vector",
           lambda       = "vector",
           dispersion   = "vector",
           hkGenes      = "vector",
           nCells       = "numeric",
           metaDataset  = "data.frame",
           metaCells    = "data.frame",
           clustersCoex = "list"
         )
)

# constructor of the COTAN CLASS
#' @export 
COTAN <- function(raw = "ANY") {
  raw <- methods::as(as.matrix(raw), "sparseMatrix")
  new("COTAN", raw = raw)
}


#' scCOTAN-class (for legacy usage)
#'
#' Define my COTAN structure
#' @slot raw ANY. To store the raw data matrix
#' @slot raw.norm ANY. To store the raw data matrix divided for the cell
#' efficiency estimated (nu)
#' @slot coex ANY. The coex matrix (sparce)
#' @slot nu vector.
#' @slot lambda vector.
#' @slot a vector.
#' @slot hk vector.
#' @slot n_cells numeric.
#' @slot meta data.frame.
#' @slot yes_yes ANY.
#' @slot clusters vector.
#' @slot cluster_data data.frame.
#'
#' @return the object class
#' @export
#' @examples
#'
#' data("ERCCraw")
#' obj <- new("scCOTAN", raw = data)
#'
setClass("scCOTAN",
         slots = c(
           raw = "ANY",
           raw.norm = "ANY",
           coex = "ANY",
           nu = "vector",
           lambda = "vector",
           a = "vector",
           hk = "vector",
           n_cells = "numeric",
           meta = "data.frame",
           yes_yes = "ANY",
           clusters = "vector",
           cluster_data = "data.frame"
         )
) -> scCOTAN



## Automatically convert an object from class "scCOTAN" into "COTAN"
setIs("scCOTAN",
      "COTAN",
      #test = function(obj) {is.null(obj@yes_yes)},
      coerce = function(obj) {
        new("COTAN",
            raw          = obj@raw,
            rawNorm      = obj@raw.norm,
            coex         = obj@coex,
            cellsCoex    = NULL,
            nu           = obj@nu,
            lambda       = obj@lambda,
            dispersion   = obj@a,
            hkGenes      = obj@hk,
            nCells       = obj@n_cells,
            metaDataset  = obj@meta,
            metaCells    = as.data.frame(obj@clusters, names(obj@raw)),
            clustersCoex = as.list(obj@cluster_data) )
      },
      replace = function(obj, value) {
        obj@raw          <- value@raw
        obj@rawNorm      <- value@raw.norm
        obj@coex         <- value@coex
        obj@cellCoex     <- NULL
        obj@clustersCoex <- as.list(value@cluster_data)
        obj@nu           <- value@nu
        obj@lambda       <- value@lambda
        obj@dispersion   <- value@a
        obj@hkGenes      <- value@hk
        obj@n_cells      <- value@n_cells
        obj@metaDataset  <- value@meta
        obj@metaCells    <- as.data.frame(value@clusters, names(value@raw))
        obj@clustersCoex <- as.list(value@cluster_data)
        obj})