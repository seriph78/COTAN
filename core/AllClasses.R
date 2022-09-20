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
