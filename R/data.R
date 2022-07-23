#' COTAN object
#'
#' The COTAN object for the ERCC dataset.
#'
#' @format A structure with:
#' \describe{
#'   \item{raw}{the raw dataset: 88 fake genes for 1015 fake cells}
#'   \item{raw.norm}{raw divided for nu}
#'   \item{coex}{}
#'   \item{nu}{UDE}
#'   \item{lambda}{ average gene expression}
#'   \item{a}{}
#'   \item{hk}{genes expressed in all cells}
#'   \item{n_cells}{final number of cells}
#'   \item{meta}{meta data}
#'}
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/ercc?}
"ERCC.cotan"

#' Raw sample dataset
#'
#' A subsample of a real sc-RNAseq dataset
#'
#' @format A data frame with 2000 genes and 815 cells:
#'
#' @source GEO GSM2861514
"raw.dataset"

#' Raw ERCC dataset
#'
#' ERCC dataset
#'
#' @format A data frame
#'
#' @source ERCC
"ERCCraw"

#' test.dataset.col
#'
#' Data set sampled from the raw one selecting 5000 random cells and 2000 genes:
#' half randomly and half randomly between the top differentially expressed genes
#' (GDI score > than 2). This is a dataset sampled from GEO GSE121380 Col1-6 
#' after Seurat cleaning.
#'
#'
#' @format A sparce matrix with cells as col names and genes as row names.
"test.dataset.col"

