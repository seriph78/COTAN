#' COTAN object
#'
#' The COTAN object for the ERCC dataset.
#'
#' @format A structure with:
#' \describe{
#'   \item{raw}{the raw dataset}
#'   \item{genesCoex}{a symmetric matrix expressing the gene/gene co-expression}
#'   \item{cellsCoex}{a symmetric matrix expressing the cell/cell co-expression}
#'   \item{nu}{UDE}
#'   \item{lambda}{average gene expression}
#'   \item{dispersion}{negative binomial dispersion}
#'   \item{metaDataset}{meta data for the dataset}
#'   \item{metaGenes}{meta data for the genes}
#'   \item{metaCells}{meta data for the cells}
#'   \item{clustersCoex}{list of co-expression matrices. One per each clusterization in metaCells columns}
#' }
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

