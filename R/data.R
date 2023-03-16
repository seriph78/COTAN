#' Raw sample dataset
#'
#' A subsample of a real sc-RNAseq dataset
#'
#' @format A data frame with 2000 genes and 815 cells:
#'
#' @source GEO GSM2861514
#'
"raw.dataset"

#' Raw ERCC dataset
#'
#' ERCC dataset
#'
#' @format A data frame
#'
#' @source ERCC
#'
"ERCCraw"

#' test.dataset.col
#'
#' Data set sampled from the raw one selecting 5000 random cells and 2000 genes:
#' half randomly and half randomly between the top differentially expressed
#' genes (GDI score > than 2). This is a dataset sampled from GEO GSE121380
#' Col1-6 after Seurat cleaning.
#'
#' @format A sparse matrix with cells as col names and genes as row names.
#'
"test.dataset.col"
