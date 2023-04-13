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

#' test.dataset
#'
#' Artificial data set obtained by sampling target negative binomial
#' distributions on a set of 600 genes on 2 two cells clusters of 600 cells
#' each. Each clusters has its own set of parameters for the distributions even,
#' but a fraction of the genes has the same expression in both clusters.
#'
#' @format A data frame
#'
"test.dataset"
"test.dataset.clusters1"
"test.dataset.clusters2"

#' sampled.dataset
#'
#' Data set sampled from the raw one selecting 5000 random cells and 2000 genes:
#' half randomly and half randomly between the top differentially expressed
#' genes (GDI score > than 2). This is a dataset sampled from GEO GSE121380
#' Col1-6 after Seurat cleaning.
#'
#' @format A data frame.
#'
"sampled.dataset"
