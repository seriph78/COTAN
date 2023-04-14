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
#' @format a `data.frame`
#'
#' @source ERCC
#'
"ERCCraw"

#' test.dataset
#'
#' @description `test.dataset` is an artificial data set obtained by sampling
#'   target negative binomial distributions on a set of `600` genes on `2` two
#'   cells *clusters* of `600` cells each. Each *clusters* has its own set of
#'   parameters for the distributions even, but a fraction of the genes has the
#'   same expression in both *clusters*.
#'
#' @format a `data.frame`
#'
#' @rdname test.dataset
#'
"test.dataset"

#' @details `test.dataset.clusters1` is the clusterization obtained running
#'   `cellsUniformClustering()` on the `test.dataset`
#'
#' @format a `character array`
#'
#' @rdname test.dataset
#'
"test.dataset.clusters1"

#' @details `test.dataset.clusters2` is the clusterization obtained running
#'   `mergeUniformCellsClusters()` on the `test.dataset` using the
#'   `test.dataset.clusters1`
#'
#' @format a `character array`
#'
#' @rdname test.dataset
#'
"test.dataset.clusters2"
