#' @title Data-sets
#'
#' @description Simple data-sets included in the package
#'
#' @name Datasets
#'
NULL

# Raw sample dataset

#' @details `raw.dataset` is a sub-sample of a real *sc-RNAseq* data-set
#'
#' @format `raw.dataset` is a data frame with \eqn{2000} genes and \eqn{815}
#'   cells
#'
#' @usage data(raw.dataset)
#'
#' @source GEO GSM2861514
#'
#' @rdname Datasets
#'
"raw.dataset"

# Raw ERCC dataset

#' @details `ERCCRaw` dataset
#'
#' @format `ERCCRaw` is a `data.frame`
#'
#' @usage data(ERCCraw)
#'
#' @source ERCC
#'
#' @rdname Datasets
#'
"ERCCraw"

#' @details `test.dataset` is an artificial data set obtained by sampling target
#'   negative binomial distributions on a set of \eqn{600} genes on \eqn{2} two
#'   cells *clusters* of \eqn{600} cells each. Each *clusters* has its own set
#'   of parameters for the distributions even, but a fraction of the genes has
#'   the same expression in both *clusters*.
#'
#' @format `test.dataset` is a `data.frame` with \eqn{600} genes and \eqn{1200}
#'   cells
#'
#' @usage data(test.dataset)
#'
#' @rdname Datasets
#'
"test.dataset"

#' @details `test.dataset.clusters1` is the *clusterization* obtained running
#'   `cellsUniformClustering()` on the `test.dataset`
#'
#' @format `test.dataset.clusters1` is a `character array`
#'
#' @usage data(test.dataset.clusters1)
#'
#' @rdname Datasets
#'
"test.dataset.clusters1"

#' @details `test.dataset.clusters2` is the *clusterization* obtained running
#'   `mergeUniformCellsClusters()` on the `test.dataset` using the previous
#'   *clusterization*
#'
#' @format `test.dataset.clusters2` is a `character array`
#'
#' @usage data(test.dataset.clusters2)
#'
#' @rdname Datasets
#'
"test.dataset.clusters2"
