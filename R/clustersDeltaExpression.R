
#'
#' @details `clustersDeltaExpression()` is a legacy function now superseded by
#'   [DEAOnClusters()]. It estimates the change in genes' expression inside the
#'   *cluster* compared to the average situation in the data set.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName`
#'
#' @returns `clustersDeltaExpression()` returns a `data.frame` with the
#'   \eqn{\nu} weighted discrepancy of the expression of each gene within the
#'   *cluster* against the corresponding model expectations
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom Rfast rowsums
#'
#' @importFrom assertthat assert_that
#'
#' @rdname COTAN_Legacy
#'
clustersDeltaExpression <- function(objCOTAN, clName = "", clusters = NULL) {
  logThis("clustersDeltaExpression - START", logLevel = 2L)

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  clustersList <- toClustersList(clusters)

  zeroOne <- getZeroOneProj(objCOTAN)
  probOne <- 1.0 - getProbabilityOfZero(objCOTAN)

  deltaExpression <- data.frame()
  for (cl in names(clustersList)) {
    logThis(paste0("Handling cluster '", cl, "'"),
            appendLF = FALSE, logLevel = 3L)

    cellsPos <- getCells(objCOTAN) %in% clustersList[[cl]]

    observedYI <- rowSums(zeroOne[, cellsPos])
    expectedYI <- rowsums(probOne[, cellsPos], parallel = TRUE)
    sumUDECluster <- sum(getNu(objCOTAN)[cellsPos])

    logThis(paste0(" with mean UDE ", sumUDECluster / sum(cellsPos)),
            logLevel = 3L)

    deltaE <- (observedYI - expectedYI) / sumUDECluster

    deltaExpression <- setColumnInDF(deltaExpression,
                                     colToSet = deltaE,
                                     colName = cl,
                                     rowNames = getGenes(objCOTAN))
  }

  logThis("clustersDeltaExpression - DONE", logLevel = 2L)

  return(deltaExpression)
}
