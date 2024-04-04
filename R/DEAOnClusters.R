
#'
#'
#' @details `DEAOnClusters()` is used to run the Differential Expression
#'   analysis using the `COTAN` contingency tables on each *cluster* in the
#'   given *clusterization*
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName`
#'
#' @return `DEAOnClusters()` returns the co-expression `data.frame` for the
#'   genes in each *cluster*
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom assertthat assert_that
#'
#' @rdname HandlingClusterizations
#'
DEAOnClusters <- function(objCOTAN, clName = "", clusters = NULL) {
  logThis("Differential Expression Analysis - START", logLevel = 2L)

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  assert_that(estimatorsAreReady(objCOTAN),
              msg = paste("Estimators lambda, nu, dispersion are not ready:",
                          "Use proceeedToCoex() to prepare them"))

  clustersList <- toClustersList(clusters)

  zeroOne <- getZeroOneProj(objCOTAN)

  probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

  numCells <- getNumCells(objCOTAN)

  coexDF <- data.frame()

  for (cl in names(clustersList)) {
    gc()
    logThis("*", appendLF = FALSE, logLevel = 1L)
    logThis(paste0(" analysis of cluster: '", cl, "' - START"), logLevel = 3L)

    cellsIn <- getCells(objCOTAN) %in%  clustersList[[cl]]

    numCellsIn  <- sum(cellsIn)
    numCellsOut <- numCells - numCellsIn

    if (numCellsIn == 0L) {
      warning("Cluster '", cl, "' has no cells assigned to it!")
    }

    observedYI <- rowSums(zeroOne[, cellsIn, drop = FALSE])

    expectedNI <- rowSums(probZero[,  cellsIn, drop = FALSE])
    expectedNO <- rowSums(probZero[, !cellsIn, drop = FALSE])
    expectedYI <- numCellsIn  - expectedNI
    expectedYO <- numCellsOut - expectedNO

    coex <- (observedYI  - expectedYI) / sqrt(numCells) *
              sqrt(1.0 / pmax(1.0, expectedNI) +
                   1.0 / pmax(1.0, expectedNO) +
                   1.0 / pmax(1.0, expectedYI) +
                   1.0 / pmax(1.0, expectedYO))

    rm(expectedYO, expectedYI, expectedNO, expectedNI)
    rm(observedYI)
    gc()

    coexDF <- setColumnInDF(coexDF, colToSet = coex,
                            colName = cl, rowNames = rownames(zeroOne))

    logThis(paste0("* analysis of cluster: '", cl, "' - DONE"), logLevel = 3L)
  }
  logThis("", logLevel = 1L)

  logThis("Differential Expression Analysis - DONE", logLevel = 2L)

  return(coexDF)
}


#'
#'
#' @details `pValueFromDEA()` is used to convert to *p-value* the Differential
#'   Expression analysis using the `COTAN` contingency tables on each *cluster*
#'   in the given *clusterization*
#'
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#' @param numCells the number of overall cells in all *clusters*
#'
#' @return `pValueFromDEA()` returns a `data.frame` with the *p-values*
#'   corresponding to the given *coex*
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom stats pchisq
#'
#' @importFrom assertthat assert_that
#'
#' @rdname HandlingClusterizations
#'
pValueFromDEA <- function(coexDF, numCells) {

  pValue <- pchisq(as.matrix(coexDF^2L * numCells),
                             df = 1L, lower.tail = FALSE)

  if (anyNA(pValue)) {
    warning("Got some NA in the p-value",
            toString(which(is.na(pValue), arr.ind = TRUE)))
  }

  return(as.data.frame(pValue))
}
