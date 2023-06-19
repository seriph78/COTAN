
#'
#'
#' @details `DEAOnClusters()` is used to run the Differential Expression
#'   analysis using the `COTAN` contingency tables on each *cluster* in the
#'   given *clusterization*
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters a `vector` or `factor` with the cell *clusterization* to be
#'   analyzed. If empty, the function will use the last *clusterization* stored
#'   in the `COTAN` object
#'
#' @return `DEAOnClusters()` returns the coexpression `data.frame` for the genes
#'   in each *cluster*
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
DEAOnClusters <- function(objCOTAN, clusters = NULL) {
  if (is_empty(clusters)) {
    # pick the last clusterization
    clusters <- getClusterizationData(objCOTAN)[["clusters"]]
  } else if (!inherits(clusters, "factor")) {
    clusters <- factor(clusters)
  }

  assert_that(length(clusters) == getNumCells(objCOTAN),
              msg = "Passed/retrieved clusterization has the wrong size")

  assert_that(all(names(clusters) %in% getCells(objCOTAN)),
              msg = paste("Some cells in the clusterization",
                          "are not part of the 'COTAN' object"))

  logThis("Differential Expression Analysis - START", logLevel = 2L)

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

    observedYI <- rowSums(zeroOne[, cellsIn])

    expectedNI <- rowSums(probZero[,  cellsIn])
    expectedNO <- rowSums(probZero[, !cellsIn])
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
#' @param coexDF the co-expression `data.frame` for the genes in each *cluster*
#' @param numCells the number of cells in all *clusters*
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
