
#'
#'
#' @details `DEAOnClusters()` is used to run the Differential Expression
#'   analysis using the `COTAN` contingency tables on each *cluster* in the
#'   given *clusterization*
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters a `vector` or `factor` with the cell *clusterization* to
#'   be analyzed. If empty, the function will use the last *clusterization*
#'   stored in the `COTAN` object
#'
#' @return `DEAOnClusters()` returns a `list` with two objects:
#'   * "coex"    - the coexpression `data.frame` for the genes in each *cluster*
#'   * "p-value" - the corresponding p-values `data.frame`
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

  clustersCoex <- data.frame()
  clustersPVal <- data.frame()

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

    pValue <- pchisq(coex^2L * numCells, df = 1L, lower.tail = FALSE)

    if (anyNA(pValue)) {
      warning("Got some NA in the p-value",
              toString(which(is.na(pValue), arr.ind = TRUE)))
    }

    rm(expectedYO, expectedYI, expectedNO, expectedNI)
    rm(observedYI)
    gc()

    clustersCoex <- setColumnInDF(clustersCoex, colToSet = coex,
                                  colName = cl, rowNames = rownames(zeroOne))
    clustersPVal <- setColumnInDF(clustersPVal, colToSet = pValue,
                                  colName = cl, rowNames = rownames(zeroOne))

    logThis(paste0("* analysis of cluster: '", cl, "' - DONE"), logLevel = 3L)
  }
  logThis("", logLevel = 1L)

  logThis("Differential Expression Analysis - DONE", logLevel = 2L)

  return(list("coex" = clustersCoex, "p-value" = clustersPVal))
}
