
#' clustersDeltaExpression
#'
#' @description This function estimates the change in genes' expression inside
#'   the cluster compared to the average situation in the data set.
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters The clusterization. If none given the last available
#'   clusterization will be used, as it is probably the most significant!
#'
#' @returns a `data.frame` with the weighted discrepancy of the expression of
#'   each gene within the cluster against model expectations
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = FALSE,
#'                                          outDir = tempdir())
#'
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#'
#' deltaExpression <- clustersDeltaExpression(objCOTAN, clusters)
#'
#' @rdname clustersDeltaExpression
#'
clustersDeltaExpression <- function(objCOTAN, clusters = NULL) {
  logThis("clustersDeltaExpression - START", logLevel = 2L)

  if (is_empty(clusters)) {
    clusters <- getClusterizationData(objCOTAN)[["clusters"]]
  }
  assert_that(!is_empty(clusters),
              msg = "No clusterization given or present in the COTAN object")
  assert_that(!is_empty(names(clusters)),
              msg = "No names attached to the given clusterization")
  assert_that(setequal(names(clusters), getCells(objCOTAN)),
              msg = "Non compatible clusterization")

  clustersList <- toClustersList(clusters)

  zeroOne <- getZeroOneProj(objCOTAN)
  probOne <- 1.0 - funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

  deltaExpression <- data.frame()
  for (cl in names(clustersList)) {
    logThis(paste0("Handling cluster '", cl, "'"),
            appendLF = FALSE, logLevel = 3L)

    cellsPos <- getCells(objCOTAN) %in% clustersList[[cl]]

    observedYI <- rowSums(zeroOne[, cellsPos])
    expectedYI <- rowSums(probOne[, cellsPos])
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
