
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
#' @examples
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
#'                                          GEO = "S"
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = TRUE,
#'                                          outDir = tempdir())
#'
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = TRUE,
#'                                    outDir = tempdir())
#'
#' deltaExpression <- clustersDeltaExpression(objCOTAN, clusters)
#'
#' @rdname clustersDeltaExpression
#'
clustersDeltaExpression <- function(objCOTAN, clusters = NULL) {
  logThis("clustersDeltaExpression - START", logLevel = 2)

  if (is_empty(clusters)) {
    clusters <- getClusterizationData(objCOTAN)[["clusters"]]
  }
  stopifnot("No clusterization given or present in the COTAN object" <- !is_empty(clusters))
  stopifnot("No names attached to the given clusterization" <- !is_empty(names(clusters)))
  stopifnot("Non compatible clusterization" <- setequal(names(clusters), getCells(objCOTAN)))

  clustersList <- toClustersList(clusters)

  zeroOne <- getZeroOneProj(objCOTAN)
  probOne <- 1 - funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

  deltaExpression <- data.frame()
  for (cl in names(clustersList)) {
    logThis(paste0("Handling cluster '", cl, "'"), appendLF = FALSE, logLevel = 3)

    cellsPos <- getCells(objCOTAN) %in% clustersList[[cl]]

    observedYI <- rowSums(zeroOne[, cellsPos])
    expectedYI <- rowSums(probOne[, cellsPos])
    sumUDECluster <- sum(getNu(objCOTAN)[cellsPos])

    logThis(paste0("with mean UDE ", sumUDECluster / sum(cellsPos)), logLevel = 3)

    deltaE <- (observedYI - expectedYI) / sumUDECluster

    deltaExpression <- setColumnInDF(deltaExpression, colToSet = deltaE,
                                     colName = cl, rowNames = getGenes(objCOTAN))
  }

  logThis("clustersDeltaExpression - DONE", logLevel = 2)

  return(deltaExpression)
}
