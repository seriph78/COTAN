
#' DEAOnClusters
#'
#' @description This function is used to run the differential expression
#'   analysis using the `COTAN` contingency tables on each cluster in the given
#'   clusterization
#'
#' @details The formulae for the coex are similar to those used in the
#'   [calculateCoex()] method, with the **role** of the second gene taken by the
#'   In/Out status of the cells with respect to each cluster.
#' @seealso [calculatePValue()] for details about the associated p-values
#'
#' @param objCOTAN A COTAN object
#' @param clusterization a `vector` or `factor` with the cell clusterization to
#'   be analysed. If empty, the function will use the last clusterization stored
#'   in the `objCOTAN`
#'
#' @return a `list` with two objects:
#'   * 'coex'    - the coexpression `data.frame` for the genes in each cluster,
#'   * 'p-value' - the corresponding p-values `data.frame`.
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
#' list[coexDF, pvalDF] <- DEAOnClusters(objCOTAN)
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "clusters",
#'                               clusters = clusters, coexDF = coexDF)
#'
#' @rdname DEAOnClusters
#'
DEAOnClusters <- function(objCOTAN, clusterization = NULL) {
  if (is_empty(clusterization)) {
    clusterization <- getClusterizationData(objCOTAN)[["clusters"]]
  }

  stopifnot("Passed/retrieved clusterization has the wrong size" <-
              (length(clusterization) == getNumCells(objCOTAN)))

  stopifnot("Some cells in the clusterization are not part of the 'COTAN' object" <-
              all(names(clusterization) %in% getCells(objCOTAN)))

  logThis("Differential Expression Analysis - START", logLevel = 2)

  clustersList <- toClustersList(clusterization)

  noHKFlags <- flagNotHousekeepingGenes(objCOTAN)

  zeroOne <- getZeroOneProj(objCOTAN)[noHKFlags, ]

  muEstimator <- calculateMu(objCOTAN)[noHKFlags, ]

  probZero <- funProbZero(getDispersion(objCOTAN)[noHKFlags], muEstimator)

  numCells <- getNumCells(objCOTAN)

  exp_yes <- rowSums(zeroOne)
  exp_no  <- numCells - exp_yes

  clusters_coex <- data.frame()
  clusters_pval <- data.frame()

  for (cl in names(clustersList)) {
    gc()
    logThis("*", appendLF = FALSE, logLevel = 1)
    logThis(paste0(" analysis of cluster: '", cl, "' - START"), logLevel = 3)

    cellsIn <- getCells(objCOTAN) %in%  clustersList[[cl]]

    numCellsIn  <- sum(cellsIn)
    numCellsOut <- numCells - numCellsIn

    observedYI  <- rowSums(zeroOne[, cellsIn])

    estimatedNI <- rowSums(probZero[, cellsIn])
    estimatedNO <- rowSums(probZero[,!cellsIn])
    estimatedYI <- numCellsIn  - estimatedNI
    estimatedYO <- numCellsOut - estimatedNO

    coex <- (observedYI  - estimatedYI) / sqrt(numCells) *
              sqrt(1 / pmax(1, estimatedNI) +
                   1 / pmax(1, estimatedNO) +
                   1 / pmax(1, estimatedYI) +
                   1 / pmax(1, estimatedYO))

    pValue <- pchisq(coex^2 * numCells, df <- 1, lower.tail <- FALSE)

    if (any(is.na(pValue))) {
      warning(paste("Got some NA in the p-value",
                    which(is.na(S), arr.ind = T), collapse = ", "))
    }

    rm(estimatedYO, estimatedYI, estimatedNO, estimatedNI)
    rm(observedYI)
    gc()

    clusters_coex <- setColumnInDF(clusters_coex, colToSet = coex,
                                   colName = cl, rowNames = rownames(zeroOne))
    clusters_pval <- setColumnInDF(clusters_pval, colToSet = pValue,
                                   colName = cl, rowNames = rownames(zeroOne))

    logThis(paste0("* analysis of cluster: '", cl, "' - DONE"), logLevel = 3)
  }
  logThis("", logLevel = 1)

  logThis("Differential Expression Analysis - DONE", logLevel = 2)

  return( list("coex" = clusters_coex, "p-value" = clusters_pval) )
}

