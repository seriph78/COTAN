
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
#' @param method *p-value* multi-test adjustment method. Defaults to
#'   `"none"` (i.e. no adjustment)
#'
#' @return `pValueFromDEA()` returns a `data.frame` containing the *p-values*
#'   corresponding to the given `COEX` adjusted for *multi-test*
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom stats pchisq
#' @importFrom stats p.adjust
#'
#' @importFrom assertthat assert_that
#'
#' @rdname HandlingClusterizations
#'
pValueFromDEA <- function(coexDF, numCells, method = "none") {

  pValue <- pchisq(as.matrix(coexDF^2L * numCells),
                             df = 1L, lower.tail = FALSE)

  if (anyNA(pValue)) {
    warning("Got some NA in the p-value",
            toString(which(is.na(pValue), arr.ind = TRUE)))
  }

  if (method != "none") {
    for (cl in colnames(pValue)) {
      pValue[, cl] <- p.adjust(pValue[, cl], method = method, n = nrow(coexDF))
    }
  }

  if (anyNA(pValue)) {
    warning("Got some NA in the adjusted p-value",
            toString(which(is.na(pValue), arr.ind = TRUE)))
  }

  return(as.data.frame(pValue))
}


#
#
#' @details `logFoldChangeOnClusters()` is used to get the log difference of the
#'   expression levels for each *cluster* in the given *clusterization* against
#'   the rest of the data-set
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName`
#'
#' @return `logFoldChangeOnClusters()` returns the log-expression-change
#'   `data.frame` for the genes in each *cluster*
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
logFoldChangeOnClusters <- function(objCOTAN, clName = "", clusters = NULL) {
  logThis("Log Fold Change Analysis - START", logLevel = 2L)

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  clustersList <- toClustersList(clusters)

  numCells <- getNumCells(objCOTAN)

  normData <- getNormalizedData(objCOTAN)

  lfcDF <- data.frame()

  for (cl in names(clustersList)) {
    logThis("*", appendLF = FALSE, logLevel = 1L)
    logThis(paste0(" analysis of cluster: '", cl, "' - START"), logLevel = 3L)

    cellsIn <- getCells(objCOTAN) %in%  clustersList[[cl]]

    numCellsIn  <- sum(cellsIn)
    numCellsOut <- numCells - numCellsIn

    if (numCellsIn == 0L) {
      warning("Cluster '", cl, "' has no cells assigned to it!")
    }

    # log of average expression inside/outside the cluster
    logAverageIn  <-
      log10(pmax(rowSums(normData[,  cellsIn, drop = FALSE]), 1L) /
              max(numCellsIn,  1L))
    logAverageOut <-
      log10(pmax(rowSums(normData[, !cellsIn, drop = FALSE]), 1L) /
              max(numCellsOut, 1L))

    lfc <- logAverageIn - logAverageOut

    rm(logAverageIn, logAverageOut)
    gc()

    lfcDF <- setColumnInDF(lfcDF, colToSet = lfc,
                           colName = cl, rowNames = rownames(normData))

    logThis(paste0("* analysis of cluster: '", cl, "' - DONE"), logLevel = 3L)
  }
  logThis("", logLevel = 1L)

  logThis("Log Fold Change Analysis - DONE", logLevel = 2L)

  return(lfcDF)
}


#'
#'
#' @details `distancesBetweenClusters()` is used to obtain a model independent
#'   distance between the clusters. It is the `"euclidian"` distance between the
#'   per clusters averages of the corresponding columns of the
#'   *Zero/One projection* for the raw data matrix
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName`
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#' @param useDEA Boolean indicating whether to use the *DEA* to define the
#'   distance; alternatively it will use the average *Zero-One* counts, that is
#'   faster but less precise.
#' @param distance type of distance to use. Default is `"cosine"` for *DEA* and
#'   `"euclidean"` for *Zero-One*. Can be chosen among those supported by
#'   [parallelDist::parDist()]
#'
#' @return `distancesBetweenClusters()` returns a `dist` object
#'
#' @export
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom parallelDist parDist
#'
#' @importFrom assertthat assert_that
#'
#' @rdname HandlingClusterizations
#'
distancesBetweenClusters <- function(objCOTAN, clName = "",
                                     clusters = NULL, coexDF = NULL,
                                     useDEA = TRUE, distance = NULL) {
  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  clList <- toClustersList(clusters)

  if (isTRUE(useDEA)) {
    if (is.null(distance)) {
      distance <- "cosine"
    }

    if (is_empty(coexDF) && (clName %in% getClusterizations(objCOTAN))) {
        coexDF <- getClusterizationData(objCOTAN, clName = clName)[["coex"]]
    }
    if (is_empty(coexDF)) {
      coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)
    }

    # merge small cluster based on distances
    return(parDist(t(as.matrix(coexDF)), method = distance,
                   diag = TRUE, upper = TRUE))
  } else {
    if (is.null(distance)) {
      distance <- "euclidean"
    }

    zeroOne <- getZeroOneProj(objCOTAN)

    zeroOneClAvg <- data.frame(row.names = getGenes(objCOTAN))
    for (cluster in clList) {
      zeroOneClAvg <- cbind(zeroOneClAvg, rowMeans(zeroOne[, cluster]))
    }

    # ensure no zeros in the matrix
    zeroOneClAvg[zeroOneClAvg == 0.0] <- 1.0e-6
    colnames(zeroOneClAvg) <- names(clList)

    return(parDist(t(zeroOneClAvg), method = distance,
                   diag = TRUE, upper = TRUE))
  }
}
