
runSingleDEA <- function(clName, cellsInList,
                         probZero, rowSumsProbZero, zeroOne) {
  logThis("*", appendLF = FALSE, logLevel = 1L)

  cellsIn <- cellsInList[[clName]]
  if (!any(cellsIn)) {
    warning("Cluster '", clName, "' has no cells assigned to it!")
  }
  logThis(paste0(" DEA on cluster '", clName,
                 "' with ", sum(cellsIn), " cells"), logLevel = 3L)

  numCells <- ncol(zeroOne)
  numCellsIn  <- sum(cellsIn)

  observedYI <- rowSums(zeroOne [ , cellsIn, drop = FALSE])

  expectedNI <- rowsums(probZero[ , cellsIn, drop = FALSE], parallel = TRUE)
  expectedNO <- rowSumsProbZero - expectedNI
  expectedYI <- numCellsIn  - expectedNI
  expectedYO <- (numCells - numCellsIn) - expectedNO

  clCoex <- (observedYI  - expectedYI) / sqrt(numCells) *
    sqrt(1.0 / pmax(1.0, expectedNI) +
         1.0 / pmax(1.0, expectedNO) +
         1.0 / pmax(1.0, expectedYI) +
         1.0 / pmax(1.0, expectedYO))

  return(clCoex)
}


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
#' @importFrom rlang is_null
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom Rfast rowsums
#' @importFrom Matrix rowSums
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
              msg = paste("Estimators `lambda`, `nu`, `dispersion` are not",
                          "ready: Use proceedToCoex() to prepare them"))

  clustersList <- toClustersList(clusters)

  zeroOne <- getZeroOneProj(objCOTAN)

  probZero <- getProbabilityOfZero(objCOTAN)

  rowSumsProbZero <- rowsums(probZero, parallel = TRUE)

  cellsInList <- lapply(clustersList, function(cl) getCells(objCOTAN) %in% cl)

  coexCls <- lapply(names(cellsInList), runSingleDEA,
                    cellsInList = cellsInList,
                    probZero = probZero, rowSumsProbZero = rowSumsProbZero,
                    zeroOne = zeroOne)

  logThis("", logLevel = 1L)

  coexDF <- as.data.frame(coexCls)
  colnames(coexDF) <- names(cellsInList)

  logThis("Differential Expression Analysis - DONE", logLevel = 2L)

  return(coexDF)
}

#' @details `clusterGeneContingencyTables()` returns the observed and expected
#'   contingency tables for a given gene and a given set of cells (a cluster).
#'   The implementation runs the same algorithms used to calculate the full
#'   observed/expected contingency tables used for DEA, but restricted to only
#'   the relevant gene and cluster, thus much faster and less memory intensive
#'
#' @param objCOTAN a `COTAN` object
#' @param gene a gene
#' @param cells a sub-set of the cells
#'
#' @return `clusterGeneContingencyTables()` returns a list containing the
#'   observed and expected contingency tables
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
clusterGeneContingencyTables <- function(objCOTAN, gene, cells) {
  allGenes <- getGenes(objCOTAN)
  allCells <- getCells(objCOTAN)

  assert_that(c(gene) %in% allGenes, msg = "the given gene is unknown")
  assert_that(all(cells %in% allCells),
              msg = "some of the given cells are unknown")

  dimnames <- list(c(paste(gene, "yes", sep = "."),
                     paste(gene, "no", sep = ".")),
                   c("cells.in", "cells.out"))

  # reorder and standardize cells list
  cellsIn <- allCells %in% cells

  # observed

  zeroOne <- sign(getRawData(objCOTAN)[gene, ])

  numCells    <- length(zeroOne)
  numCellsIn  <- sum(cellsIn)
  numCellsOut <- numCells - numCellsIn

  observedY  <- sum(zeroOne)

  observedYI <- sum(zeroOne[cellsIn])
  observedYO <- observedY   - observedYI
  observedNI <- numCellsIn  - observedYI
  observedNO <- numCellsOut - observedYO

  observedCT <- matrix(c(observedYI, observedNI, observedYO, observedNO),
                       ncol = 2L, nrow = 2L, dimnames = dimnames)

  # estimated
  lambda <- suppressWarnings(getLambda(objCOTAN))
  assert_that(!is_empty(lambda),
              msg = "`lambda` must not be empty, estimate it")

  nu <- suppressWarnings(getNu(objCOTAN))
  assert_that(!is_empty(nu),
              msg = "`nu` must not be empty, estimate it")

  dispersion <- suppressWarnings(getDispersion(objCOTAN))
  assert_that(!is_empty(dispersion),
              msg = "`dispersion` must not be empty, estimate it")

  probZero <- funProbZero(dispersion[gene], lambda[gene] * nu)

  if (anyNA(probZero)) {
    warning("Some NA in estimated probability of zero matrix")
  }

  expectedN  <- sum(probZero)

  expectedNI <- sum(probZero[cellsIn])
  expectedNO <- expectedN   - expectedNI
  expectedYI <- numCellsIn  - expectedNI
  expectedYO <- numCellsOut - expectedNO

  expectedCT <- matrix(c(expectedYI, expectedNI, expectedYO, expectedNO),
                       ncol = 2L, nrow = 2L, dimnames = dimnames)

  return(list("observed" = observedCT, "expected" = expectedCT))
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
#' @param adjustmentMethod *p-value* multi-test adjustment method, see
#'   [stats::p.adjust.methods()]. Defaults to `"none"` (i.e. no adjustment)
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
pValueFromDEA <- function(coexDF, numCells, adjustmentMethod) {
  pValue <- pchisq(as.matrix(coexDF^2L * numCells),
                             df = 1L, lower.tail = FALSE)

  if (anyNA(pValue)) {
    warning("Got some NA in the p-value",
            toString(which(is.na(pValue), arr.ind = TRUE)))
  }

  if (adjustmentMethod != "none") {
    for (cl in colnames(pValue)) {
      pValue[, cl] <- p.adjust(pValue[, cl],
                               method = adjustmentMethod, n = nrow(coexDF))
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
#' @param floorLambdaFraction Indicates the lower bound to the average count
#'   sums inside or outside the cluster for each gene as fraction of the
#'   relevant `lambda` parameter. Default is \eqn{5\%}
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
logFoldChangeOnClusters <- function(objCOTAN, clName = "", clusters = NULL,
                                    floorLambdaFraction = 0.05) {
  logThis("Log Fold Change Analysis - START", logLevel = 2L)

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  clustersList <- toClustersList(clusters)

  numCells <- getNumCells(objCOTAN)

  normData <- getNuNormData(objCOTAN)

  lambda <- suppressWarnings(getLambda(objCOTAN))
  assert_that(!is_empty(lambda),
              msg = "`lambda` must not be empty, estimate it")

  floorAverage <- lambda * floorLambdaFraction

  allSums <- rowSums(normData)

  lfcDF <- data.frame()

  for (cl in names(clustersList)) {
    logThis("*", appendLF = FALSE, logLevel = 1L)
    logThis(paste0(" Analysis of cluster: '", cl, "'"), logLevel = 3L)

    cellsIn <- getCells(objCOTAN) %in%  clustersList[[cl]]

    numCellsIn  <- sum(cellsIn)
    numCellsOut <- numCells - numCellsIn

    if (numCellsIn == 0L) {
      warning("Cluster '", cl, "' has no cells assigned to it!")
    }

    # log of average expression inside/outside the cluster
    inSums <- rowSums(normData[,  cellsIn, drop = FALSE])

    averageIn  <- inSums             / max(numCellsIn,  1L)
    averageOut <- (allSums - inSums) / max(numCellsOut, 1L)

    logAverageIn  <- log10(pmax(averageIn,  floorAverage))
    logAverageOut <- log10(pmax(averageOut, floorAverage))

    lfc <- logAverageIn - logAverageOut

    rm(averageIn, averageOut, logAverageIn, logAverageOut)
    gc()

    lfcDF <- setColumnInDF(lfcDF, colToSet = lfc,
                           colName = cl, rowNames = rownames(normData))
  }
  logThis("", logLevel = 1L)

  logThis("Log Fold Change Analysis - DONE", logLevel = 2L)

  return(lfcDF)
}


#'
#'
#' @details `distancesBetweenClusters()` is used to obtain a distance between
#'   the clusters. Depending on the value of the `useDEA` flag will base the
#'   distance on the *DEA* columns or the averages of the *Zero-One* matrix.
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
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_empty
#' @importFrom rlang is_null
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

  if (isTRUE(useDEA) && is_empty(coexDF) && !estimatorsAreReady(objCOTAN)) {
    logThis("cannot calculate DEA - falling back to case 'useDEA = FALSE'",
            logLevel = 1L)
    useDEA <- FALSE
  }

  if (isTRUE(useDEA)) {
    if (is_null(distance)) {
      distance <- "cosine"
    }

    if (is_empty(coexDF) && (clName %in% getClusterizations(objCOTAN))) {
        coexDF <- getClusterizationData(objCOTAN, clName = clName)[["coex"]]
    }

    if (is_empty(coexDF)) {
      coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)
    }

    # merge small cluster based on distances
    return(calcDist(t(as.matrix(coexDF)), method = distance,
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

    return(calcDist(t(zeroOneClAvg), method = distance,
                    diag = TRUE, upper = TRUE))
  }
}
