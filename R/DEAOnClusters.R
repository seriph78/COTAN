
runSingleDEA <- function(clName, probZero, zeroOne, cellsInList) {
  cellsIn <- cellsInList[[clName]]
  if (!any(cellsIn)) {
    warning("Cluster '", clName, "' has no cells assigned to it!")
  }

  numCells <- ncol(zeroOne)
  numCellsIn  <- sum(cellsIn)
  numCellsOut <- numCells - numCellsIn

  observedYI <- rowSums(zeroOne[, cellsIn, drop = FALSE])

  expectedNI <- rowsums(probZero[,  cellsIn, drop = FALSE], parallel = TRUE)
  expectedNO <- rowsums(probZero[, !cellsIn, drop = FALSE], parallel = TRUE)
  expectedYI <- numCellsIn  - expectedNI
  expectedYO <- numCellsOut - expectedNO

  clCoex <- (observedYI  - expectedYI) / sqrt(numCells) *
    sqrt(1.0 / pmax(1.0, expectedNI) +
         1.0 / pmax(1.0, expectedNO) +
         1.0 / pmax(1.0, expectedYI) +
         1.0 / pmax(1.0, expectedYO))

  return(clCoex)
}

runDEA <- function(clNames, probZero, zeroOne, cellsInList, cores) {
  if (FALSE) {
    res <- parallel::mclapply(
      clNames,
      runSingleDEA,
      probZero = probZero,
      zeroOne = zeroOne,
      cellsInList = cellsInList,
      mc.cores = cores)

    # spawned errors are stored as try-error classes
    resError <- unlist(lapply(res, inherits, "try-error"))
    if (any(resError)) {
      stop(paste(res[which(resError)[[1L]]]), call. = FALSE)
    }
    return(res)
  } else {
    return(lapply(
      clNames,
      runSingleDEA,
      probZero = probZero,
      zeroOne = zeroOne,
      cellsInList = cellsInList))
  }
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
#' @param cores number of cores to use. Default is 1.
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
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @importFrom parallelly supportsMulticore
#' @importFrom parallelly availableCores
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom Rfast rowsums
#'
#' @rdname HandlingClusterizations
#'
DEAOnClusters <- function(objCOTAN, clName = "", clusters = NULL, cores = 1L) {
  logThis("Differential Expression Analysis - START", logLevel = 2L)

  # as sharing the enviroment with the forked processes takes a LOT of memory
  # it is much better to use only a few separate processes at a time:
  # dividing bt 4 makes possible to re-use the same number of cores
  # as the one passed in for the dispersion estimator
  cores <- handleMultiCore(min(4L, cores %/% 4L + 1L))

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

  cellsInList <- lapply(clustersList, function(cl) {getCells(objCOTAN) %in% cl})

  numSplits <- length(cellsInList)
  splitStep <- max(1L, cores)

  coexCls <- list()

  pBegin <- 1L
  while (pBegin <= numSplits) {
    logThis("*", appendLF = FALSE, logLevel = 1L)

    pEnd <- min(pBegin + splitStep - 1L, numSplits)

    logThis(paste0(" Executing ", (pEnd - pBegin + 1L), " DEA batches from",
                   " [", pBegin, ":", pEnd, "]"),
            logLevel = 3L)

    res <- NULL
    resError <- "No errors"
    failCount <- 0L
    while (!is_null(resError) && failCount < 3L) {
      failCount <- failCount + 1L
      c(res, resError) %<-%
        tryCatch(list(runDEA(clNames = names(cellsInList)[pBegin:pEnd],
                             probZero = probZero, zeroOne = zeroOne,
                             cellsInList = cellsInList, cores = cores), NULL),
                 error = function(e) {
                   logThis(paste("In DEA batches -", e), logLevel = 2L)
                   list(NULL, e) })
    }

    assert_that(is_null(resError),
                msg = paste("DEA batches failed", failCount,
                            "times with", resError))

    coexCls <- append(coexCls, res)

    pBegin <- pEnd + 1L
  }
  logThis("", logLevel = 1L)

  coexDF <- as.data.frame(coexCls)
  colnames(coexDF) <- names(cellsInList)

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

  normData <- getNormalizedData(objCOTAN)

  assert_that(!is_empty(getLambda(objCOTAN)),
              msg = "lambda must not be empty, estimate it")

  floorAverage <- getLambda(objCOTAN) * floorLambdaFraction

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
    averageIn  <-
      rowSums(normData[,  cellsIn, drop = FALSE]) / max(numCellsIn,  1L)
    averageOut <-
      rowSums(normData[, !cellsIn, drop = FALSE]) / max(numCellsOut, 1L)

    logAverageIn  <- log10(pmax(averageIn,  floorAverage))
    logAverageOut <- log10(pmax(averageOut, floorAverage))

    lfc <- logAverageIn - logAverageOut

    rm(averageIn, averageOut, logAverageIn, logAverageOut)
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
#' @param cores number of cores to use. Default is 1.
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
#' @importFrom rlang is_empty
#' @importFrom rlang is_null
#'
#' @rdname HandlingClusterizations
#'
distancesBetweenClusters <- function(objCOTAN, clName = "",
                                     clusters = NULL, coexDF = NULL,
                                     useDEA = TRUE, cores = 1L,
                                     distance = NULL) {
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
      coexDF <- DEAOnClusters(objCOTAN, clusters = clusters, cores = cores)
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
