
#' @details `calculateGenesCE()` is used to calculate the discrepancy between
#'   the expected probability of zero and the observed zeros across all cells
#'   for each gene as *cross-entropy*: \eqn{-\sum_{c}{\mathbb{1}_{X_c == 0}
#'    \log(p_c) - \mathbb{1}_{X_c != 0} \log(1 - p_c)}} where \eqn{X_c} is the
#'   observed count and \eqn{p_c} the probability of zero
#'
#' @param objCOTAN a `COTAN` object
#'
#' @return `calculateGenesCE()` returns a named `array` with the *cross-entropy*
#'   of each gene
#'
#' @importFrom Rfast rowsums
#'
#' @export
#'
#' @rdname GenesStatistics
#'
calculateGenesCE <- function(objCOTAN) {
  # estimate Probabilities of 0 with internal function funProbZero
  probZero <- getProbabilityOfZero(objCOTAN)

  zeroOne <- as.matrix(getZeroOneProj(objCOTAN))

  feg <- getNumOfExpressingCells(objCOTAN) == getNumCells(objCOTAN)

  minusEntrM <- matrix(0.0, nrow = getNumGenes(objCOTAN),
                       ncol = getNumCells(objCOTAN))
  minusEntrM[!feg, ] <- (zeroOne[!feg, ] * log(1.0 - probZero[!feg, ])) +
                        ((1.0 - zeroOne[!feg, ]) * log(probZero[!feg, ]))

  return(set_names(-rowsums(minusEntrM, parallel = TRUE) /
                     getNumCells(objCOTAN), getGenes(objCOTAN)))
}


# local utility wrapper for parallel estimation of dispersion
runGDICalc <- function(genesBatches, S, topRows, cores) {

  ##  worker: operate on one batch of columns
  worker <- function(geneBatch, S, topRows) {
    tryCatch({
      ## slice columns (single read avoids copying the rest)
      subS <- as.matrix(S[, geneBatch, drop = FALSE])

      ## 1. sort each column descending
      subS <- colSort(subS, descending = TRUE)

      ## 2. keep only the top rows
      subS <- subS[seq_len(topRows), , drop = FALSE]

      ## 3. p-value of χ²₁ for each entry
      subS <- stats::pchisq(subS, df = 1L, lower.tail = FALSE)

      ## 4. mean p-value per column → transform
      gdi <- log(-log(colMeans(subS)))

      rm(subS)

      return(gdi)
    }, error = function(e) {
      structure(list(e), class = "try-error")
    })
  }

  if (cores != 1L) {
    ##  tiny sandbox env to keep each worker lean
    mini <- new.env(parent = baseenv())
    mini$colSort <- Rfast::colSort
    mini$colmeans <- Rfast::colmeans
    mini$logThis <- logThis
    mini$worker <- worker

    environment(mini$worker) <- mini

    ##  fork once, stream tasks
    res <- parallel::mclapply(genesBatches,
                              worker,
                              S = S,
                              topRows = topRows,
                              mc.cores  = cores,
                              mc.preschedule = FALSE)

    # spawned errors are stored as try-error classes
    resError <- unlist(lapply(res, inherits, "try-error"))
    if (any(resError)) {
      stop(res[[which(resError)[[1L]]]], call. = FALSE)
    }

    return(res)
  } else {
    res <- lapply(genesBatches,
                  worker,
                  S = S,
                  topRows = topRows)

    return(res)
  }
}


#' @details `calculateGDIGivenS()` produces a `vector` with the `GDI` for each
#'   column based on the `S` matrix (*Pearson's *\eqn{\chi^{2}}* test*)
#'
#' @param S a `matrix` object
#' @param rowsFraction The fraction of rows that will be averaged to calculate
#'   the `GDI`. Defaults to \eqn{5\%}
#' @param cores number of cores to use. Default is 1.
#' @param chunkSize number of elements to solve in batch in a single core.
#'   Default is 1024.
#'
#' @returns `calculateGDIGivenS()` returns a `vector` with the *GDI* data for
#'   each column of the input
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom stats pchisq
#'
#' @importFrom Rfast colSort
#' @importFrom Rfast colmeans
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @rdname GenesStatistics
#'

calculateGDIGivenS  <- function(S,
                                rowsFraction = 0.05,
                                cores        = 1L,
                                chunkSize    = 1024L) {
  # Beware S might not be square!
  assertthat::assert_that(length(dim(S)) == 2L,
                          rowsFraction > 0, rowsFraction <= 1)

  logThis("Calculate `GDI`: START", logLevel = 2L)

  cores <- handleMultiCore(cores)

  genes <- colnames(S)

  ## rows to retain per column
  topRows <- max(1L, round(nrow(S) * rowsFraction))

  ##  split genes into batches
  spIdx <- parallel::splitIndices(length(genes),
                                  ceiling(length(genes) / chunkSize))

  spGenes <- lapply(spIdx, function(x) genes[x])

  cores <- min(cores, length(spGenes))

  logThis(paste0("Executing ", length(spGenes), " genes batches"),
          logLevel = 3L)

  gdiList <- runGDICalc(genesBatches = spGenes, S = S,
                        topRows = topRows, cores = cores)

  ## glue result and reorder to original gene order
  GDI <- unlist(gdiList, use.names = FALSE, recursive = TRUE)

  names(GDI) <- genes

  logThis("Calculate `GDI`: DONE", logLevel = 2L)

  return(GDI)
}

#' @details `calculateGDIGivenCorr()` produces a `vector` with the `GDI` for
#'   each column based on the given correlation matrix, using the *Pearson's
#'   *\eqn{\chi^{2}}* test*
#'
#' @param corr a `matrix` object, possibly a subset of the columns of the full
#'   symmetric matrix
#' @param numDegreesOfFreedom a `int` that determines the number of degree of
#'   freedom to use in the \eqn{\chi^{2}} test
#' @param rowsFraction The fraction of rows that will be averaged to calculate
#'   the `GDI`. Defaults to \eqn{5\%}
#'
#' @returns `calculateGDIGivenCorr()` returns a `vector` with the *GDI* data for
#'   each column of the input
#'
#' @export
#'
#' @rdname GenesStatistics
#'

calculateGDIGivenCorr <-
  function(corr, numDegreesOfFreedom, rowsFraction = 0.05) {
  return(calculateGDIGivenS(corr^2L * numDegreesOfFreedom,
                            rowsFraction = rowsFraction))
}


#' @details `calculateGDI()` produces a `data.frame` with the *GDI* for each
#'   gene based on the `COEX` matrix
#'
#' @param objCOTAN a `COTAN` object
#' @param statType Which statistics to use to compute the p-values. By default
#'   it will use the `S` (*Pearson's *\eqn{\chi^{2}}* test*) otherwise the `G`
#'   (*G-test*)
#' @param rowsFraction The fraction of rows that will be averaged to calculate
#'   the `GDI`. Defaults to \eqn{5\%}
#' @param cores number of cores to use. Default is 1.
#' @param chunkSize number of elements to solve in batch in a single core.
#'   Default is 1024.
#'
#' @returns `calculateGDI()` returns a `data.frame` with:
#'  * `"sum.raw.norm"` the sum of the normalized data rows
#'  * `"GDI"` the *GDI* data
#'  * `"exp.cells"` the percentage of cells expressing the gene
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom tibble column_to_rownames
#'
#' @rdname GenesStatistics
#'

calculateGDI <- function(objCOTAN,
                         statType     = "S",
                         rowsFraction = 0.05,
                         cores        = 1L,
                         chunkSize    = 1024L) {
  logThis("Calculate GDI dataframe: START", logLevel = 2L)

  if (statType == "S") {
    logThis("Using S", logLevel = 3L)
    S <- calculateS(objCOTAN)
  } else if (statType == "G") {
    logThis("Using G", logLevel = 3L)
    S <- calculateG(objCOTAN)
  } else {
    stop("Unrecognised stat type: must be either 'S' or 'G'")
  }

  GDI <- calculateGDIGivenS(S, rowsFraction = rowsFraction,
                            cores = cores, chunkSize = chunkSize)
  GDI <- set_names(as.data.frame(GDI), "GDI")
  gc()

  sumRawNorm <- log(rowSums(getNuNormData(objCOTAN)))
  GDI <- merge(GDI, as.data.frame(list(sumRawNorm), col.names = "sum.raw.norm"),
               by = "row.names", all.x = TRUE)
  GDI <- column_to_rownames(GDI, var = "Row.names")
  rm(sumRawNorm)
  gc()

  expCells <- getNumOfExpressingCells(objCOTAN) / getNumCells(objCOTAN) * 100.0
  GDI <- merge(GDI, as.data.frame(list(expCells), col.names = "exp.cells"),
               by = "row.names", all.x = TRUE)
  GDI <- column_to_rownames(GDI, var = "Row.names")
  rm(expCells)
  gc()

  GDI <- GDI[, c("sum.raw.norm", "GDI", "exp.cells")]

  logThis("Calculate GDI dataframe: DONE", logLevel = 2L)

  return(GDI)
}


#' @details `calculatePValue()` computes the p-values for genes in the `COTAN`
#'   object. It can be used genome-wide or by setting some specific genes of
#'   interest. By default it computes the *p-values* using the `S` statistics
#'   (\eqn{\chi^{2}})
#'
#' @param objCOTAN a `COTAN` object
#' @param statType Which statistics to use to compute the p-values. By default
#'   it will use the "S" (Pearson's \eqn{\chi^{2}} test) otherwise the "G"
#'   (G-test)
#' @param geneSubsetCol an array of genes. It will be put in columns. If left
#'   empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows. If left empty
#'   the function will do it genome-wide.
#'
#' @return `calculatePValue()` returns a *p-value* `matrix` as `dspMatrix`
#'
#' @export
#'
#' @importFrom Matrix forceSymmetric
#' @importFrom Matrix pack
#'
#' @importClassesFrom Matrix dspMatrix
#'
#' @rdname GenesStatistics
#'
calculatePValue <- function(objCOTAN, statType = "S",
                            geneSubsetCol = vector(mode = "character"),
                            geneSubsetRow = vector(mode = "character")) {
  geneSubsetCol <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetCol)
  geneSubsetRow <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetRow)

  if (statType == "S") {
    logThis("Using S", logLevel = 3L)
    S <- calculateS(objCOTAN,
                    geneSubsetCol = geneSubsetCol,
                    geneSubsetRow = geneSubsetRow)
  } else if (statType == "G") {
    logThis("Using G", logLevel = 3L)
    S <- calculateG(objCOTAN,
                    geneSubsetCol = geneSubsetCol,
                    geneSubsetRow = geneSubsetRow)
  } else {
    stop("Unrecognised stat type: must be either 'S' or 'G'")
  }

  logThis("calculating PValues: START", logLevel = 2L)

  allCols <- identical(geneSubsetCol, getGenes(objCOTAN))
  allRows <- identical(geneSubsetRow, getGenes(objCOTAN))

  strCol <- if (allCols) "genome wide" else "on a set of genes"
  strRow <- if (allRows) "genome wide" else "on a set of genes"
  logThis(paste("Get p-values", strCol, "on columns and", strRow, "on rows"),
          logLevel = 2L)

  pValues <- pchisq(as.matrix(S), df = 1L, lower.tail = FALSE)

  if (allCols && allRows) {
    pValues <- pack(forceSymmetric(pValues))
  }

  logThis("calculating PValues: DONE", logLevel = 2L)

  return(pValues)
}


#' @details `calculatePDI()` computes the p-values for genes in the `COTAN`
#'   object using [calculatePValue()] and takes their
#'   \eqn{\log{({-\log{(\cdot)}})}} to calculate the genes' *Pair Differential
#'   Index*
#'
#' @param objCOTAN a `COTAN` object
#' @param statType Which statistics to use to compute the p-values. By default
#'   it will use the "S" (Pearson's \eqn{\chi^{2}} test) otherwise the "G"
#'   (G-test)
#' @param geneSubsetCol an array of genes. It will be put in columns. If left
#'   empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows. If left empty
#'   the function will do it genome-wide.
#'
#' @return `calculatePDI()` returns a *Pair Differential Index* `matrix` as
#'   `dspMatrix`
#'
#' @export
#'
#' @rdname GenesStatistics
#'
calculatePDI <- function(objCOTAN, statType = "S",
                         geneSubsetCol = vector(mode = "character"),
                         geneSubsetRow = vector(mode = "character")) {
  pValues <- calculatePValue(objCOTAN, statType = statType,
                             geneSubsetCol = geneSubsetCol,
                             geneSubsetRow = geneSubsetRow)
  return(log(-log(pValues)))
}
