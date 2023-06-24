
#' Calculations of genes statistics
#'
#' @description A collectino fo functions returning various statistics
#'   associated to the genes. In particular the *discrepancy* between the
#'   expected probabilities of zero and their actual occurrences, both at single
#'   gene level or looking at genes' pairs
#'
#' @name GenesStatistics
NULL


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
#' @export
#'
#' @rdname GenesStatistics
#'
calculateGenesCE <- function(objCOTAN) {
  # estimate Probabilities of 0 with internal function funProbZero
  probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

  zeroOne <- as.matrix(getZeroOneProj(objCOTAN))

  feg <- rowsums(zeroOne) == getNumCells(objCOTAN)

  minusEntrM <- matrix(0.0, nrow = getNumGenes(objCOTAN),
                       ncol = getNumCells(objCOTAN))
  minusEntrM[!feg, ] <- (zeroOne[!feg, ] * log(1.0 - probZero[!feg, ])) +
                        ((1.0 - zeroOne[!feg, ]) * log(probZero[!feg, ]))

  return(set_names(-rowMeans(minusEntrM), getGenes(objCOTAN)))
}



#' @details `calculateGDI()` produces a `data.frame` with the *GDI* for each
#'   gene based on the `COEX` matrix
#'
#' @param objCOTAN a `COTAN` object
#' @param statType Which statistics to use to compute the p-values. By default
#'   it will use the `S` (*Pearson's *\eqn{\chi^{2}}* test*) otherwise the `G`
#'   (*G-test*)
#'
#' @returns `calculateGDI()` returns a `data.frame` with the *GDI* data
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom stats pchisq
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#'
#' @rdname GenesStatistics
#'
calculateGDI <- function(objCOTAN, statType = "S") {
  if (statType == "S") {
    logThis("Using S", logLevel = 3L)
    S <- calculateS(objCOTAN)
  } else if (statType == "G") {
    logThis("Using G", logLevel = 3L)
    S <- calculateG(objCOTAN)
  } else {
    stop("Unrecognised stat type: must be either 'S' or 'G'")
  }

  logThis("Calculate GDI dataframe: START", logLevel = 2L)

  top5pcRows <- as.integer(max(1L:round(getNumGenes(objCOTAN) / 20.0,
                                        digits = 0L)))

  pValueSorted <- apply(S, 2L, sort, decreasing = TRUE)
  pValueSorted <- pValueSorted[1L:top5pcRows, , drop = FALSE]
  pValueSorted <- pchisq(as.matrix(pValueSorted), df = 1L, lower.tail = FALSE)

  GDI <- log(-log(colMeans(pValueSorted)))
  GDI <- set_names(as.data.frame(GDI), "GDI")
  rm(pValueSorted)
  gc()

  sumRawNorm <- log(rowSums(getNormalizedData(objCOTAN)))
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
#' @return `calculatePDI()` returns a *PDI* `matrix` as `dspMatrix`
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
