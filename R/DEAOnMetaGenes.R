#'
#' @details `DEAOnMetaGenes()` is used to run the *Differential Expression
#'   Analysis* using the `COTAN` contingency tables on each given *meta-gene*
#'
#' @param objCOTAN a `COTAN` object
#' @param metaGenes a `list` of *meta-genes*, such as the result of
#'   [defineMetaGenes()].
#'
#' @return `DEAOnMetaGenes()` returns a `list` with two objects:
#'   * "coex"    - the co-expression `data.frame` for the cells against each
#'                 *meta-gene*
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
#' @examples
#' # pca <- pcaGenesCoex()
#' # metaGenes <- defineMetaGenes(pca, )
#' # metaGenesPlot(pca[, 1L], pca[, 2L], metaGenes)
#' # DEAOnMetaGenes()
#'
#' @rdname MetaGenes
#'
DEAOnMetaGenes <- function(objCOTAN, metaGenes) {


#  assert_that(length(clusters) == getNumCells(objCOTAN),
#              msg = "Passed/retrieved clusterization has the wrong size")

  assert_that(all(unlist(metaGenes) %in% getGenes(objCOTAN)),
              msg = paste("Some genes in the meta-genes",
                          "are not part of the 'COTAN' object"))

  logThis("Differential Expression Analysis - START", logLevel = 2L)

  zeroOne <- getZeroOneProj(objCOTAN)

  probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

  numGenes <- getNumGenes(objCOTAN)

  metaGenesCoex <- data.frame()
  metaGenesPVal <- data.frame()

  for (center in names(metaGenes)) {
    gc()

    mgName <- paste0("MG_", center)

    logThis("*", appendLF = FALSE, logLevel = 1L)
    logThis(paste0(" analysis of meta-gene: '", mgName, "' - START"),
            logLevel = 3L)

    genesIn <- getGenes(objCOTAN) %in% metaGenes[[center]]

    numGenesIn  <- sum(genesIn)
    numGenesOut <- numGenes - numGenesIn

    if (numGenesIn == 0L) {
      warning("Meta-gene '", center, "' has no genes assigned to it!")
    }

    observedYI <- colSums(zeroOne[genesIn, ])

    expectedNI <- colSums(probZero[ genesIn, ])
    expectedNO <- colSums(probZero[!genesIn, ])
    expectedYI <- numGenesIn  - expectedNI
    expectedYO <- numGenesOut - expectedNO

    coex <- (observedYI  - expectedYI) / sqrt(numGenes) *
              sqrt(1.0 / pmax(1.0, expectedNI) +
                   1.0 / pmax(1.0, expectedNO) +
                   1.0 / pmax(1.0, expectedYI) +
                   1.0 / pmax(1.0, expectedYO))

    pValue <- pchisq(coex^2L * numGenes, df = 1L, lower.tail = FALSE)

    if (anyNA(pValue)) {
      warning("Got some NA in the p-value",
              toString(which(is.na(pValue), arr.ind = TRUE)))
    }

    rm(expectedYO, expectedYI, expectedNO, expectedNI)
    rm(observedYI)
    gc()

    metaGenesCoex <- setColumnInDF(metaGenesCoex, colToSet = coex,
                                   colName = mgName,
                                   rowNames = colnames(zeroOne))
    metaGenesPVal <- setColumnInDF(metaGenesPVal, colToSet = pValue,
                                   colName = mgName,
                                   rowNames = colnames(zeroOne))

    logThis(paste0("* analysis of meta-gene: '", mgName, "' - DONE"),
            logLevel = 3L)
  }
  logThis("", logLevel = 1L)

  logThis("Differential Expression Analysis - DONE", logLevel = 2L)

  return(list("coex" = metaGenesCoex, "p-value" = metaGenesPVal))
}
