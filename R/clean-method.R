#' @aliases clean
#'
#' @details `clean()` is the main method that can be used to check and clean the
#'   dataset. It will discard any genes that has less than 3 non-zero counts per
#'   thousand cells and all cells expressing less than 2 per thousand genes. It
#'   also produces and stores the estimators for nu and lambda
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `clean()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "clean",
  "COTAN",
  function(objCOTAN,
           genesLowThreshold = 0.003, cellsLowThreshold = 0.002,
           genesHighThreshold = 0.99, cellsHighThreshold = 0.99) {

    # We want to discard genes having less than given threshold
    # default: less than 3 non-zero counts per 1000 cells
    threshold <- round(getNumCells(objCOTAN) * genesThreshold, digits = 0L)
    genesToDrop <-
      getGenes(objCOTAN)[rowSums(getZeroOneProj(objCOTAN)) <= threshold]

    if (!is_empty(genesToDrop)) {
      objCOTAN <- dropGenesCells(objCOTAN, genes = genesToDrop)
    }

    # We want to discard cells having less than given threshold
    # default: less than 2 non-zero counts per 1000 genes
    threshold <- round(getNumGenes(objCOTAN) * cellsThreshold, digits = 0L)
    cellsToDrop <-
      getCells(objCOTAN)[colSums(getZeroOneProj(objCOTAN)) <= threshold]

    if (!is_empty(cellsToDrop)) {
      objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)
    }

    logThis(paste0("Genes/cells selection done:",
                   " dropped [", length(genesToDrop), "] genes",
                   " and [", length(cellsToDrop), "] cells"),
            logLevel = 1L)

    logThis(paste0("Working on [", getNumGenes(objCOTAN), "]",
                   " genes and [", getNumCells(objCOTAN), "] cells"),
            logLevel = 2L)

    objCOTAN <- estimateLambdaLinear(objCOTAN)
    objCOTAN <- estimateNuLinear(objCOTAN)
    objCOTAN <- findFullyExpressedGenes(objCOTAN)
    objCOTAN <- findFullyExpressingCells(objCOTAN)

    gc()

    return(objCOTAN)
  }
)
