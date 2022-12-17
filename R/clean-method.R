#' clean
#'
#' @description Main function that can be used to check and clean the dataset.
#'   It also produces and stores the estimators for nu and lambda
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns The updated `COTAN` object
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#'
#' @rdname clean
#'
setMethod(
  "clean",
  "COTAN",
  function(objCOTAN) {

    # We want to discard genes having less than 3 non-zero counts per 1000 cells
    threshold <- round(getNumCells(objCOTAN) * 3 / 1000, digits = 0)
    genesToDrop <- getGenes(objCOTAN)[rowSums(getZeroOneProj(objCOTAN)) <= threshold]

    if (!is_empty(genesToDrop)) {
      objCOTAN <- dropGenesCells(objCOTAN, genes = genesToDrop)
    }

    # We want to discard cells having less than 2 non-zero counts per 1000 genes
    threshold <- round(getNumGenes(objCOTAN) * 2 / 1000, digits = 0)
    cellsToDrop <- getCells(objCOTAN)[colSums(getZeroOneProj(objCOTAN)) <= threshold]

    if (!is_empty(cellsToDrop)) {
      objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)
    }

    logThis(paste0("Genes/cells selection done:",
                   " dropped [", length(genesToDrop), "] genes",
                   " and [", length(cellsToDrop), "] cells"),
            logLevel = 1)

    logThis(paste0("Working on [", getNumGenes(objCOTAN), "]",
                   " genes and [", getNumCells(objCOTAN), "] cells"),
            logLevel = 2)

    objCOTAN <- estimateLambdaLinear(objCOTAN)
    objCOTAN <- estimateNuLinear(objCOTAN)
    objCOTAN <- findHousekeepingGenes(objCOTAN)
    objCOTAN <- findFullyExpressedCells(objCOTAN)

    gc()

    return(objCOTAN)
  }
)
