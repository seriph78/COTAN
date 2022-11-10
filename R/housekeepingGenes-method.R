#' determines the housekeepingGenes vector of a COTAN object
#' @export
setMethod(
  "housekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    cells <- as.matrix(objCOTAN@raw)

    # determine positive UMI
    cells[cells > 0] <- 1
    #cells[cells <= 0] <- 0 as cells cannot be negative

    # name of the genes with positive UMI count in every single cell
    nCells <- getNumCells(objCOTAN)
    if (nCells == 0) {
      stop("zero cells given")
    }
    
    objCOTAN@hkGenes <- names(which(rowSums(cells) == nCells))

    return(objCOTAN)
  }
)
