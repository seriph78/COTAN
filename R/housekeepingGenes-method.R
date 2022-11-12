#' determines the housekeepingGenes vector of a COTAN object
#' @export
setMethod(
  "housekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    nCells <- getNumCells(objCOTAN)

    if (nCells == 0) {
      stop("zero cells given")
    }

    # determine positive UMI
    cells <- getZeroOneProj(objCOTAN)
    
    # name of the genes with positive UMI count in every single cell
    objCOTAN@hkGenes <- names(which(rowSums(cells) == nCells))

    return(objCOTAN)
  }
)
