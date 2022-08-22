# determines the houseKeepingGenes vector of a COTAN object
#' @export

setMethod(
  "houseKeepingGenes",
  "COTAN",
  function(objCOTAN) {
    cells <- as.matrix(objCOTAN@raw)

    # determine positive UMI
    cells[cells > 0] <- 1
    #cells[cells <= 0] <- 0 as cells cannot be negative

    # name of the genes with positive UMI count in every single cell
    if(is.na(objCOTAN@nCells) | objCOTAN@nCells == 0 ) {
      stop("nCells not initialized or zero")
    }
    
    objCOTAN@hKGenes <- names(which(rowSums(cells) == objCOTAN@nCells))
    
    return(objCOTAN)
  }
)
