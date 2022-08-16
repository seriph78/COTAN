# determines the houseKeepingGenes vector of a COTAN object
#' @export

setMethod(
  "houseKeepingGenes",
  "COTAN",
  function(objCOTAN) {
    cells <- as.matrix(objCOTAN@raw)

    # determine positive UMI
    cells[cells > 0] <- 1
    cells[cells <= 0] <- 0

    # name of the genes with positive UMI count in every single cell
    objCOTAN@hKGenes <- names(which(rowSums(cells) == length(colnames(cells))))
    
    return(objCOTAN)
  }
)
