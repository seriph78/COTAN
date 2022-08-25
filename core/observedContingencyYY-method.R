#' calculate observed yes/yes contingency table
setMethod(
  "observedContingencyYY", 
  "COTAN",
  function(objCOTAN) {
    cells <- as.matrix(objCOTAN@raw)
    # Cells matrix : formed by row data matrix changed to 0-1 matrix
    cells[cells > 0] <- 1
    
    YY <- as.matrix(cells) %*% t(as.matrix(cells))
    YY <- as(as.matrix(YY), "sparseMatrix")
    return(YY)
  }
)