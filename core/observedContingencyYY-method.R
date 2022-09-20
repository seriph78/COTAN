#' calculate observed yes/yes field of contingency table
setMethod(
  "observedContingencyYY", 
  "COTAN",
  function(objCOTAN, cells) {
    zeroOne <- as.matrix(objCOTAN@raw)
    # zeroOne matrix: formed by row data matrix changed to 0-1 matrix
    zeroOne[zeroOne > 0] <- 1
    
    if(cells){
      # for cells 
      YY <-  t(as.matrix(zeroOne)) %*% as.matrix(zeroOne)
    } else{
      # for genes
      YY <- as.matrix(zeroOne) %*% t(as.matrix(zeroOne))
    }
    
    YY <- as(as.matrix(YY), "sparseMatrix")
    
    return(YY)
  }
)