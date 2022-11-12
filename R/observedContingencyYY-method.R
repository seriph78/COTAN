#' calculate observed yes/yes field of contingency table
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells, 
#' otherwise for the genes
#' @importFrom Matrix t
setMethod(
  "observedContingencyYY", 
  "COTAN",
  function(objCOTAN, cells) {
    zeroOne <- as.matrix(objCOTAN@raw)
    # zeroOne matrix: formed by row data matrix changed to 0-1 matrix
    zeroOne[zeroOne > 0] <- 1
    zeroOne <- as.matrix(zeroOne)
    
    if(isTRUE(cells)){
      # for cells 
      YY <-  t(zeroOne) %*% zeroOne
    } else{
      # for genes
      YY <- zeroOne %*% t(zeroOne)
    }

    YY <- as(as.matrix(YY), "sparseMatrix")
    
    return(YY)
  }
)