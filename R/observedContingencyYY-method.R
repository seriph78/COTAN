
#' calculate observed yes/yes field of contingency table
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells, 
#' otherwise for the genes
#' @importFrom Matrix t
#' @importClassesFrom Matrix symmetricMatrix
#' @rdname observedContingencyYY
setMethod(
  "observedContingencyYY", 
  "COTAN",
  function(objCOTAN, cells) {
    zeroOne <- getZeroOneProj(objCOTAN)

    if(isTRUE(cells)){
      # for cells 
      YY <- t(zeroOne) %*% zeroOne
    } else{
      # for genes
      YY <- zeroOne %*% t(zeroOne)
    }

    return(as(YY, "symmetricMatrix"))
  }
)