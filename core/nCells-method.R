#' @export
setMethod(
  "nCells",
  "COTAN",
  function(objCOTAN){
    objCOTAN@nCells <- length(colnames(objCOTAN@raw))
    
    return(objCOTAN)
  }
)