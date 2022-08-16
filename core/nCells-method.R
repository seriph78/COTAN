#' @export
setMethod(
  "nCells",
  "COTAN",
  function(objCOTAN){
    objCOTAN@nCells <- length(rownames(objCOTAN@raw))
  }
  
  return(objCOTAN)
)