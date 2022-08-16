#' @export
setMethod(
  "rawNorm",
  "COTAN",
  function(objCOTAN) {
    if(is_empty(objCOTAN@nu)){
      stop("nu must not be empty, estimate it")
    }
    
    objCOTAN@rawNorm <- objCOTAN@raw / objCOTAN@nu
    
    return(objCOTAN)
  }
)
