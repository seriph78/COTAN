# estimate vector mu
setMethod(
  "estimateMu",
  "COTAN",
  function(objCOTAN){
    muEstimator <- objCOTAN@lambda %*% t(objCOTAN@nu)
    rownames(muEstimator) <- rownames(objCOTAN@raw)
    return(muEstimator)    
  }
)