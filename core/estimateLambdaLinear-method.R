# linear estimator of lambda
#' @param objCOTAN a COTAN object
#' @export
setMethod(
  "estimateLambdaLinear",
  "COTAN",
  function(objCOTAN) {
    objCOTAN@lambda <- Matrix::rowMeans(objCOTAN@raw, dims=1, na.rm=TRUE)

    return(objCOTAN)
  }
)
