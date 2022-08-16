# linear estimator of lambda
#' @export
setMethod(
  "estimateLambdaLinear",
  "COTAN",
  function(objCOTAN) {
    objCOTAN@lambda <- Matrix::rowMeans(objCOTAN@raw, dims=1, na.rm=TRUE)

    return(objCOTAN)
  }
)
