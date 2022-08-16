# linear estimator of nu
#' @export
setMethod(
  "estimateNuLinear",
  "COTAN",
  function(objCOTAN) {
    # average of all raw values
    global_mean <- Matrix::mean(objCOTAN@raw, na.rm = TRUE)

    # raw column averages divided by global_mean
    objCOTAN@nu <- Matrix::colMeans(objCOTAN@raw, dims = 1, na.rm = TRUE)

    objCOTAN@nu <- objCOTAN@nu / global_mean

    return(objCOTAN)
  }
)
