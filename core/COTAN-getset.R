#' @export
setMethod(
  "getLambda",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@lambda)) {
      warning("lambda is empty")
    }

    return(objCOTAN@lambda)
  }
)

#' @export
setMethod(
  "getNu",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@nu)) {
      warning("nu is empty")
    }

    return(objCOTAN@nu)
  }
)
