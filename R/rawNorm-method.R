#' rawNorm
#'
#' This function initializes the normalized count table.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the updated COTAN object
#' @importFrom rlang is_empty
#' @importFrom Matrix t
#' @export
#' @rdname rawNorm
setMethod(
  "rawNorm",
  "COTAN",
  function(objCOTAN) {
    if(is_empty(objCOTAN@raw)) {
      stop("empty raw")
    }

    if (is_empty(objCOTAN@nu)) {
      stop("nu must not be empty, estimate it")
    }

    objCOTAN@rawNorm <- t(t(objCOTAN@raw) * (1/(as.vector(objCOTAN@nu))))

    return(objCOTAN)
  }
)
