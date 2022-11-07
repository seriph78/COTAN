#' getNu
#' 
#' This function extract the nu array.
#'
#' @param object A COTAN object
#'
#' @return the nu array.
#' @export
#' @rdname getNu
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

#' getLambda
#' 
#' This function extract the lambda array (mean expression for each gene).
#'
#' @param object A COTAN object
#'
#' @return the lambda array.
#' @export
#' @rdname getLambda
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

#' getDispersion
#' 
#' This function extract the a array.
#'
#' @param object A COTAN object
#'
#' @return the a array.
#' @export
#' @rdname getDispersion
setMethod(
  "getDispersion",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@dispersion)) {
      warning("dispersion is empty")
    }
    
    return(objCOTAN@dispersion)
  }
)

#' getHousekeepingGenes
#' 
#' This function return the genes expressed in all cells in the dataset.
#'
#' @param object A COTAN object
#'
#' @return an array containing all genes expressed in all cells
#' @export
#' @rdname getHousekeepingGenes
setMethod(
  "getHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@hkGenes)) {
      warning("hkGenes is empty")
    }
    
    return(objCOTAN@hkGenes)
  }
)


#' getGenes
#'
#' This function extract all genes in the dataset.
#'
#' @param object A COTAN object
#'
#' @return a gene array
#' @export
#' @rdname getGenes
setMethod(
  "getGenes", "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@raw)) {
      warning("raw is empty")
    }
    
    return(rownames(objCOTAN@raw))
  }
)


