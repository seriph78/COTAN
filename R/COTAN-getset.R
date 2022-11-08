
#' getRawData
#'
#' This function extract the raw count table.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the raw count dataframe
#' @export
#' @rdname getRawData
setMethod(
  "getRawData",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@raw)) {
      warning("raw is empty")
    }
    
    return(objCOTAN@raw)
  }
)

#' getNormalizedData
#'
#' This function extract the normalized count table.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the normalized count dataframe (divided by nu).
#' @export
#' @rdname getNormalizedData
setMethod(
  "getNormalizedData",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@rawNorm)) {
      warning("rawNorm is empty")
    }
    
    return(objCOTAN@rawNorm)
  }
)

#' getNu
#' 
#' This function extract the nu array.
#'
#' @param objCOTAN A COTAN object
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
#' @param objCOTAN A COTAN object
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
#' @param objCOTAN A COTAN object
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
#' @param objCOTAN A COTAN object
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
#' @param objCOTAN A COTAN object
#'
#' @return a gene array
#' @export
#' @rdname getGenes
setMethod(
  "getGenes",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@raw)) {
      warning("raw is empty")
    }
    
    return(rownames(objCOTAN@raw))
  }
)


#' getMetadataDataset
#'
#' This function extract the meta-data stored for the dataset.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the meta-data dataframe
#' @export
#' @rdname getMetadataDataset
setMethod(
  "getMetadataDataset",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@raw)) {
      warning("raw is empty")
    }
    
    meta <- object@meta
    return(meta)
  }
)

