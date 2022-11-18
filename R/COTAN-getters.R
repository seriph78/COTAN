# COTAN objects (read-only) accessors

#' getRawData
#'
#' This function extract the raw count table.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the raw count dataframe
#' @importFrom rlang is_empty
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

#' getNumCells
#'
#' This function extracts the number of cells in the sample (𝑚)
#'
#' @param objCOTAN A COTAN object
#'
#' @return the number of cells in the sample (𝑚).
#' @importFrom rlang is_empty
#' @export
#' @rdname getNumCells
setMethod(
  "getNumCells",
  "COTAN",
  function(objCOTAN){
    return(ncol(objCOTAN@raw))
  }
)

#' getNumGenes
#'
#' This function extracts the number of genes in the sample (𝑛)
#'
#' @param objCOTAN A COTAN object
#'
#' @return the number of genes in the sample (𝑛).
#' @export
#' @rdname getNumGenes
setMethod(
  "getNumGenes",
  "COTAN",
  function(objCOTAN){
    return(nrow(objCOTAN@raw))
  }
)


#' getCells
#'
#' This function extract all cells in the dataset.
#'
#' @param objCOTAN A COTAN object
#'
#' @return a cell array
#' @importFrom rlang is_empty
#' @export
#' @rdname getCells
setMethod(
  "getCells",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@raw)) {
      warning("raw is empty")
    }

    return(colnames(objCOTAN@raw))
  }
)


#' getGenes
#'
#' This function extract all genes in the dataset.
#'
#' @param objCOTAN A COTAN object
#'
#' @return a gene array
#' @importFrom rlang is_empty
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


#' getZeroOneProj
#'
#' This function extract the raw count table where any
#' positive number has been replaced with 1
#'
#' @param objCOTAN A COTAN object
#'
#' @return the raw count projected to 0/1
#' @importFrom rlang is_empty
#' @export
#' @rdname getZeroOneProj
setMethod(
  "getZeroOneProj",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@raw)) {
      warning("raw is empty")
    }

    return(sign(objCOTAN@raw))
  }
)


#' getCellsSize
#'
#' This function extracts the cell raw library size.
#'
#' @param object A COTAN object
#'
#' @return an array with the library sizes
#' @importFrom rlang is_empty
#' @export
#' @rdname getCellsSize
setMethod(
  "getCellsSize", "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@raw)) {
      warning("raw is empty")
    }

    return(colSums(objCOTAN@raw))
  }
)


#' getNormalizedData
#'
#' This function extracts the normalized count table.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the normalized count dataframe (divided by nu).
#' @importFrom rlang is_empty
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
#' @importFrom rlang is_empty
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
#' @importFrom rlang is_empty
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
#' @importFrom rlang is_empty
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
#' @importFrom rlang is_empty
#' @export
#' @rdname getHousekeepingGenes
setMethod(
  "getHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@hkGenes)
  }
)


#' flagNotHousekeepingGenes
#'
#' This function returns a Boolean vector with TRUE for those genes that are
#' not housekeeping.
#'
#' @param objCOTAN A COTAN object
#'
#' @return an array of Booleans with TRUE for genes that are not housekeeping
#'
#' @importFrom rlang is_empty
#'
#' @export
#' @rdname flagNotHousekeepingGenes
setMethod(
  "flagNotHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    return( !(rownames(objCOTAN@raw) %in% objCOTAN@hkGenes ) )
  }
)


#' getMetadataDataset
#'
#' This function extract the meta-data stored for the data-set.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the meta-data data.frame
#' @importFrom rlang is_empty
#' @export
#' @rdname getMetadataDataset
setMethod(
  "getMetadataDataset",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(objCOTAN@metaDataset)) {
      warning("metaDataset is empty")
    }

    return(objCOTAN@metaDataset)
  }
)


#' getCoex
#'
#' This function extract a complete or a partial coex matrix
#' from the COTAN object.
#'
#' @param object A COTAN object
#' @param asMatrix A Boolean to query whether return the coex as a matrix
#' or leave it in which-ever form it already is
#' @param genes A vector of gene names. It will exclude any gene not on the list.
#' By defaults the function will keep all genes.
#'
#' @return the coex values
#'
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' coex <- getCoex(ERCC.cotan, asMatrix = TRUE)
#'
#' @rdname getCoex
setMethod(
  "getCoex",
  "COTAN",
  function(objCOTAN, asMatrix = TRUE, genes = "all") {
    if (!asMatrix) {
      stopifnot("Get coex as list. Genes subset not supported yet." = genes == "all")
      return(objCOTAN@coex)
    }
    else {
      return( vec2mat_rfast(objCOTAN@coex, genes = genes) )
    }
  }
)


#' getCellsCoex
#'
#' This function extract a complete or a partial cellsCoex matrix
#' from the COTAN object.
#'
#' @param object A COTAN object
#' @param asMatrix A Boolean to query whether return the cellsCoex as a matrix
#' or leave it in which-ever form it already is
#' @param cells A vector of cell names. It will exclude any cell not on the list.
#' By defaults the function will keep all cells.
#'
#' @return the cellsCoex values
#'
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' coex <- getCellsCoex(ERCC.cotan, asMatrix = TRUE)
#'
#' @rdname getCellsCoex
setMethod(
  "getCellsCoex",
  "COTAN",
  function(objCOTAN, asMatrix = TRUE, cells = "all") {
    if (!asMatrix) {
      stopifnot("Get coex as list. Cells subset not supported yet." = cells == "all")
      return(objCOTAN@coex)
    }
    else {
      return( vec2mat_rfast(objCOTAN@cellsCoex, genes = cells) )
    }
  }
)

