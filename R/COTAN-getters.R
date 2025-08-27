# ------- `COTAN` raw data accessors --------

#' @title Raw data `COTAN` accessors
#'
#' @description These methods extract information out of a just created `COTAN`
#'   object. The accessors have **read-only** access to the object.
#'
#' @name RawDataGetters
NULL

#' @aliases getRawData
#'
#' @details `getRawData()` extracts the raw count table.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getRawData()` returns the raw count sparse matrix
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#'
#' rawData <- getRawData(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getRawData",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@raw)
  }
)


#' @aliases getNumCells
#'
#' @details `getNumCells()` extracts the number of cells in the sample (\eqn{m})
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNumCells()` returns the number of cells in the sample (\eqn{m})
#'
#' @export
#'
#' @examples
#' numCells <- getNumCells(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getNumCells",
  "COTAN",
  function(objCOTAN) {
    return(ncol(objCOTAN@raw))
  }
)


#' @aliases getNumGenes
#'
#' @details `getNumGenes()` extracts the number of genes in the sample (\eqn{n})
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNumGenes()` returns the number of genes in the sample (\eqn{n})
#'
#' @export
#'
#' @examples
#' numGenes <- getNumGenes(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getNumGenes",
  "COTAN",
  function(objCOTAN) {
    return(nrow(objCOTAN@raw))
  }
)


#' @aliases getCells
#'
#' @details `getCells()` extract all cells in the dataset.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getCells()` returns a character array with the cells' names
#'
#' @export
#'
#' @examples
#' cellsNames <- getCells(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getCells",
  "COTAN",
  function(objCOTAN) {
    return(colnames(objCOTAN@raw))
  }
)


#' @aliases getGenes
#'
#' @details `getGenes()` extract all genes in the dataset.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getGenes()` returns a character array with the genes' names
#'
#' @export
#'
#' @examples
#' genesNames <- getGenes(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getGenes",
  "COTAN",
  function(objCOTAN) {
    return(rownames(objCOTAN@raw))
  }
)


#' @aliases getZeroOneProj
#'
#' @details `getZeroOneProj()` extracts the raw count table where any
#'   positive number has been replaced with `1`
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getZeroOneProj()` returns the raw count matrix projected to `0` or
#'   `1`
#'
#' @export
#'
#' @examples
#' zeroOne <- getZeroOneProj(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getZeroOneProj",
  "COTAN",
  function(objCOTAN) {
    return(sign(objCOTAN@raw))
  }
)


#' @aliases getCellsSize
#'
#' @details `getCellsSize()` extracts the cell raw library size.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @return `getCellsSize()` returns an array with the library sizes
#'
#' @export
#'
#' @importFrom Matrix colSums
#'
#' @examples
#' cellsSize <- getCellsSize(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getCellsSize",
  "COTAN",
  function(objCOTAN) {
    return(colSums(objCOTAN@raw))
  }
)


#' @aliases getNumExpressedGenes
#'
#' @details `getNumExpressedGenes()` extracts the number of genes expressed for
#'   each cell. Exploits a feature of [Matrix::CsparseMatrix-class]
#'
#' @param objCOTAN a `COTAN` object
#'
#' @return `getNumExpressedGenes()` returns an array with the library sizes
#'
#' @export
#'
#' @examples
#' numExpGenes <- getNumExpressedGenes(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getNumExpressedGenes",
  "COTAN",
  function(objCOTAN) {
    # We expliot CsparseMatrix feature
    ret <- diff(objCOTAN@raw@p)
    names(ret) <- getCells(objCOTAN)
    return(ret)
  }
)


#' @aliases getGenesSize
#'
#' @details `getGenesSize()` extracts the genes raw library size.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @return `getGenesSize()` returns an array with the library sizes
#'
#' @export
#'
#' @importFrom Matrix rowSums
#'
#' @examples
#' genesSize <- getGenesSize(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getGenesSize",
  "COTAN",
  function(objCOTAN) {
    return(rowSums(objCOTAN@raw))
  }
)


#' @aliases getNumOfExpressingCells
#'
#' @details `getNumOfExpressingCells()` extracts, for each gene, the number of
#'   cells that are expressing it. Exploits a feature of
#'   [Matrix::CsparseMatrix-class]
#'
#' @param objCOTAN a `COTAN` object
#'
#' @return `getNumOfExpressingCells()` returns an array with the library sizes
#'
#' @export
#'
#' @examples
#' numExpCells <- getNumOfExpressingCells(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getNumOfExpressingCells",
  "COTAN",
  function(objCOTAN) {
    # We exploit CsparseMatrix feature
    ret <- diff(t(objCOTAN@raw)@p)
    names(ret) <- getGenes(objCOTAN)
    return(ret)
  }
)


# ---------- Handle meta-data --------

#' @title Handling *meta-data* in `COTAN` objects
#'
#' @description Much of the information stored in the `COTAN` object is
#'   compacted into three `data.frame`s:
#'   * `"metaDataset"` - contains all general information about the data-set
#'   * `"metaGenes"` - contains genes' related information along the `lambda`
#'     and `dispersion` vectors and the fully-expressed flag
#'   * `"metaCells"` - contains cells' related information along the `nu`
#'     vector, the fully-expressing flag, the *clusterizations* and the
#'     *conditions*
#'
#' @name HandleMetaData
NULL

#' @aliases getMetadataDataset
#'
#' @details `getMetadataDataset()` extracts the meta-data stored for the
#'   current data-set.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getMetadataDataset()` returns the meta-data `data.frame`
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#'
#' objCOTAN <- initializeMetaDataset(objCOTAN, GEO = "test_GEO",
#'                                   sequencingMethod = "distribution_sampling",
#'                                   sampleCondition = "reconstructed_dataset")
#'
#' objCOTAN <- addElementToMetaDataset(objCOTAN, "Test",
#'                                     c("These are ", "some values"))
#'
#' dataSetInfo <- getMetadataDataset(objCOTAN)
#'
#' @rdname HandleMetaData
#'
setMethod(
  "getMetadataDataset",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaDataset)
  }
)


#' @aliases getMetadataElement
#'
#' @details `getMetadataElement()` extracts the value associated with the
#'   given tag if present or an empty string otherwise.
#'
#' @param objCOTAN a `COTAN` object
#' @param tag The tag associated to the wanted value
#'
#' @returns `getMetadataElement()` returns a string with the relevant value
#'
#' @export
#'
#' @examples
#' numInitialCells <- getMetadataElement(objCOTAN, datasetTags()[["cells"]])
#'
#' @rdname HandleMetaData
#'
setMethod(
  "getMetadataElement",
  "COTAN",
  function(objCOTAN, tag) {
    meta <- getMetadataDataset(objCOTAN)
    rowPos <- getMetaInfoRow(meta, tag)
    if (rowPos == 0L) {
      return("")
    } else {
      return(meta[, -1L, drop = FALSE][rowPos, ])
    }
  }
)


#' @aliases getMetadataGenes
#'
#' @details `getMetadataGenes()` extracts the meta-data stored for the genes
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getMetadataGenes()` returns the genes' meta-data `data.frame`
#'
#' @export
#'
#' @examples
#' metaGenes <- getMetadataGenes(objCOTAN)
#'
#' @rdname HandleMetaData
#'
setMethod(
  "getMetadataGenes",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaGenes)
  }
)


#' @aliases getMetadataCells
#'
#' @details `getMetadataCells()` extracts the meta-data stored for the cells
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getMetadataCells()` returns the cells' meta-data `data.frame`
#'
#' @export
#'
#' @examples
#' metaCells <- getMetadataCells(objCOTAN)
#'
#' @rdname HandleMetaData
#'
setMethod(
  "getMetadataCells",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaCells)
  }
)


#' @aliases getDims
#'
#' @details `getDims()` extracts the sizes of all slots of the `COTAN` object
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getDims()` returns a named `list` with the sizes of the slots
#'
#' @export
#'
#' @examples
#' allSizes <- getDims(objCOTAN)
#'
#' @rdname HandleMetaData
#'
setMethod(
  "getDims",
  "COTAN",
  function(objCOTAN) {
    return(list("raw"          = dim(objCOTAN@raw),
                "genesCoex"    = dim(objCOTAN@genesCoex),
                "cellsCoex"    = dim(objCOTAN@cellsCoex),
                "metaDataset"  = nrow(objCOTAN@metaDataset),
                "metaGenes"    = ncol(objCOTAN@metaGenes),
                "metaCells"    = ncol(objCOTAN@metaCells),
                "clustersCoex" = length(objCOTAN@clustersCoex)))
  }
)


# ------- `COTAN` estimated data accessors ------

#' @aliases getNu
#'
#' @details `getNu()` extracts the `nu` array (normalized cells' counts
#'   averages)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNu()` returns the `nu` array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "getNu",
  "COTAN",
  function(objCOTAN) {
    nu <- getMetadataCells(objCOTAN)[["nu"]]

    if (is_empty(nu)) {
      warning("`nu` is empty")
    } else {
      names(nu) <- getCells(objCOTAN)
    }

    return(nu)
  }
)


#' @aliases getLambda
#'
#' @details `getLambda()` extracts the `lambda` array (mean expression for each
#'   gene)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getLambda()` returns the `lambda` array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "getLambda",
  "COTAN",
  function(objCOTAN) {
    lambda <- getMetadataGenes(objCOTAN)[["lambda"]]

    if (is_empty(lambda)) {
      warning("`lambda` is empty")
    } else {
      names(lambda) <- getGenes(objCOTAN)
    }

    return(lambda)
  }
)


#' @aliases getDispersion
#'
#' @details `getDispersion()` extracts the `dispersion` array (one value for
#'   each gene)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getDispersion()` returns the `dispersion` array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "getDispersion",
  "COTAN",
  function(objCOTAN) {
    dispersion <- getMetadataGenes(objCOTAN)[["dispersion"]]

    if (is_empty(dispersion)) {
      warning("`dispersion` is empty")
    } else {
      names(dispersion) <- getGenes(objCOTAN)
    }

    return(dispersion)
  }
)


#' @details `estimatorsAreReady()` checks whether the estimators arrays
#'   `lambda`, `nu`, `dispersion` are available
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `estimatorsAreReady()` returns a boolean specifying whether all
#'   three arrays are non-empty
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname ParametersEstimations
#'
estimatorsAreReady <- function(objCOTAN) {
  anyEmptyArrays <- is_empty(getLambda(objCOTAN)) ||
                    is_empty(getNu(objCOTAN)) ||
                    is_empty(getDispersion(objCOTAN))
  if (anyEmptyArrays) {
    logThis(paste0("Estimators are not ready - array sizes: `lambda` ",
                   length(getLambda(objCOTAN)), ", `nu` ",
                   length(getNu(objCOTAN)), ", `dispersion` ",
                   length(getDispersion(objCOTAN))), logLevel = 2L)
  }
  return(!anyEmptyArrays)
}


#' @aliases getMu
#' @aliases calculateMu
#'
#' @details `getMu()` calculates the vector \eqn{\mu = \lambda \times
#'   \nu^T}
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getMu()` returns the `mu` matrix
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
getMu <- function(objCOTAN) {
  lambda <- suppressWarnings(getLambda(objCOTAN))
  assert_that(!is_empty(lambda),
              msg = "`lambda` must not be empty, estimate it")

  nu <- suppressWarnings(getNu(objCOTAN))
  assert_that(!is_empty(nu),
              msg = "`nu` must not be empty, estimate it")

  return(lambda %o% nu)
}


#' @details `getNuNormData()` extracts the \eqn{\nu}*-normalized* count table
#'   (i.e. where each column is divided by `nu`) and returns it
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNuNormData()` returns the \eqn{\nu}*-normalized* count
#'   `data.frame`
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' nuNorm <- getNuNormData(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
getNuNormData <- function(objCOTAN) {
  nu <- suppressWarnings(getNu(objCOTAN))
  assert_that(!is_empty(nu),
              msg = "`nu` must not be empty, estimate it")

  return(t(t(getRawData(objCOTAN)) * (1.0 / nu)))
}

#' @details `getLogNormData()` extracts the *log-normalized* count table (i.e.
#'   where each column is divided by the [getCellsSize()]), takes its `log10`
#'   and returns it.
#'
#' @param objCOTAN a `COTAN` object
#' @param retLog When `TRUE` returns
#'
#' @returns `getLogNormData()` returns a `data.frame` after applying the formula
#'   \eqn{\log_{10}{(10^4 * x + 1)}} to the raw counts normalized by
#'   *cells-size*
#'
#' @export
#'
#' @examples
#' logNorm <- getLogNormData(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
getLogNormData <- function(objCOTAN) {
  normData <- t(t(getRawData(objCOTAN)) * (1.0e4 / getCellsSize(objCOTAN)))

  return(log1p(normData) / log(10.0))
}

#' @details `getNormalizedData()` is deprecated: please use [getNuNormData()] or
#'   [getLogNormData()] directly as appropriate
#'
#' @param objCOTAN a `COTAN` object
#' @param retLog When `TRUE` calls [getLogNormData()], calls [getNuNormData()]
#'
#' @returns `getNormalizedData()` returns a `data.frame`
#'
#' @export
#'
#' @examples
#' logNorm <- getNormalizedData(objCOTAN, retLog = TRUE)
#'
#' @rdname ParametersEstimations
#'
getNormalizedData <- function(objCOTAN, retLog = FALSE) {
  if (retLog) {
    return(getLogNormData(objCOTAN))
  } else {
    return(getNuNormData(objCOTAN))
  }
}


#' @details `getProbabilityOfZero()` gives for each cell and each gene the
#'   probability of observing zero reads
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getProbabilityOfZero()` returns a `data.frame` with the
#'   probabilities of zero
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' probZero <- getProbabilityOfZero(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
getProbabilityOfZero <- function(objCOTAN) {
  dispersion <- suppressWarnings(getDispersion(objCOTAN))
  assert_that(!is_empty(dispersion),
              msg = "`dispersion` must not be empty, estimate it")

  # estimate Probabilities of 0 with internal function funProbZero
  return(funProbZero(dispersion, getMu(objCOTAN)))
}


# ----------- Raw data cleaning ------------

#' @title Raw data cleaning
#'
#' @description These methods are to be used to clean the raw data. That is drop
#'   any number of genes/cells that are too sparse or too present to allow
#'   proper calibration of the `COTAN` model.
#'
#'   We call genes that are expressed in all cells *Fully-Expressed* while cells
#'   that express all genes in the data are called *Fully-Expressing*. In case
#'   it has been made quite easy to exclude the flagged genes/cells in the user
#'   calculations.
#'
#' @name RawDataCleaning
NULL

#' @aliases flagNotFullyExpressedGenes
#'
#' @details `flagNotFullyExpressedGenes()` returns a Boolean array with TRUE for
#'   those genes that are not fully-expressed.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `flagNotFullyExpressedGenes()` returns a Booleans array with TRUE
#'   for genes that are not fully-expressed
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @export
#'
#' @examples
#' library(zeallot)
#'
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#'
#' genes.to.rem <- getGenes(objCOTAN)[grep('^MT', getGenes(objCOTAN))]
#' cells.to.rem <- getCells(objCOTAN)[which(getCellsSize(objCOTAN) == 0)]
#' objCOTAN <- dropGenesCells(objCOTAN, genes.to.rem, cells.to.rem)
#'
#' objCOTAN <- clean(objCOTAN)
#'
#' objCOTAN <- findFullyExpressedGenes(objCOTAN)
#' goodPos <- flagNotFullyExpressedGenes(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "flagNotFullyExpressedGenes",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getMetadataGenes(objCOTAN)[["feGenes"]])) {
      return(set_names(rep(TRUE, getNumGenes(objCOTAN)), getGenes(objCOTAN)))
    } else {
      return(!getMetadataGenes(objCOTAN)[["feGenes"]])
    }
  }
)


#' @aliases flagNotFullyExpressingCells
#'
#' @details `flagNotFullyExpressingCells()`returns a Boolean vector with TRUE
#'   for those cells that are not expressing all genes
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `flagNotFullyExpressingCells()` returns an array of Booleans with
#'   TRUE for cells that are not expressing all genes
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @export
#'
#' @examples
#' objCOTAN <- findFullyExpressingCells(objCOTAN)
#' goodPos <- flagNotFullyExpressingCells(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "flagNotFullyExpressingCells",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getMetadataCells(objCOTAN)[["feCells"]])) {
      return(set_names(rep(TRUE, getNumCells(objCOTAN)), getCells(objCOTAN)))
    } else {
      return(!getMetadataCells(objCOTAN)[["feCells"]])
    }
  }
)


#' @aliases getFullyExpressedGenes
#'
#' @details `getFullyExpressedGenes()` returns the genes expressed in all cells
#'   of the dataset
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getFullyExpressedGenes()` returns an array containing all genes
#'   that are expressed in all cells
#'
#' @export
#'
#' @examples
#' feGenes <- getFullyExpressedGenes(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "getFullyExpressedGenes",
  "COTAN",
  function(objCOTAN) {
    return(getGenes(objCOTAN)[!flagNotFullyExpressedGenes(objCOTAN)])
  }
)


#' @aliases getFullyExpressingCells
#'
#' @details `getFullyExpressingCells()` returns the cells that did express
#'   all genes of the dataset
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getFullyExpressingCells()` returns an array containing all cells
#'   that express all genes
#'
#' @export
#'
#' @examples
#' feCells <- getFullyExpressingCells(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "getFullyExpressingCells",
  "COTAN",
  function(objCOTAN) {
    return(getCells(objCOTAN)[!flagNotFullyExpressingCells(objCOTAN)])
  }
)

# ------- `COTAN` coex data accessors ------

#' Calculating the COEX matrix for genes and cells
#'
#' @name CalculatingCOEX
NULL

#' @aliases getGenesCoex
#'
#' @details `getGenesCoex()` extracts a complete (or a partial after genes
#'   dropping) genes' `COEX` matrix from the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes The given genes' names to select the wanted `COEX` columns. If
#'   missing all columns will be returned. When not empty a proper result is
#'   provided by calculating the partial `COEX` matrix on the fly
#' @param zeroDiagonal When TRUE the `COEX` of any element with itself is set to
#'   zero
#' @param ignoreSync When `TRUE` ignores whether the `lambda`/`nu`/`dispersion`
#'   have been updated since the `COEX` matrix was calculated.
#'
#' @returns `getGenesCoex()` returns the genes' `COEX` values
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- initializeMetaDataset(objCOTAN, GEO = "test_GEO",
#'                                   sequencingMethod = "distribution_sampling",
#'                                   sampleCondition = "reconstructed_dataset")
#' objCOTAN <- clean(objCOTAN)
#'
#' objCOTAN <- estimateLambdaLinear(objCOTAN)
#' objCOTAN <- estimateDispersionViaSolver(objCOTAN, cores = 6L)
#'
#' ## Now the `COTAN` object is ready to calculate the genes' `COEX`
#'
#' ## mu <- getMu(objCOTAN)
#' ## observedY <- observedContingencyTablesYY(objCOTAN, asDspMatrices = TRUE)
#' obs <- observedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' ## expectedN <- expectedContingencyTablesNN(objCOTAN, asDspMatrices = TRUE)
#' exp <- expectedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#'
#' stopifnot(isCoexAvailable(objCOTAN))
#' genesCoex <- getGenesCoex(objCOTAN)
#' genesSample <- sample(getNumGenes(objCOTAN), 10)
#' partialGenesCoex <- calculatePartialCoex(objCOTAN, genesSample,
#'                                          actOnCells = FALSE)
#'
#' stopifnot(all(1e-6 >
#'                 abs(partialGenesCoex -
#'                       getGenesCoex(objCOTAN,
#'                                    getGenes(objCOTAN)[sort(genesSample)],
#'                                    zeroDiagonal = FALSE))))
#'
#' ## S <- calculateS(objCOTAN)
#' ## G <- calculateG(objCOTAN)
#' ## pValue <- calculatePValue(objCOTAN)
#' gdiDF <- calculateGDI(objCOTAN)
#' objCOTAN <- storeGDI(objCOTAN, genesGDI = gdiDF)
#'
#' ## Touching any of the `lambda`/`nu`/`dispersion` parameters invalidates the
#' ## `COEX` matrix and derivatives, so it can be dropped it from the `COTAN`
#' ## object
#' objCOTAN <- dropGenesCoex(objCOTAN)
#' stopifnot(!isCoexAvailable(objCOTAN))
#'
#'
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 6L)
#'
#' ## Now the `COTAN` object is ready to calculate the cells' `COEX`
#' ## In case one needs to calculate both, it is more sensible to run the above
#' ## before any `COEX` evaluation
#'
#' g1 <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 1)]
#' g2 <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 1)]
#' tables <- contingencyTables(objCOTAN, g1 = g1, g2 = g2)
#' tables
#'
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#' stopifnot(isCoexAvailable(objCOTAN, actOnCells = TRUE, ignoreSync = TRUE))
#' cellsCoex <- getCellsCoex(objCOTAN, zeroDiagonal = FALSE)
#'
#' cellsSample <- sample(getNumCells(objCOTAN), 10)
#' partialCellsCoex <- calculatePartialCoex(objCOTAN, cellsSample,
#'                                          actOnCells = TRUE)
#'
#' stopifnot(all(1e-6 >
#'                 abs(partialCellsCoex - cellsCoex[, sort(cellsSample)])))
#'
#' objCOTAN <- dropCellsCoex(objCOTAN)
#' stopifnot(!isCoexAvailable(objCOTAN, actOnCells = TRUE))
#'
#' @rdname CalculatingCOEX
#'
setMethod(
  "getGenesCoex",
  "COTAN",
  function(objCOTAN, genes = vector(mode = "character"),
           zeroDiagonal = TRUE, ignoreSync = FALSE) {
    coexIsReady <- isCoexAvailable(objCOTAN, ignoreSync = ignoreSync)

    if (is_empty(genes)) {
      assert_that(coexIsReady,
                  msg = paste0("Cannot return genes' COEX as the matrix was ",
                               "never calculated or is out-of-sync with the ",
                               "estimators."))

      ret <- objCOTAN@genesCoex
      if (isTRUE(zeroDiagonal)) {
        diag(ret) <- 0.0
        ret <- pack(forceSymmetric(ret))
      }
      return(ret)
    } else {
      genes <- handleNamesSubsets(getGenes(objCOTAN), genes)

      ret <- objCOTAN@genesCoex

      if (!coexIsReady) {
        warning("Missing or out-of-sync genes' COEX:",
                "calculating the required subset now!")
        ret <- calculatePartialCoex(objCOTAN, columnsSubset = genes,
                                    actOnCells = FALSE)
      } else {
        ret <- ret[, genes, drop = FALSE]
      }

      if (isTRUE(zeroDiagonal)) {
        diag(ret[genes, , drop = FALSE]) <- 0.0
      }

      return(ret)
    }
  }
)


#' @aliases getCellsCoex
#'
#' @details `getCellsCoex()` extracts a complete (or a partial after cells
#'   dropping) cells' `COEX` matrix from the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param cells The given cells' names to select the wanted `COEX` columns. If
#'   missing all columns will be returned. When not empty a proper result is
#'   provided by calculating the partial `COEX` matrix on the fly
#' @param zeroDiagonal When `TRUE` sets the diagonal to zero.
#' @param ignoreSync When `TRUE` ignores whether the `lambda`/`nu`/`dispersion`
#'   have been updated since the `COEX` matrix was calculated.
#'
#' @returns `getCellsCoex()` returns the cells' `COEX` values
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
setMethod(
  "getCellsCoex",
  "COTAN",
  function(objCOTAN, cells = vector(mode = "character"),
           zeroDiagonal = TRUE, ignoreSync = FALSE) {
    coexIsReady <- isCoexAvailable(objCOTAN, actOnCells = TRUE,
                                   ignoreSync = ignoreSync)

    if (is_empty(cells)) {
      assert_that(coexIsReady,
                  msg = paste0("Cannot return cells' COEX as the matrix was ",
                               "never calculated or is out-of-sync with the ",
                               "estimators."))

      ret <- objCOTAN@cellsCoex
      if (isTRUE(zeroDiagonal)) {
        diag(ret) <- 0.0
        ret <- pack(forceSymmetric(ret))
      }
      return(ret)
    } else {
      cells <- handleNamesSubsets(getCells(objCOTAN), cells)

      ret <- objCOTAN@cellsCoex

      if (!coexIsReady) {
        warning("Missing or out-of-sync cells' COEX:",
                "calculating the required subset now!")
        ret <- calculatePartialCoex(objCOTAN, columnsSubset = cells,
                                    actOnCells = TRUE)
      } else {
        ret <- ret[, cells, drop = FALSE]
      }

      if (isTRUE(zeroDiagonal)) {
        diag(ret[cells, , drop = FALSE]) <- 0.0
      }

      return(ret)
    }
  }
)


#' @aliases isCoexAvailable
#'
#' @details `isCoexAvailable()` allows to query whether the relevant `COEX`
#'   matrix from the `COTAN` object is available to use
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param ignoreSync When `TRUE` ignores whether the `lambda`/`nu`/`dispersion`
#'   have been updated since the `COEX` matrix was calculated.
#'
#' @returns `isCoexAvailable()` returns whether relevant `COEX` matrix has been
#'   calculated and, in case, if it is still aligned to the estimators.
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
setMethod(
  "isCoexAvailable",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE, ignoreSync = FALSE) {
    if (isTRUE(actOnCells)) {
      isInSync <-
        as.logical(getMetadataElement(objCOTAN, datasetTags()[["csync"]]))

      return(!is_empty(objCOTAN@cellsCoex) &&
               (isTRUE(ignoreSync) || isTRUE(isInSync)))
    } else {
      isInSync <-
        as.logical(getMetadataElement(objCOTAN, datasetTags()[["gsync"]]))

      return(!is_empty(objCOTAN@genesCoex) &&
               (isTRUE(ignoreSync) || isTRUE(isInSync)))
    }
  }
)


#' @aliases getGDI
#'
#' @details `getGDI()` extracts the genes' `GDI` array as it was stored by the
#'   method [storeGDI()]
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getGDI()` returns the genes' `GDI`` array if available or `NULL`
#'   otherwise
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname GenesStatistics
#'
setMethod(
  "getGDI",
  "COTAN",
  function(objCOTAN) {
    gdi <- getMetadataGenes(objCOTAN)[["GDI"]]

    if (is_empty(gdi)) {
      warning("no GDI data has been stored")
    } else {
      names(gdi) <- getGenes(objCOTAN)
    }

    return(gdi)
  }
)


#' @aliases getClusterizations
#'
#' @details `getClusterizations()` extracts the list of the *clusterizations*
#'   defined in the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param dropNoCoex When `TRUE` drops the names from the *clusterizations* with
#'   empty associated `COEX` `data.frame`
#' @param keepPrefix When `TRUE` returns the internal name of the
#'   *clusterization*: the one with the `CL_` prefix.
#'
#' @returns `getClusterizations()` returns a vector of *clusterization* names,
#'   usually without the `CL_` prefix
#'
#' @importFrom rlang is_empty
#'
#' @importFrom methods validObject
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 6L, calcCoex = TRUE,
#'                           optimizeForSpeed = TRUE, saveObj = FALSE)
#'
#' data("test.dataset.clusters1")
#' clusters <- test.dataset.clusters1
#'
#' coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)
#'
#' groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030",
#'                             "g-000150", "g-000160", "g-000170"),
#'                      G2 = c("g-000300", "g-000330", "g-000450",
#'                             "g-000460", "g-000470"),
#'                      G3 = c("g-000510", "g-000530", "g-000550",
#'                             "g-000570", "g-000590"))
#'
#' geneClusters <- rep(1:3, each = 240)[1:600]
#' names(geneClusters) <- getGenes(objCOTAN)
#'
#' umapPlot <- UMAPPlot(coexDF, clusters = NULL, elements = groupMarkers)
#' plot(umapPlot)
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "first_clusterization",
#'                               clusters = clusters, coexDF = coexDF)
#'
#' lfcDF <- logFoldChangeOnClusters(objCOTAN, clusters = clusters)
#' umapPlot2 <- UMAPPlot(lfcDF, clusters = geneClusters)
#' plot(umapPlot2)
#'
#' if (FALSE) {
#'   objCOTAN <- estimateNuLinearByCluster(objCOTAN, clusters = clusters)
#' }
#'
#' clSummaryPlotAndData <-
#'   clustersSummaryPlot(objCOTAN, clName = "first_clusterization",
#'                       plotTitle = "first clusterization")
#' plot(clSummaryPlotAndData[["plot"]])
#'
#' if (FALSE) {
#'   objCOTAN <- dropClusterization(objCOTAN, "first_clusterization")
#' }
#'
#' clusterizations <- getClusterizations(objCOTAN, dropNoCoex = TRUE)
#' stopifnot(length(clusterizations) == 1)
#'
#' cellsUmapPlotAndDF <- cellsUMAPPlot(objCOTAN, dataMethod = "LogLikelihood",
#'                                     useCoexEigen = TRUE, numComp = 25L,
#'                                     clName = "first_clusterization")
#' plot(cellsUmapPlotAndDF[["plot"]])
#'
#' enrichment <- geneSetEnrichment(clustersCoex = coexDF,
#'                                 groupMarkers = groupMarkers)
#'
#' clHeatmapPlotAndData <- clustersMarkersHeatmapPlot(objCOTAN, groupMarkers)
#' clHeatmapPlotAndData[["heatmapPlot"]]
#'
#' conditions <- as.integer(substring(getCells(objCOTAN), 3L))
#' conditions <- factor(ifelse(conditions <= 600, "L", "H"))
#' names(conditions) <- getCells(objCOTAN)
#'
#' objCOTAN <- addCondition(objCOTAN, condName = "High/Low",
#'                          conditions = conditions)
#'
#' clHeatmapPlotAndData2 <-
#'   clustersMarkersHeatmapPlot(objCOTAN, groupMarkers, kCuts = 2,
#'                              condNameList = list("High/Low"))
#' clHeatmapPlotAndData2[["heatmapPlot"]]
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "getClusterizations",
  "COTAN",
  function(objCOTAN, dropNoCoex = FALSE, keepPrefix = FALSE) {
    validObject(objCOTAN)

    clsCoex <- getClustersCoex(objCOTAN)

    if (isTRUE(dropNoCoex)) {
      emptyClsCoex <- vapply(clsCoex, is_empty, logical(1L))
      out <- names(clsCoex[!emptyClsCoex])
    } else {
      out <- names(clsCoex)
    }

    # drop the internal 'CL_' prefix
    if (isFALSE(keepPrefix)) {
      out <- substring(out, 4L)
    }

    return(out)
  }
)

#' @aliases getClusterizationName
#'
#' @details `getClusterizationName()` normalizes the given *clusterization* name
#'   or, if none were given, returns the name of last available *clusterization*
#'   in the `COTAN` object. It can return the *clusterization* **internal name**
#'   if needed
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available
#'   *clusterization* will be returned, as it is probably the most significant!
#' @param keepPrefix When `TRUE` returns the internal name of the
#'   *clusterization*: the one with the `CL_` prefix.
#'
#' @returns `getClusterizationName()` returns the normalized *clusterization*
#'   name or `NULL` if no *clusterizations* are present
#'
#' @export
#'
#' @examples
#' clName <- getClusterizationName(objCOTAN)
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "getClusterizationName",
  "COTAN",
  function(objCOTAN, clName = "", keepPrefix = FALSE) {
    allClustNames <- getClusterizations(objCOTAN, keepPrefix = TRUE)

    if (isEmptyName(clName)) {
      # pick last clusterization
      outName <- allClustNames[length(allClustNames)]
    } else {
      outName <- clName
      if (!startsWith(clName, "CL_")) {
        outName <- paste0("CL_", clName)
      }
      assert_that(outName %in% allClustNames,
                  msg = paste0("Given cluster name '", clName,
                               "' is not among the stored clusterizations"))
    }

    # drop the internal 'CL_' prefix
    if (isFALSE(keepPrefix)) {
      outName <- substring(outName, 4L)
    }

    return(outName)
  }
)

#' @aliases getClusterizationData
#'
#' @details `getClusterizationData()` extracts the asked *clusterization* and
#'   its associated `COEX` `data.frame` from the `COTAN` object
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be returned, as it is probably the most
#'   significant!
#'
#' @returns `getClusterizationData()` returns a `list` with 2 elements:
#'   * `"clusters"` the named cluster labels array
#'   * `"coex"` the associated `COEX` `data.frame`. This will be an **empty**
#'     `data.frame` when not specified for the relevant *clusterization*
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @examples
#' clusterDataList <- getClusterizationData(objCOTAN, clName = clName)
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "getClusterizationData",
  "COTAN",
  function(objCOTAN, clName = "") {
    internalName <- getClusterizationName(objCOTAN, clName = clName,
                                          keepPrefix = TRUE)

    # clName can still be empty if no clusterization was store in the objCOTAN
    assert_that(!isEmptyName(internalName),
                msg = "No clusterizations are present in the 'COTAN' object")

    clusters <- factor(set_names(getMetadataCells(objCOTAN)[[internalName]],
                                 getCells(objCOTAN)))

    coexDF <- getClustersCoex(objCOTAN)[[internalName]]

    return(list("clusters" = clusters, "coex" = coexDF))
  }
)

#'
#' @details `getClusters()` extracts the asked *clusterization* from the `COTAN`
#'   object
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be returned, as it is probably the most
#'   significant!
#'
#' @returns `getClusters()` returns the named cluster labels array
#'
#' @export
#'
#' @examples
#' clusters <- getClusters(objCOTAN, clName = clName)
#'
#' @rdname HandlingClusterizations
#'
getClusters <- function(objCOTAN, clName = "") {
  return(getClusterizationData(objCOTAN, clName = clName)[["clusters"]])
}

#' @aliases getClustersCoex
#'
#' @details `getClustersCoex()` extracts the full `clusterCoex` member `list`
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getClustersCoex()` returns the list with a `COEX` `data.frame` for
#'   each *clusterization*. When not empty, each `data.frame` contains a `COEX`
#'   column for each *cluster*.
#'
#' @export
#'
#' @examples
#' allClustersCoexDF <- getClustersCoex(objCOTAN)
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "getClustersCoex",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@clustersCoex)
  }
)


# ------- `COTAN` conditions data accessors ------

#' @title Handling cells' *conditions* and related functions
#'
#' @description These functions manage the *conditions*.
#'
#'   A *condition* is a set of **labels** that can be assigned to cells:
#'   one **label** per cell. This is especially useful in cases when the
#'   `data-set` is the result of merging multiple experiments' raw data
#'
#' @name HandlingConditions
NULL

#' @aliases getAllConditions
#'
#' @details `getAllConditions()` extracts the list of the *conditions* defined
#'   in the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param keepPrefix When `TRUE` returns the internal name of the
#'   *condition*: the one with the `COND_` prefix.
#'
#' @returns `getAllConditions()` returns a vector of *conditions* names,
#'   usually without the `COND_` prefix
#'
#' @importFrom rlang is_empty
#'
#' @importFrom methods validObject
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#'
#' cellLine <- rep(c("A", "B"), getNumCells(objCOTAN) / 2)
#' names(cellLine) <- getCells(objCOTAN)
#' objCOTAN <- addCondition(objCOTAN, condName = "Line", conditions = cellLine)
#'
#' if (FALSE) {
#'   objCOTAN <- dropCondition(objCOTAN, "Genre")
#' }
#'
#' conditionsNames <- getAllConditions(objCOTAN)
#'
#' @rdname HandlingConditions
#'
setMethod(
  "getAllConditions",
  "COTAN",
  function(objCOTAN, keepPrefix = FALSE) {
    validObject(objCOTAN)

    cNames <- colnames(getMetadataCells(objCOTAN))
    areCondNames <- vapply(cNames, startsWith, logical(1L), "COND_")
    out <- cNames[areCondNames]

    # drop the internal 'COND_' prefix
    if (isFALSE(keepPrefix)) {
      out <- substring(out, 6L)
    }

    return(out)
  }
)

#' @aliases getConditionName
#'
#' @details `getConditionName()` normalizes the given *condition* name or, if
#'   none were given, returns the name of last available *condition* in the
#'   `COTAN` object. It can return the *condition* **internal name** if needed
#'
#' @param objCOTAN a `COTAN` object
#' @param condName The name of the *condition*. If not given an empty string
#'   will be returned!
#' @param keepPrefix When `TRUE` returns the internal name of the
#'   *condition*: the one with the `COND_` prefix.
#'
#' @returns `getConditionName()` returns the normalized *condition* name or
#'   `NULL` if no *conditions* are present
#'
#' @export
#'
#' @examples
#' condName <- getConditionName(objCOTAN)
#'
#' @rdname HandlingConditions
#'
setMethod(
  "getConditionName",
  "COTAN",
  function(objCOTAN, condName = "", keepPrefix = FALSE) {
    allCondNames <- getAllConditions(objCOTAN, keepPrefix = TRUE)

    if (isEmptyName(condName)) {
      # return a default name
      outName <- "COND_"
    } else {
      outName <- condName
      if (!startsWith(condName, "COND_")) {
        outName <- paste0("COND_", condName)
      }
      assert_that(outName %in% allCondNames,
                  msg = paste0("Given condition name '", condName,
                               "' is not among the stored conditions"))
    }

    # drop the internal 'COND_' prefix
    if (isFALSE(keepPrefix)) {
      outName <- substring(outName, 6L)
    }

    return(outName)
  }
)

#' @aliases getCondition
#'
#' @details `getCondition()` extracts the asked *condition* from the `COTAN`
#'   object
#'
#' @param objCOTAN a `COTAN` object
#' @param condName The name of the *condition*. If not given a dummy constant
#'   **label** *condition* will be returned.
#'
#' @returns `getCondition()` returns a named `factor` with the *condition*
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @examples
#' condition <- getCondition(objCOTAN, condName = condName)
#' isa(condition, "factor")
#'
#' @rdname HandlingConditions
#'
setMethod(
  "getCondition",
  "COTAN",
  function(objCOTAN, condName = "") {
    internalName <- getConditionName(objCOTAN, condName = condName,
                                     keepPrefix = TRUE)

    if (internalName == "COND_") {
      # no corresponding condition
      conditions <- factor(rep_len("NoCond", getNumCells(objCOTAN)))
    } else {
      conditions <- getMetadataCells(objCOTAN)[[internalName]]
    }

    return(set_names(conditions, getCells(objCOTAN)))
  }
)


#' @details `normalizeNameAndLabels()` takes a pair of name/labels and
#'   normalize them based on the available information in the `COTAN` object
#'
#' @param objCOTAN a `COTAN` object
#' @param name the name of the *clusterization*/*condition*. If not given the
#'   last available *clusterization* will be used, or no *conditions*
#' @param labels a *clusterization*/*condition* to use. If given it will take
#'   precedence on the one indicated by `name`
#' @param isCond a Boolean to indicate whether the function is dealing with
#'   *clusterizations* [`FALSE`] or *conditions* [`TRUE`]
#'
#' @returns `normalizeNameAndLabels()` returns a `list` with:
#'   * `"name"` the relevant name
#'   * `"labels"` the relevant *clusterization*/*condition*
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#' nameAndCond <- normalizeNameAndLabels(objCOTAN, name = condName,
#'                                       isCond = TRUE)
#' isa(nameAndCond[["labels"]], "factor")
#'
#' @rdname HandlingConditions
#'
normalizeNameAndLabels <- function(objCOTAN, name = "",
                                   labels = NULL, isCond = FALSE) {
  if (is_empty(labels)) {
    if (isFALSE(isCond)) {
      name <- getClusterizationName(objCOTAN, clName = name)
      labels <- getClusters(objCOTAN, clName = name)
    } else {
      name <- getConditionName(objCOTAN, condName = name)
      labels <- getCondition(objCOTAN, condName = name)
    }
  } else {
    assert_that(!is_empty(names(labels)),
                msg = "No names attached to the given labels")
    assert_that(setequal(names(labels), getCells(objCOTAN)),
                msg = "Non compatible labels")

    if (isEmptyName(name)) {
      name <- "DummyName"
    }

    if (!inherits(labels, "factor")) {
      labels <- factor(labels)
    }
  }
  return(list("name" = name, "labels" = labels))
}
