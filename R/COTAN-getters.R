# ------- `COTAN` raw data accessors --------

#' Raw data `COTAN` accessors
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

#' Handling *meta-data* in `COTAN` objects
#'
#' @description Much of the information stored in the `COTAN` object is
#'   compacted into three `data.frame`s:
#'   * "metaDataset" - contains all general information about the data-set
#'   * "metaGenes" - contains genes' related information along the `lambda` and
#'     `dispersion` vectors and the housekeeping flag
#'   * "metaCells" - contains cells' related information along the `nu` vector,
#'     the fully expressed flag and the clusterizations
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
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' numInitialCells <- getMetadataElement(objCOTAN, "cells")
#'
#' @rdname HandleMetaData
#'
setMethod(
  "getMetadataElement",
  "COTAN",
  function(objCOTAN, tag) {
    meta <- getMetadataDataset(objCOTAN)

    if (is_empty(meta) || !(tag %in% meta[[1L]])) {
      out <- ""
    } else {
      rowPos <- which(meta[[1L]] %in% tag)
      out <- meta[, -1L, drop = FALSE][rowPos, ]
    }

    return(out)
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

#' @aliases getNormalizedData
#'
#' @details `getNormalizedData()` extracts the *normalized* count table (i.e.
#'   divided by `nu`)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNormalizedData()` returns the normalized count `data.frame`
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' rawNorm <- getNormalizedData(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "getNormalizedData",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getNu(objCOTAN))) {
      stop("nu must not be empty, estimate it")
    }

    return(t(t(getRawData(objCOTAN)) * (1.0 / getNu(objCOTAN))))
  }
)


#' @aliases getNu
#'
#' @details `getNu()` extracts the nu array (normalized cells' counts averages)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNu()` returns the nu array
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
      warning("nu is empty")
    } else {
      names(nu) <- getCells(objCOTAN)
    }

    return(nu)
  }
)


#' @aliases getLambda
#'
#' @details `getLambda()` extracts the lambda array (mean expression for each
#'   gene)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getLambda()` returns the lambda array
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
      warning("lambda is empty")
    } else {
      names(lambda) <- getGenes(objCOTAN)
    }

    return(lambda)
  }
)


#' @aliases getDispersion
#'
#' @details `getDispersion()` extracts the dispersion array (a)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getDispersion()` returns the dispersion array
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
      warning("dispersion is empty")
    } else {
      names(dispersion) <- getGenes(objCOTAN)
    }

    return(dispersion)
  }
)


# ----------- Raw data cleaning ------------

#' Raw data cleaning
#'
#' @description These methods are to be used to clean the raw data. That is drop
#'   any number of genes/cells that are too sparse or too present to allow
#'   proper calibration of the `COTAN` model.
#'
#'   We call genes that are expressed in all cells *House Keeping* while cells
#'   that express all genes in the data are called *Fully Expressed*. In case it
#'   has been made quite easy to exclude the flagged genes/cells in the user
#'   calculations.
#'
#' @name RawDataCleaning
NULL

#' @aliases flagNotHousekeepingGenes
#'
#' @details `flagNotHousekeepingGenes()` returns a Boolean array with TRUE for
#'   those genes that are not housekeeping.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `flagNotHousekeepingGenes()` returns a Booleans array with TRUE for
#'   genes that are not housekeeping
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#'
#' genes.to.rem <- getGenes(objCOTAN)[grep('^MT', getGenes(objCOTAN))]
#' cells.to.rem <- getCells(objCOTAN)[which(getCellsSize(objCOTAN) == 0)]
#' objCOTAN <- dropGenesCells(objCOTAN, genes.to.rem, cells.to.rem)
#'
#' objCOTAN <- clean(objCOTAN)
#'
#' objCOTAN <- findHousekeepingGenes(objCOTAN)
#' goodPos <- flagNotHousekeepingGenes(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "flagNotHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getMetadataGenes(objCOTAN)[["hkGenes"]])) {
      return(set_names(rep(TRUE, getNumGenes(objCOTAN)), getGenes(objCOTAN)))
    } else {
      return(!getMetadataGenes(objCOTAN)[["hkGenes"]])
    }
  }
)


#' @aliases flagNotFullyExpressedCells
#'
#' @details `flagNotFullyExpressedCells()`returns a Boolean vector with TRUE for
#'   those cells that are not fully expressed
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `flagNotFullyExpressedCells()` returns an array of Booleans with
#'   TRUE for cells that are not fully expressed
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @export
#'
#' @examples
#' objCOTAN <- findFullyExpressedCells(objCOTAN)
#' goodPos <- flagNotFullyExpressedCells(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "flagNotFullyExpressedCells",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getMetadataCells(objCOTAN)[["feCells"]])) {
      return(set_names(rep(TRUE, getNumCells(objCOTAN)), getCells(objCOTAN)))
    } else {
      return(!getMetadataCells(objCOTAN)[["feCells"]])
    }
  }
)


#' @aliases getHousekeepingGenes
#'
#' @details `getHousekeepingGenes()` returns the genes expressed in all cells of
#'   the dataset
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getHousekeepingGenes()` returns an array containing all genes
#'   that are expressed in all cells
#'
#' @export
#'
#' @examples
#' hkGenes <- getHousekeepingGenes(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "getHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    return(getGenes(objCOTAN)[!flagNotHousekeepingGenes(objCOTAN)])
  }
)


#' @aliases getFullyExpressedCells
#'
#' @details `getFullyExpressedCells()` returns the cells that did express
#'   all genes of the dataset
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getFullyExpressedCells()` returns an array containing all cells
#'   that express all genes
#'
#' @export
#'
#' @examples
#' feCells <- getFullyExpressedCells(objCOTAN)
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "getFullyExpressedCells",
  "COTAN",
  function(objCOTAN) {
    return(getCells(objCOTAN)[!flagNotFullyExpressedCells(objCOTAN)])
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
#' @param genes A vector of gene names. It will exclude any gene not on the
#'   list. By defaults the function will keep all genes.
#' @param zeroDiagonal When TRUE sets the diagonal to zero.
#' @param ignoreSync When `TRUE` ignores whether the `lambda`/`nu`/`dispersion`
#'   have been updated since the `COEX` matrix was calculated.
#'
#' @returns `getGenesCoex()` returns the genes' `COEX` values
#'
#' @importFrom rlang is_empty
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
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#'
#' ## Now the `COTAN` object is ready to calculate the genes' `COEX`
#'
#' ## mu <- calculateMu(objCOTAN)
#' ## observedY <- observedContingencyTablesYY(objCOTAN, asDspMatrices = TRUE)
#' obs <- observedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' ## expectedN <- expectedContingencyTablesNN(objCOTAN, asDspMatrices = TRUE)
#' exp <- expectedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' genesCoex <- getGenesCoex(objCOTAN)
#'
#' ## S <- calculateS(objCOTAN)
#' ## G <- calculateG(objCOTAN)
#' ## pValue <- calculatePValue(objCOTAN)
#' GDI <- calculateGDI(objCOTAN)
#'
#' ## Touching any of the lambda/nu/dispersino parameters invalidates the `COEX`
#' ## matrix and derivatives, so it can be dropped it from the `COTAN` object
#' objCOTAN <- dropGenesCoex(objCOTAN)
#'
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#'
#' ## Now the `COTAN` object is ready to calculate the cells' `COEX`
#' ## In case one need to caclualte both it is more sensible to run the above
#' ## before any `COEX` evaluation
#'
#' g1 <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 1)]
#' g2 <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 1)]
#' tables <- contingencyTables(objCOTAN, g1 = g1, g2 = g2)
#' tables
#'
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#' cellsCoex <- getCellsCoex(objCOTAN)
#'
#' objCOTAN <- dropCellsCoex(objCOTAN)
#'
#' @rdname CalculatingCOEX
#'
setMethod(
  "getGenesCoex",
  "COTAN",
  function(objCOTAN, genes = c(), zeroDiagonal = TRUE, ignoreSync = FALSE) {
    if (is_empty(objCOTAN@genesCoex)) {
      stop("Cannot return genes' coex as the matrix was never calculated")
    }

    if (isFALSE(ignoreSync) &&
        isFALSE(as.logical(getMetadataElement(objCOTAN,
                                              datasetTags()[["gsync"]])))) {
      stop("Cannot return genes' coex as the matrix is",
           " out of sync with the other parameters.",
           " It is still ossible to access the data directly!")
    }

    ret <- objCOTAN@genesCoex

    if (isTRUE(zeroDiagonal)) {
      ret@x[cumsum(seq_len(nrow(ret)))] <- 0.0
    }

    if (!is_empty(genes)) {
      ret <- ret[, getGenes(objCOTAN) %in% genes, drop = FALSE]
    }

    return(ret)
  }
)


#' @aliases getCellsCoex
#'
#' @details `getCellsCoex()` extracts a complete (or a partial after cells
#'   dropping) cells' `COEX` matrix from the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param cells A vector of cell names. It will exclude any cell not on the
#'   list. By defaults the function will keep all cells.
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
  function(objCOTAN, cells = c(), zeroDiagonal = TRUE, ignoreSync = FALSE) {
    if (is_empty(objCOTAN@cellsCoex)) {
      stop("Cannot return cells' coex as the matrix was never calculated")
    }

    if (isFALSE(ignoreSync) &&
        isFALSE(as.logical(getMetadataElement(objCOTAN,
                                              datasetTags()[["csync"]])))) {
      stop("Cannot return cells' coex as the matrix is",
           " out of sync with the other parameters.",
           " It is still ossible to access the data directly!")
    }

    ret <- objCOTAN@cellsCoex

    if (isTRUE(zeroDiagonal)) {
      ret@x[cumsum(seq_len(nrow(ret)))] <- 0.0
    }

    if (!is_empty(cells)) {
      ret <- ret[, getCells(objCOTAN) %in% cells, drop = FALSE]
    }

    return(ret)
  }
)


# ------- `COTAN` clusterization data accessors ------

#' Handling cells' *clusterization* and related functions
#'
#' @description These functions manage the *clusterizations* and their
#'   associated *cluster* `COEX` `data.frame`s.
#'
#'   A *clusterization* is any partition of the cells where to each cell it is
#'   assigned a **label**; a group of cells with the same label is called
#'   *cluster*.
#'
#'   For each cluster is also possible to define a `COEX` value for each gene,
#'   indicating its increased or decreased expression in the *cluster* compared
#'   to the whole background. A `data.frame` with these values listed in a
#'   column for each *cluster* is stored separately for each *clusterization* in
#'   the `clustersCoex` member.
#'
#'   The formulae for this *In/Out* `COEX` are similar to those used in the
#'   [calculateCoex()] method, with the **role** of the second gene taken by the
#'   *In/Out* status of the cells with respect to each *cluster*.
#'
#' @name HandlingClusterizations
NULL

#' @aliases getClusterizations
#'
#' @details `getClusterizations()` extracts the list of the *clusterizations*
#'   defined in the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param dropNoCoex When `TRUE` drops the names from the clusterizations with
#'   empty associated coex `data.frame`
#' @param keepPrefix When `TRUE` returns the internal name of the
#'   clusterization: the one with the `CL_` prefix.
#'
#' @returns `getClusterizations()` returns a vector of *clusterization* names,
#'   usually without the `CL_` prefix
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#'
#' data("test.dataset.clusters1")
#' clusters <- test.dataset.clusters1
#'
#' coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)[["coex"]]
#'
#' groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
#'                      G2 = c("g-000300", "g-000330"),
#'                      G3 = c("g-000510", "g-000530", "g-000550",
#'                             "g-000570", "g-000590"))
#'
#' umapPlot <- UMAPPlot(coexDF, clusters = NULL, elements = groupMarkers)
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "first_clusterization",
#'                               clusters = clusters, coexDF = coexDF)
#'
#' ##objCOTAN <- dropClusterization(objCOTAN, "first_clusterization")
#'
#' clusterizations <- getClusterizations(objCOTAN, dropNoCoex = TRUE)
#'
#' enrichment <- geneSetEnrichment(clustersCoex = coexDF,
#'                                 groupMarkers = groupMarkers)
#'
#' ##clHeatmapPlot <- clustersMarkersHeatmapPlot(objCOTAN, groupMarkers)
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
      emptyClsCoex <- vapply(clsCoex, is_empty, logical(1))
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
#' @param clName The name of the clusterization. If not given the last available
#'   clusterization will be returned, as it is probably the most significant!
#' @param keepPrefix When `TRUE` returns the internal name of the
#'   clusterization: the one with the `CL_` prefix.
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

    emptyName <- !(length(clName) && any(nzchar(clName)))
    if (emptyName) {
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
#' @param clName The name of the clusterization. If not given the last available
#'   clusterization will be returned, as it is probably the most significant!
#'
#' @returns `getClusterizationData()` returns a `list` with 2 elements:
#'   * "clusters" the named cluster labels array
#'   * "coex" the associated `COEX` `data.frame`; it will be **empty** if not
#'     defined
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
    emptyName <- !(length(internalName) && any(nzchar(internalName)))
    assert_that(!emptyName,
                msg = "No clusterizations are present in the 'COTAN' object")

    clusters <- set_names(getMetadataCells(objCOTAN)[[internalName]],
                          getCells(objCOTAN))

    return(list("clusters" = clusters,
                "coex" = getClustersCoex(objCOTAN)[[internalName]]))
  }
)


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

