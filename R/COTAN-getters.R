# ------- `COTAN` raw data accessors --------

#' Raw data `COTAN` accessors
#'
#' @description These methods extract information out of a just created `COTAN`
#'   object. The accessors have **read-only** access to the object.
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


#' @details `getNumCells()` extracts the number of cells in the sample (𝑚)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNumCells()` returns the number of cells in the sample (𝑚).
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

#' @details `getNumGenes()` extracts the number of genes in the sample (𝑛)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `getNumGenes()` returns the number of genes in the sample (𝑛).
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
#' dataSetInfo <- getMetadataDataset(objCOTAN)
#'
#' @rdname RawDataGetters
#'
setMethod(
  "getMetadataDataset",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaDataset)
  }
)


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
#' objCOTAN <- initializeMetaDataset(objCOTAN, GEO = "code",
#'                                   sequencingMethod = "10X",
#'                                   sampleCondition = "mouse dataset")
#'
#' GEO <- getMetadataElement(objCOTAN, "GEO")
#'
#' @rdname RawDataGetters
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


# ------- `COTAN` estimated data accessors ------

#' getNormalizedData
#'
#' @description This function extracts the normalized count table.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns the normalized count dataframe (i.e. divided by nu).
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- estimateNuLinear(objCOTAN)
#' rawNorm <- getNormalizedData(objCOTAN)
#'
#' @rdname getNormalizedData
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


#' getMetadataGenes
#'
#' @description This function extract the meta-data stored for the genes.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns the meta-data data.frame
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- estimateLambdaLinear(objCOTAN)
#' metaGenes <- getMetadataGenes(objCOTAN)
#'
#' @rdname getMetadataGenes
#'
setMethod(
  "getMetadataGenes",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaGenes)
  }
)


#' getMetadataCells
#'
#' @description This function extract the meta-data stored for the cells.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns the meta-data data.frame
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- estimateNuLinear(objCOTAN)
#' metaCells <- getMetadataCells(objCOTAN)
#'
#' @rdname getMetadataCells
#'
setMethod(
  "getMetadataCells",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaCells)
  }
)


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
#'   has been made quite easy to excelude the flagged genes/cells in the user
#'   calculations.
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

#' getGenesCoex
#'
#' @description This function extract a complete (or a partial after genes
#'   dropping) genes' coex matrix from the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes A vector of gene names. It will exclude any gene not on the
#'   list. By defaults the function will keep all genes.
#' @param zeroDiagonal When TRUE sets the diagonal to zero.
#' @param ignoreSync When `TRUE` ignores whether the `lambda`/`nu`/`dispersion`
#'   have been updated since the coex matrix was calculated.
#'
#' @return the genes' coex values
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
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' genesCoex <- getGenesCoex(objCOTAN)
#'
#' @rdname getGenesCoex
#'
setMethod(
  "getGenesCoex",
  "COTAN",
  function(objCOTAN, genes = c(), zeroDiagonal = TRUE, ignoreSync = FALSE) {
    if (!is_empty(objCOTAN@genesCoex) && isFALSE(ignoreSync) &&
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


#' getCellsCoex
#'
#' @description This function extract a complete (or a partial after cells
#'   dropping) cells' coex matrix from the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param cells A vector of cell names. It will exclude any cell not on the
#'   list. By defaults the function will keep all cells.
#' @param zeroDiagonal When `TRUE` sets the diagonal to zero.
#' @param ignoreSync When `TRUE` ignores whether the `lambda`/`nu`/`dispersion`
#'   have been updated since the coex matrix was calculated.
#'
#' @returns the cells' coex values
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
#' objCOTAN <- estimateNuBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#' cellsCoex <- getCellsCoex(objCOTAN)
#'
#' @rdname getCellsCoex
#'
setMethod(
  "getCellsCoex",
  "COTAN",
  function(objCOTAN, cells = c(), zeroDiagonal = TRUE, ignoreSync = FALSE) {
    if (!is_empty(objCOTAN@cellsCoex) && isFALSE(ignoreSync) &&
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

#' getClusterizations
#'
#' @description This function extract the list of clusterizations defined in the
#'   `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param dropNoCoex When `TRUE` drops the names from the clusterizations with
#'   empty associated coex `data.frame`
#' @param keepPrefix When `TRUE` returns the internal name of the
#'   clusterization: the one with the 'CL_' prefix.
#'
#' @returns a vector of clusterizations names, usually without the 'CL_' prefix
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' clusterizations <- getClusterizations(objCOTAN)
#'
#' @rdname getClusterizations
#'
setMethod(
  "getClusterizations",
  "COTAN",
  function(objCOTAN, dropNoCoex = FALSE, keepPrefix = FALSE) {
    validObject(objCOTAN)

    clsCoex <- getClustersCoex(objCOTAN)

    if (isTRUE(dropNoCoex)) {
      out <- names(clsCoex[!sapply(clsCoex, is_empty)])
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


#' getClusterizationData
#'
#' @description This function extract the asked clusterization column and its
#'   coex data.frame from the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the clusterization. If not given the last available
#'   clusterization will be returned, as it is probably the most significant!
#'
#' @returns a list with 'clusters' and 'coex'
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = test.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = FALSE,
#'                                          outDir = tempdir())
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#' coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)[["coex"]]
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "clusters",
#'                               clusters = clusters, coexDF = coexDF)
#'
#' clDataList <- getClusterizationData(objCOTAN, clName = "clusters")
#' clusters <- clDataList[[1]]
#' coexDF   <- clDataList[[2]]
#'
#' @rdname getClusterizationData
#'
setMethod(
  "getClusterizationData",
  "COTAN",
  function(objCOTAN, clName = NULL) {
    if (is_empty(clName)) {
      clName <-
        getClusterizations(objCOTAN)[length(getClusterizations(objCOTAN))]
    }
    if (is_empty(clName)) {
      stop("No clusterizations are present in the 'COTAN' object")
    }
    # clName can still be empty if no clusterization was store in the objCOTAN
    internalName <- clName
    if (!startsWith(internalName, "CL_")) {
      internalName <- paste0("CL_", clName)
    }

    if (!internalName %in% getClusterizations(objCOTAN, keepPrefix = TRUE)) {
      stop("Asked clusterization", clName, "not present in the 'COTAN' object")
    }

    clusters <- set_names(getMetadataCells(objCOTAN)[[internalName]],
                          getCells(objCOTAN))

    return(list("clusters" = clusters,
                "coex" = getClustersCoex(objCOTAN)[[internalName]]))
  }
)


#' getClustersCoex
#'
#' @description This function extract the complete clusterCoex list
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns the list with a coex data.frame for each clusterization When not
#'   empty, each data.frame contains a coex column for each cluster.
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' allClustersCoexDF <- getClustersCoex(objCOTAN)
#'
#' @rdname getClustersCoex
#'
setMethod(
  "getClustersCoex",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@clustersCoex)
  }
)


#' getDims
#'
#' This function extracts the sizes of all slots of the `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#'
#' @return a named `list` with the sizes of the slots.
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = test.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = FALSE,
#'                                          outDir = tempdir())
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#' coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)[["coex"]]
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "clusters",
#'                               clusters = clusters, coexDF = coexDF)
#' allSizes <- getDims(objCOTAN)
#'
#' @rdname getDims
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
