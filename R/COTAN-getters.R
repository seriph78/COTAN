# `COTAN` objects (read-only) accessors

#' getRawData
#'
#' @description This function extracts the raw count table.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the raw count sparse matrix
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' rawData <- getRawData(objCOTAN)
#'
#' @rdname getRawData
setMethod(
  "getRawData",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@raw)
  }
)


#' getNumCells
#'
#' @description This function extracts the number of cells in the sample (𝑚)
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the number of cells in the sample (𝑚).
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' numCells <- getNumCells(objCOTAN)
#'
#' @rdname getNumCells
#'
setMethod(
  "getNumCells",
  "COTAN",
  function(objCOTAN){
    return(ncol(objCOTAN@raw))
  }
)

#' getNumGenes
#'
#' @description This function extracts the number of genes in the sample (𝑛)
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the number of genes in the sample (𝑛).
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' numGenes <- getNumGenes(objCOTAN)
#'
#' @rdname getNumGenes
#'
setMethod(
  "getNumGenes",
  "COTAN",
  function(objCOTAN){
    return(nrow(objCOTAN@raw))
  }
)


#' getCells
#'
#' @description This function extract all cells in the dataset.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns a character array with the cells' names
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' cellsNames <- getCells(objCOTAN)
#'
#' @rdname getCells
#'
setMethod(
  "getCells",
  "COTAN",
  function(objCOTAN) {
    return(colnames(objCOTAN@raw))
  }
)


#' getGenes
#'
#' @description This function extract all genes in the dataset.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns a character array with the genes' names
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' genesNames <- getGenes(objCOTAN)
#'
#' @rdname getGenes
#'
setMethod(
  "getGenes",
  "COTAN",
  function(objCOTAN) {
    return(rownames(objCOTAN@raw))
  }
)


#' getZeroOneProj
#'
#' @description This function extract the raw count table where any positive
#'   number has been replaced with 1
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the raw count matrix projected to 0/1
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' zeroOne <- getZeroOneProj(objCOTAN)
#'
#' @rdname getZeroOneProj
#'
setMethod(
  "getZeroOneProj",
  "COTAN",
  function(objCOTAN) {
    return(sign(objCOTAN@raw))
  }
)


#' getCellsSize
#'
#' This function extracts the cell raw library size.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @return an array with the library sizes
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' cellsSize <- getCellsSize(objCOTAN)
#'
#' @rdname getCellsSize
#'
setMethod(
  "getCellsSize",
  "COTAN",
  function(objCOTAN) {
    return(colSums(objCOTAN@raw))
  }
)


#' getGenesSize
#'
#' This function extracts the genes raw library size.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @return an array with the library sizes
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' genesSize <- getGenesSize(objCOTAN)
#'
#' @rdname getGenesSize
#'
setMethod(
  "getGenesSize",
  "COTAN",
  function(objCOTAN) {
    return(rowSums(objCOTAN@raw))
  }
)


#' getNormalizedData
#'
#' @description This function extracts the normalized count table.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the normalized count dataframe (i.e. divided by nu).
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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

    return( t(t(getRawData(objCOTAN)) * (1/getNu(objCOTAN))) )
  }
)


#' getMetadataDataset
#'
#' @description This function extract the meta-data stored for the data-set.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the meta-data data.frame
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' dataSetInfo <- getMetadataDataset(objCOTAN)
#'
#' @rdname getMetadataDataset
#'
setMethod(
  "getMetadataDataset",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaDataset)
  }
)


#' getMetadataElement
#'
#' @description This function extracts the value associated with the given tag
#'   if present or an empty string otherwise.
#'
#' @param objCOTAN A `COTAN` object
#' @param tag The tag associated to the wanted value
#'
#' @returns A string with the relevant value
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- initializeMetaDataset(objCOTAN, )
#' GEO <- getMetadataElement(objCOTAN, "GEO")
#'
#' @rdname getMetadataElement
#'
setMethod(
  "getMetadataElement",
  "COTAN",
  function(objCOTAN, tag) {
    meta <- getMetadataDataset(objCOTAN)

    if (is_empty(meta) || !(tag %in% meta[[1]])){
      out <- ""
    } else {
      rowPos <- which(meta[[1]] %in% tag)
      out <- meta[, -1, drop = FALSE][rowPos, ]
    }

    return(out)
  }
)

#' getMetadataGenes
#'
#' @description This function extract the meta-data stored for the genes.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the meta-data data.frame
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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
#' @param objCOTAN A `COTAN` object
#'
#' @returns the meta-data data.frame
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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


#' getClustersCoex
#'
#' @description This function extract the complete clusterCoex list
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the list with a coex data.frame for each clusterization When not
#'   empty, each data.frame contains a coex column for each cluster.
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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


#' getNu
#'
#' @description This function extract the nu array.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the nu array.
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- estimateNuLinear(objCOTAN)
#' nu <- getNu(objCOTAN)
#'
#' @rdname getNu
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


#' getLambda
#'
#' This function extract the lambda array (mean expression for each gene).
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the lambda array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- estimateLambdaLinear(objCOTAN)
#' lambda <- getLambda(objCOTAN)
#'
#' @rdname getLambda
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


#' getDispersion
#'
#' @description This function extract the a array.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns the dispersion array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#' dispersion <- getDispersion(objCOTAN)
#'
#' @rdname getDispersion
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


#' flagNotHousekeepingGenes
#'
#' @description This function returns a Boolean vector with TRUE for those genes
#'   that are not housekeeping.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns an array of Booleans with TRUE for genes that are not housekeeping
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- findHousekeepingGenes(objCOTAN)
#' goodPos <- flagNotHousekeepingGenes(objCOTAN)
#'
#' @rdname flagNotHousekeepingGenes
#'
setMethod(
  "flagNotHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getMetadataGenes(objCOTAN)[["hkGenes"]])) {
      return(set_names(rep(TRUE, getNumGenes(objCOTAN)), getGenes(objCOTAN)))
    }
    else {
      return(!getMetadataGenes(objCOTAN)[["hkGenes"]])
    }
  }
)


#' flagNotFullyExpressedCells
#'
#' @description This function returns a Boolean vector with TRUE for those cells
#'   that are not fully expressed
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns an array of Booleans with TRUE for cells that are not fully
#'   expressed
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- findFullyExpressedCells(objCOTAN)
#' goodPos <- flagNotFullyExpressedCells(objCOTAN)
#'
#' @rdname flagNotFullyExpressedCells
#'
setMethod(
  "flagNotFullyExpressedCells",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getMetadataCells(objCOTAN)[["feCells"]])) {
      return(set_names(rep(TRUE, getNumCells(objCOTAN)), getCells(objCOTAN)))
    }
    else {
      return(!getMetadataCells(objCOTAN)[["feCells"]])
    }
  }
)


#' getHousekeepingGenes
#'
#' @description This function returns the genes expressed in all cells of the
#'   dataset.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns an array containing all genes expressed in all cells
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- findHousekeepingGenes(objCOTAN)
#' hkGenes <- getHousekeepingGenes(objCOTAN)
#'
#' @rdname getHousekeepingGenes
#'
setMethod(
  "getHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    return(getGenes(objCOTAN)[!flagNotHousekeepingGenes(objCOTAN)])
  }
)


#' getFullyExpressedCells
#'
#' @description This function return the cells that did express all genes of the
#'   dataset.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns an array containing all genes expressed in all cells
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- findFullyExpressedCells(objCOTAN)
#' feCells <- getFullyExpressedCells(objCOTAN)
#'
#' @rdname getFullyExpressedCells
#'
setMethod(
  "getFullyExpressedCells",
  "COTAN",
  function(objCOTAN) {
    return(getCells(objCOTAN)[!flagNotFullyExpressedCells(objCOTAN)])
  }
)


#' getGenesCoex
#'
#' @description This function extract a complete (or a partial after genes
#'   dropping) genes' coex matrix from the `COTAN` object.
#'
#' @param objCOTAN A `COTAN` object
#' @param genes A vector of gene names. It will exclude any gene not on the
#'   list. By defaults the function will keep all genes.
#'
#' @return the genes' coex values
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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
  function(objCOTAN, genes = c()) {
    if (!is_empty(objCOTAN@genesCoex) &&
        isFALSE(as.logical(getMetadataElement(objCOTAN, datasetTags()[5]))) ) {
      stop(paste0("Cannot return genes' coex as the matrix is",
                  " out of sync with the other parameters.",
                  " It is still ossible to access the data directly!"))
    }

    if (is_empty(genes)) {
      return(objCOTAN@genesCoex)
    }
    else {
      return(objCOTAN@genesCoex[, genes, drop = FALSE])
    }
  }
)


#' getCellsCoex
#'
#' @description This function extract a complete (or a partial after cells
#'   dropping) cells' coex matrix from the `COTAN` object.
#'
#' @param objCOTAN A `COTAN` object
#' @param cells A vector of cell names. It will exclude any cell not on the
#'   list. By defaults the function will keep all cells.
#'
#' @returns the cells' coex values
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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
  function(objCOTAN, cells = c()) {
    if (!is_empty(objCOTAN@genesCoex) &&
        isFALSE(as.logical(getMetadataElement(objCOTAN, datasetTags()[6]))) ) {
      stop(paste0("Cannot return genes' coex as the matrix is",
                  " out of sync with the other parameters.",
                  " It is still ossible to access the data directly!"))
    }

    if (is_empty(cells)) {
      return(objCOTAN@cellsCoex)
    }
    else {
      return(objCOTAN@cellsCoex[, cells, drop = FALSE])
    }
  }
)


#' getClusterizations
#'
#' @description This function extract the list of clusterizations defined in the
#'   `COTAN` object.
#'
#' @param objCOTAN A `COTAN` object
#' @param dropNoCoex When TRUE drops the names from the clusterizations with
#'   empty associated coex data.frame
#'
#' @returns a vector of clusterizations names without the 'CL_' prefix
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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
    }
    else {
      out <- names(clsCoex)
    }

    # drop the internal 'CL_' prefix
    if (isFALSE(keepPrefix)) {
      out <- substring(out, 4)
    }

    return(out)
  }
)

#' getClusterizationData
#'
#' @description This function extract the asked clusterization column and its
#'   coex data.frame from the `COTAN` object.
#'
#' @param objCOTAN A `COTAN` object
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
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
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
#' clDataList <- getClusterizationData(objCOTAN, clName = "merged")
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
      clName <- getClusterizations(objCOTAN)[length(getClusterizations(objCOTAN))]
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
      stop(paste("Asked clusterization", clName,
                    "not present in the 'COTAN' object"))
    }

    clusters <- set_names(getMetadataCells(objCOTAN)[[internalName]],
                          getCells(objCOTAN))

    return( list("clusters" = clusters,
                 "coex" = getClustersCoex(objCOTAN)[[internalName]]) )
  }
)

#' getDims
#'
#' This function extracts the sizes of all slots of the `COTAN` object.
#'
#' @param objCOTAN A `COTAN` object
#'
#' @return a named `list` with the sizes of the slots.
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
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
