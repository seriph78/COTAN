# COTAN objects (read-only) accessors

#' getRawData
#'
#' @description This function extracts the raw count table.
#'
#' @param objCOTAN A COTAN object
#'
#' @returns the raw count sparse matrix
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
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
#' This function extracts the number of cells in the sample (𝑚)
#'
#' @param objCOTAN A COTAN object
#'
#' @return the number of cells in the sample (𝑚).
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' numCells <- getNumCells(objCOTAN)
#'
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
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' numGenes <- getNumGenes(objCOTAN)
#'
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
#' @return a cell character array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' cellsNames <- getCells(objCOTAN)
#'
#' @rdname getCells
setMethod(
  "getCells",
  "COTAN",
  function(objCOTAN) {
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
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' genesNames <- getGenes(objCOTAN)
#'
#' @rdname getGenes
setMethod(
  "getGenes",
  "COTAN",
  function(objCOTAN) {
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
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' zeroOne <- getZeroOneProj(objCOTAN)
#'
#' @rdname getZeroOneProj
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
#' @param objCOTAN A COTAN object
#'
#' @return an array with the library sizes
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' cellsSize <- getCellsSize(objCOTAN)
#'
#' @rdname getCellsSize
setMethod(
  "getCellsSize", "COTAN",
  function(objCOTAN) {
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
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' rawNorm <- getNormalizedData(objCOTAN)
#'
#' @rdname getNormalizedData
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
#' This function extract the meta-data stored for the data-set.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the meta-data data.frame
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' dataSetInfo <- getMetadataDataset(objCOTAN)
#'
#' @rdname getMetadataDataset
setMethod(
  "getMetadataDataset",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaDataset)
  }
)


#' getMetadataGenes
#'
#' This function extract the meta-data stored for the genes.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the meta-data data.frame
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' metaGenes <- getMetadataGenes(objCOTAN)
#'
#' @rdname getMetadataGenes
setMethod(
  "getMetadataGenes",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaGenes)
  }
)


#' getMetadataCells
#'
#' This function extract the meta-data stored for the cells.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the meta-data data.frame
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' metaCells <- getMetadataCells(objCOTAN)
#'
#' @rdname getMetadataCells
setMethod(
  "getMetadataCells",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@metaCells)
  }
)


#' getClustersCoex
#'
#' This function extract the complete clusterCoex list
#'
#' @param objCOTAN A COTAN object
#'
#' @return the list with a coex data.frame for each clusterization
#' When not empty, each data.frame contains a coex column for each cluster.
#'
#' @export
#'
#' @examples
#' allClustersCoexDF <- getClustersCoex(objCOTAN)
#'
#' @rdname getClustersCoex
setMethod(
  "getClustersCoex",
  "COTAN",
  function(objCOTAN) {
    return(objCOTAN@clustersCoex)
  }
)


#' getNu
#'
#' This function extract the nu array.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the nu array.
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- estimateNuLinear(objCOTAN)
#' nu <- getNu(objCOTAN)
#'
#' @rdname getNu
setMethod(
  "getNu",
  "COTAN",
  function(objCOTAN) {
    nu <- getMetadataCells(objCOTAN)[["nu"]]

    if (is_empty(nu)) {
      warning("nu is empty")
    }

    names(nu) <- getCells(objCOTAN)

    return(nu)
  }
)


#' getLambda
#'
#' This function extract the lambda array (mean expression for each gene).
#'
#' @param objCOTAN A COTAN object
#'
#' @return the lambda array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- estimateLambda(objCOTAN)
#' lambda <- getLambda(objCOTAN)
#'
#' @rdname getLambda
setMethod(
  "getLambda",
  "COTAN",
  function(objCOTAN) {
    lambda <- getMetadataGenes(objCOTAN)[["lambda"]]

    if (is_empty(lambda)) {
      warning("lambda is empty")
    }

    names(lambda) <- getGenes(objCOTAN)

    return(lambda)
  }
)


#' getDispersion
#'
#' This function extract the a array.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the dispersion array
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)
#' objCOTAN <- estimateDispersion(objCOTAN, cores = 12)
#' dispersion <- getDispersion(objCOTAN)
#'
#' @rdname getDispersion
setMethod(
  "getDispersion",
  "COTAN",
  function(objCOTAN) {
    dispersion <- getMetadataGenes(objCOTAN)[["dispersion"]]

    if (is_empty(dispersion)) {
      warning("dispersion is empty")
    }

    names(dispersion) <- getGenes(objCOTAN)

    return(dispersion)
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
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- findHousekeepingGenes(objCOTAN)
#' goodPos <- flagNotHousekeepingGenes(objCOTAN)
#'
#' @rdname flagNotHousekeepingGenes
setMethod(
  "flagNotHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getMetadataGenes(objCOTAN)[["hkGenes"]])) {
      return(rep(TRUE, getNumGenes(objCOTAN)))
    }
    else {
      return(getMetadataGenes(objCOTAN)[["hkGenes"]] == 0)
    }
  }
)


#' getHousekeepingGenes
#'
#' This function return the genes expressed in all cells in the dataset.
#'
#' @param objCOTAN A COTAN object
#'
#' @return an array containing all genes expressed in all cells
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- findHousekeepingGenes(objCOTAN)
#' hkGenes <- getHousekeepingGenes(objCOTAN)
#'
#' @rdname getHousekeepingGenes
setMethod(
  "getHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    return(getGenes(objCOTAN)[!flagNotHousekeepingGenes(objCOTAN)])
  }
)


#' getGenesCoex
#'
#' This function extract a complete (or a partial after genes dropping)
#' genes' coex matrix from the COTAN object.
#'
#' @param objCOTAN A COTAN object
#' @param genes A vector of gene names. It will exclude any gene not on the list.
#' By defaults the function will keep all genes.
#'
#' @return the genes' coex values
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)[["objCOTAN"]]
#' objCOTAN <- estimateDispersion(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' genesCoex <- getGenesCoex(objCOTAN)
#'
#' @rdname getGenesCoex
setMethod(
  "getGenesCoex",
  "COTAN",
  function(objCOTAN, genes = c()) {
    if (is_empty(genes)) {
      return(objCOTAN@genesCoex)
    }
    else {
      return(objCOTAN@genesCoex[, genes])
    }
  }
)


#' getCellsCoex
#'
#' This function extract a complete (or a partial after cells dropping)
#' cells' coex matrix from the COTAN object.
#'
#' @param objCOTAN A COTAN object
#' @param cells A vector of cell names. It will exclude any cell not on the list.
#' By defaults the function will keep all cells.
#'
#' @return the cells' coex values
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)[["objCOTAN"]]
#' objCOTAN <- estimateDispersion(objCOTAN, cores = 12)
#' objCOTAN <- estimateNuBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#' cellsCoex <- getCellsCoex(objCOTAN)
#'
#' @rdname getCellsCoex
setMethod(
  "getCellsCoex",
  "COTAN",
  function(objCOTAN, cells = c()) {
    if (is_empty(cells)) {
      return(objCOTAN@cellsCoex)
    }
    else {
      return(objCOTAN@cellsCoex[, cells])
    }
  }
)


#' getClusterizations
#'
#' This function extract the list of clusterizations defined in the COTAN object.
#'
#' @param objCOTAN A COTAN object
#' @param dropNoCoex When TRUE drops the names from the clusterizations with
#' empty associated coex data.frame
#'
#' @return a vactor of clusterizations names without the 'CL_' prefix
#'
#' @export
#'
#' @examples
#' clusterizations <- getClusterizations(objCOTAN)
#'
#' @rdname getClusterizations
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
#' This function extract the asked clusterization column and its coex
#' data.frame from the COTAN object.
#'
#' @param objCOTAN A COTAN object
#' @param clusterizationName The name of the clusterization.
#' If not given the last available clusterization will be returned,
#' as it is probably the most significant!
#'
#' @return a list with 'clusters' and 'coex'
#'
#' @export
#'
#' @examples
#'
#' list[cls, clsCoex] <- getClusterizationData(objCOTAN, clusterizationName = "merged")
#'
#' @rdname getClusterizationData
setMethod(
  "getClusterizationData",
  "COTAN",
  function(objCOTAN, clusterizationName = NULL) {
    if (is_empty(clusterizationName)) {
      clusterizationName <- getClusterizations(objCOTAN)[length(getClusterizations(objCOTAN))]
    }
    if (is_empty(clusterizationName)) {
      logThis("No clusterizations are present in the 'COTAN' object", logLevel = 3)
      return( list("clusters" = c(), "coex" = data.frame()) )
    }
    # clName can still be empty if no clusterization was store in the objCOTAN
    clName <- clusterizationName
    if (!startsWith(clName, "CL_")) {
      clName <- paste0("CL_", clName)
    }

    if (!clName %in% getClusterizations(objCOTAN, keepPrefix = TRUE)) {
      logThis(paste("Asked clusterization", clusterizationName,
                    "not present in the 'COTAN' object"), logLevel = 3)

      return( list("clusters" = c(), "coex" = data.frame()) )
    }

    return( list("clusters" = getMetadataCells(objCOTAN)[[clName]],
                 "coex" = getClustersCoex(objCOTAN)[[clName]]) )
  }
)
