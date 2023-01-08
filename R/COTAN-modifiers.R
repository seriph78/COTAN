# EXTREMELY IMPORTANT: DO NOT CHANGE THE LIST UP TO ELEMENT 8!!

#' datasetTags
#'
#' @returns a named character array with the standard labels used in the
#'   `metaDataset` of the `COTAN` objects.
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- addElementToMetaDataset(objCOTAN,
#'                                     datasetTags()[["cells"]],
#'                                     ncols(raw.dataset))
#' meta <- getMetadataDataset(objCOTAN)
#'
#' @rdname datasetTags
#'
datasetTags <- function() {
  return(c("GEO"   = "GEO:"                                             # 1
          ,"sc.m"  = "scRNAseq method:"                                 # 2
          ,"cond"  = "condition sample:"                                # 3
          ,"cells" = "starting n. of cells:"                            # 4
          ,"gsync" = "genes' coex is in sync:"                          # 5
          ,"csync" = "cells' coex is in sync:"                          # 6
          ,"gbad"  = "genes' pairs fraction with small expected count:" # 7
          ,"cbad"  = "cells' pairs fraction with small expected count:" # 8
          ))
}


#' updateMetaInfo
#'
#' @description Internal function: updates an information `data.frame`
#'
#' @param meta the information `data.frame` to update
#' @param tag the tag of the line
#' @param value the value or the values to associate to the tag
#'
#' @returns the updated `data.frame`
#'
#' @importFrom rlang is_empty
#'
#' @noRd
#'
updateMetaInfo <- function(meta, tag, value) {
  # all values are converted to strings
  newLine <- c(tag, paste0(value))

  if (!is_empty(meta) && (tag %in% meta[[1]])) {
    # existing tag: update the value
    rowPos <- which(meta[[1]] %in% tag)
  } else {
    # new tag: add a new entry
    rowPos <- nrow(meta) + 1
  }

  meta[rowPos, seq_along(newLine)] <- newLine

  return(meta)
}


#' initializeMetaDataset
#'
#' @description Initialize meta-data data-set
#'
#' @param objCOTAN the `COTAN` object
#' @param GEO a code reporting the GEO identification or other specific data-set
#'   code
#' @param sequencingMethod a string reporting the method used for the sequencing
#' @param sampleCondition a string reporting the specific sample condition or
#'   time point
#'
#' @returns the given `COTAN` object with updated `metaDataset`
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' obj <- COTAN(raw = raw.dataset)
#' obj <- initializeMetaDataset(obj, GEO = "code", sequencingMethod = "10X",
#'                              sampleCondition = "mouse dataset")
#'
#' @rdname initializeMetaDataset
#'
setMethod(
  "initializeMetaDataset",
  "COTAN",
  function(objCOTAN, GEO, sequencingMethod, sampleCondition) {
    logThis("Initializing `COTAN` meta-data", logLevel = 2)

    tags <- datasetTags()

    meta <- objCOTAN@metaDataset

    meta <- updateMetaInfo(meta, tags[["GEO"  ]], GEO)
    meta <- updateMetaInfo(meta, tags[["sc.m" ]], sequencingMethod)
    meta <- updateMetaInfo(meta, tags[["cond" ]], sampleCondition)
    meta <- updateMetaInfo(meta, tags[["cells"]], getNumCells(objCOTAN))
    meta <- updateMetaInfo(meta, tags[["gsync"]], FALSE)
    meta <- updateMetaInfo(meta, tags[["csync"]], FALSE)

    objCOTAN@metaDataset <- meta

    return(objCOTAN)
  }
)


#' addElementToMetaDataset
#'
#' @description This function is used to add a line of information to the
#'   information `data.frame`. If the tag was already used it will update the
#'   associated value(s) instead.
#'
#' @param objCOTAN a `COTAN` object
#' @param tag the new information tag
#' @param value a value (or an array) containing the information
#'
#' @returns the updated `COTAN` object
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- addElementToMetaDataset(objCOTAN, "Test", c("These are ", "some values"))
#' meta <- getMetadataDataset(objCOTAN)
#'
#' @rdname addElementToMetaDataset
#'
setMethod(
  "addElementToMetaDataset",
  "COTAN",
  function(objCOTAN, tag, value) {
    objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset, tag, value)

    return(objCOTAN)
  }
)


#' findHousekeepingGenes
#'
#' @description Determines the housekeeping genes inside a `COTAN` object
#'
#' @details Houskeeping genes are those genes that are expressed in all cells
#'
#' @param objCOTAN the `COTAN` object
#'
#' @returns the given `COTAN` object with updated housekeeping genes'
#'   information
#'
#' @importFrom Matrix rowSums
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- findHousekeepingGenes(objCOTAN)
#'
#' @rdname findHousekeepingGenes
#'
setMethod(
  "findHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    zeroOne <- getZeroOneProj(objCOTAN)

    # flag the genes with positive UMI count in every cell
    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes,
                                        rowSums(zeroOne) == getNumCells(objCOTAN),
                                        "hkGenes", getGenes(objCOTAN))

    return(objCOTAN)
  }
)


#' findFullyExpressedCells
#'
#' @description Determines the fully expressed cells inside a `COTAN` object
#'
#' @details Fully expressed cells are those cells where all genes are expressed
#'
#' @param objCOTAN the `COTAN` object
#'
#' @returns the given `COTAN` object with updated fully expressed cells'
#'   information
#'
#' @importFrom Matrix colSums
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- findFullyExpressedCells(objCOTAN)
#'
#' @rdname findFullyExpressedCells
#'
setMethod(
  "findFullyExpressedCells",
  "COTAN",
  function(objCOTAN) {
    zeroOne <- getZeroOneProj(objCOTAN)

    # flag the cells with positive UMI count in every gene
    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells,
                                        colSums(zeroOne) == getNumGenes(objCOTAN),
                                        "feCells", getCells(objCOTAN))

    return(objCOTAN)
  }
)


#' dropGenesCells
#'
#' @description  This function remove an array of genes and/or cells from the
#'   current `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes an array of gene names
#' @param cells an array of cell names
#'
#' @returns a completely new `COTAN` object with the new raw data obtained after
#'   the indicated genes/cells were expunged. Only the meta-data for the
#'   data-set are kept, while the rest is dropped as no more relevant with the
#'   restricted matrix
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' genes.to.rem <- getGenes(objCOTAN)[grep('^MT', getGenes(objCOTAN))]
#' cells.to.rem <- getCells(objCOTAN)[which(getCellsSize(objCOTAN) == 0)]
#' objCOTAN <- dropGenesCells(objCOTAN, genes.to.rem, cells.to.rem)
#'
#' @rdname dropGenesCells
#'
setMethod(
  "dropGenesCells",
  "COTAN",
  function(objCOTAN, genes, cells) {
    if (any(!genes %in% getGenes(objCOTAN)) ||
        any(!cells %in% getCells(objCOTAN))) {
      stop("Asked to drop genes and/or cells that were not present in the 'COTAN' object")
    }

    genesPosToKeep <- which(!(getGenes(objCOTAN) %in% genes))
    cellsPosToKeep <- which(!(getCells(objCOTAN) %in% cells))

    if (length(genesPosToKeep) == getNumGenes((objCOTAN)) &&
        length(cellsPosToKeep) == getNumCells((objCOTAN))) {
      logThis("No genes/cells where dropped", logLevel = 2)
    }

    # as all estimates would be wrong, a completely new object is created
    # with the same meta data for the data-set as the original
    output <- COTAN(objCOTAN@raw[genesPosToKeep, cellsPosToKeep])
    output@metaDataset <- getMetadataDataset(objCOTAN)

    return(output)
  }
)


#' addClusterization
#'
#' @description This function adds a clusterization to the current `COTAN`
#'   object, by adding a new column in the `metaCells data.frame` and adding a
#'   new element in the `clustersCoex list` using the passed in `coex`
#'   `data.frame`.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName the name of the clusterization to be added It cannot match an
#'   existing name; in case use [dropClusterization()] first.
#' @param clusters a (factors) array of cluster names
#' @param coexDF a `data.frame` where each column indicates the coex for each
#'   (or some) of the clusters of the clusterization
#'
#' @returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' clusters <- sample(1:10, getNumCells(objCOTAN), replace = TRUE)
#' names(clusters) <- getCells(objCOTAN)
#' objCOTAN <- addClusterization(objCOTAN, "dummy", clusters)
#'
#' @rdname addClusterization
#'
setMethod(
  "addClusterization",
  "COTAN",
  function(objCOTAN, clName, clusters, coexDF = data.frame()) {
    internalName <- clName
    if (!startsWith(internalName, "CL_")) {
      internalName <- paste0("CL_", clName)
    }

    if (internalName %in% colnames(getMetadataCells(objCOTAN))) {
      stop(paste0("A clusterization with name '", clName, "' already exists."))
    }

    if (length(clusters) != getNumCells(objCOTAN)) {
      stop(paste0("The passed clusterization has the wrong number of elements [",
                  length(clusters), "] instead of the expected number of cells [",
                  getNumCells(objCOTAN), "]."))
    }
    if (!is_empty(coexDF) && !isa(coexDF, "data.frame")) {
      stop(paste0("'clusterCoex' is supposedly composed of data.frames.",
                  " A '", class(coexDF), "' was given instead for clusterization '",
                  clName, "'."))
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, clusters,
                                        internalName, getCells(objCOTAN))

    # this add a new entry in the list for the new name!
    objCOTAN@clustersCoex[[internalName]] <- coexDF

    validObject(objCOTAN)

    return(objCOTAN)
  }
)

#' addClusterizationCoex
#'
#' @description This function adds a clusterization coex `data.frame` to the
#'   current `COTAN` object. It requires the clusterization to be already
#'   present.
#'
#' @seealso [addClusterization()]
#'
#' @param objCOTAN a `COTAN` object
#' @param clName the name of an existing clusterization
#' @param coexDF a `data.frame` where each column indicates the coex for all, or
#'   just some of, the clusters of the clusterization
#'
#' @returns the updated `COTAN` object
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' clusters <- sample(1:10, getNumCells(objCOTAN), replace = TRUE)
#' names(clusters) <- getCells(objCOTAN)
#' objCOTAN <- addClusterization(objCOTAN, "dummy", clusters)
#' coexDF <- as.data.frame(matrix(0, ncol=10, nrow=getNumGenes(objCOTAN)))
#' colnames(coexDF) <- c(1:10)
#' rownames(coexDF) <- getGenes(objCOTAN)
#' objCOTAN <- addClusterizationCoex(objCOTAN, clName = "dummy", coexDF = coexDF)
#'
#' @rdname addClusterizationCoex
#'
setMethod(
  "addClusterizationCoex",
  "COTAN",
  function(objCOTAN, clName, coexDF) {
    if (!isa(coexDF, "data.frame")) {
      stop(paste0("'clusterCoex' is supposedly composed of data.frames.",
                  " A '", class(coexDF), "' was given instead for clusterization '",
                  clName, "'."))
    }

    internalName <- clName
    if (!startsWith(internalName, "CL_")) {
      internalName <- paste0("CL_", clName)
    }

    if (!internalName %in% names(getClustersCoex(objCOTAN))) {
      stop(paste0("A clusterization with name '", clName, "' does not exists."))
    }

    # this should not add any new elements to the list!
    objCOTAN@clustersCoex[[internalName]] <- coexDF

    return(objCOTAN)
  }
)


#' dropGenesCoex
#'
#' @description This function drops the `genesCoex` from the given `COTAN`
#'   object
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns The updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
#'                                          GEO = "test_GEO",
#'                                          sequencingMethod = "test_method",
#'                                          sampleCondition = "test",
#'                                          cores = 12,
#'                                          saveObj = FALSE)
#' objCOTAN <- dropGenesCoex(objCOTAN)
#'
#' @rdname dropGenesCoex
#'
setMethod(
  "dropGenesCoex",
  "COTAN",
  function(objCOTAN) {
    if (!is_empty(objCOTAN@genesCoex)) {
      objCOTAN@genesCoex <- emptySymmetricMatrix()
    }

    return(objCOTAN)
  }
)


#' dropCellsCoex
#'
#' @description This function drops the `cellsCoex` from the given `COTAN`
#'   object
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns The updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
#'                                          GEO = "test_GEO",
#'                                          sequencingMethod = "test_method",
#'                                          sampleCondition = "test",
#'                                          cores = 12,
#'                                          saveObj = FALSE)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#' objCOTAN <- dropCellsCoex(objCOTAN)
#'
#' @rdname dropCellsCoex
#'
setMethod(
  "dropCellsCoex",
  "COTAN",
  function(objCOTAN) {
    if (!is_empty(objCOTAN@cellsCoex)) {
      objCOTAN@cellsCoex <- emptySymmetricMatrix()
    }

    return(objCOTAN)
  }
)


#' dropClusterization
#'
#' @description This function drops a clusterization from the current `COTAN`
#'   object, by removing the corresponding column in the `metaCells data.frame`
#'   and in case dropping the corresponding `coex data.frame` from the
#'   `clustersCoex list`.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName the name of an existing clusterization.
#'
#' @returns the updated `COTAN` object
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' clusters <- sample(1:10, getNumCells(objCOTAN), replace = TRUE)
#' names(clusters) <- getCells(objCOTAN)
#' objCOTAN <- addClusterization(objCOTAN, "dummy", clusters)
#' objCOTAN <- dropClusterization(objCOTAN, "dummy")
#'
#' @rdname dropClusterization
#'
setMethod(
  "dropClusterization",
  "COTAN",
  function(objCOTAN, clName) {
    internalName <- clName
    if (!startsWith(internalName, "CL_")) {
      internalName <- paste0("CL_", clName)
    }

    if (!internalName %in% names(getClustersCoex(objCOTAN))) {
      stop(paste0("A clusterization with name '", clName, "' does not exists."))
    }

    keptCols <- !colnames(objCOTAN@metaCells) %in% internalName
    objCOTAN@metaCells <- objCOTAN@metaCells[, keptCols, drop = FALSE]

    # assign NULL to drop elements from list
    objCOTAN@clustersCoex[[internalName]] <- NULL

    return(objCOTAN)
  }
)
