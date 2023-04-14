
# ------------ Metadata handling ------------

# EXTREMELY IMPORTANT: DO NOT CHANGE THE LIST UP TO ELEMENT 8!!

#' @details `datasetTags()` defines a list of short names associated to an
#'   enumeration. It also defines the relative long names as they appear in the
#'   meta-data
#'
#' @returns `datasetTags()` a named `character array` with the standard labels
#'   used in the `metaDataset` of the `COTAN` objects
#'
#' @export
#'
#' @rdname HandleMetaData
#'
datasetTags <- function() {
  return(c("GEO"   = "GEO:"                                             # 1
         , "sc.m"  = "scRNAseq method:"                                 # 2
         , "cond"  = "condition sample:"                                # 3
         , "cells" = "starting n. of cells:"                            # 4
         , "gsync" = "genes' coex is in sync:"                          # 5
         , "csync" = "cells' coex is in sync:"                          # 6
         , "gbad"  = "genes' pairs fraction with small expected count:" # 7
         , "cbad"  = "cells' pairs fraction with small expected count:" # 8
         ))
}


#' @details `updateMetaInfo()` is an internal function: updates an information
#'   `data.frame`
#'
#' @param meta the information `data.frame` to update
#' @param tag the tag of the line
#' @param value the value or the values to associate to the tag
#'
#' @returns `updateMetaInfo()` returns the updated `data.frame`
#'
#' @importFrom rlang is_empty
#'
#' @noRd
#'
updateMetaInfo <- function(meta, tag, value) {
  # all values are converted to strings
  newLine <- c(tag, paste0(value))

  if (!is_empty(meta) && (tag %in% meta[[1L]])) {
    # existing tag: update the value
    rowPos <- which(meta[[1L]] %in% tag)
  } else {
    # new tag: add a new entry
    rowPos <- nrow(meta) + 1L
  }

  meta[rowPos, seq_along(newLine)] <- newLine

  return(meta)
}


#' @details `initializeMetaDataset()` initializes meta-data data-set
#'
#' @param objCOTAN the `COTAN` object
#' @param GEO a code reporting the GEO identification or other specific data-set
#'   code
#' @param sequencingMethod a string reporting the method used for the sequencing
#' @param sampleCondition a string reporting the specific sample condition or
#'   time point
#'
#' @returns `initializeMetaDataset()` returns the given `COTAN` object with the
#'   updated `metaDataset`
#'
#' @export
#'
#' @rdname HandleMetaData
#'
setMethod(
  "initializeMetaDataset",
  "COTAN",
  function(objCOTAN, GEO, sequencingMethod, sampleCondition) {
    logThis("Initializing `COTAN` meta-data", logLevel = 2L)

    tags <- datasetTags()

    meta <- objCOTAN@metaDataset

    meta <- updateMetaInfo(meta, tags[["GEO"]],   GEO)
    meta <- updateMetaInfo(meta, tags[["sc.m"]],  sequencingMethod)
    meta <- updateMetaInfo(meta, tags[["cond"]],  sampleCondition)
    meta <- updateMetaInfo(meta, tags[["cells"]], getNumCells(objCOTAN))
    meta <- updateMetaInfo(meta, tags[["gsync"]], FALSE)
    meta <- updateMetaInfo(meta, tags[["csync"]], FALSE)

    objCOTAN@metaDataset <- meta

    return(objCOTAN)
  }
)


#' @details `addElementToMetaDataset()` is used to add a line of information to
#'   the meta-data `data.frame`. If the tag was already used it will update the
#'   associated value(s) instead
#'
#' @param objCOTAN a `COTAN` object
#' @param tag the new information tag
#' @param value a value (or an array) containing the information
#'
#' @returns `addElementToMetaDataset()` returns the updated `COTAN` object
#'
#' @export
#'
#' @rdname HandleMetaData
#'
setMethod(
  "addElementToMetaDataset",
  "COTAN",
  function(objCOTAN, tag, value) {
    objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset, tag, value)

    return(objCOTAN)
  }
)


# ------------ Raw data cleaning --------------

#' @details `findHousekeepingGenes()` determines the housekeeping genes inside
#'   the raw data
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `findHousekeepingGenes()` returns the given `COTAN` object with
#'   updated housekeeping genes' information
#'
#' @importFrom Matrix rowSums
#'
#' @export
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "findHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    zeroOne <- getZeroOneProj(objCOTAN)

    # flag the genes with positive UMI count in every cell
    objCOTAN@metaGenes <-
      setColumnInDF(objCOTAN@metaGenes,
                    rowSums(zeroOne) == getNumCells(objCOTAN),
                    "hkGenes", getGenes(objCOTAN))

    return(objCOTAN)
  }
)


#' @details `findFullyExpressedCells()` determines the fully expressed cells
#'   inside the raw data
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `findFullyExpressedCells()` returns the given `COTAN` object with
#'   updated fully expressed cells' information
#'
#' @importFrom Matrix colSums
#'
#' @export
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "findFullyExpressedCells",
  "COTAN",
  function(objCOTAN) {
    zeroOne <- getZeroOneProj(objCOTAN)

    # flag the cells with positive UMI count in every gene
    objCOTAN@metaCells <-
      setColumnInDF(objCOTAN@metaCells,
                    colSums(zeroOne) == getNumGenes(objCOTAN),
                    "feCells", getCells(objCOTAN))

    return(objCOTAN)
  }
)


#' @details `dropGenesCells()` removes an array of genes and/or cells from the
#'   current `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes an array of gene names
#' @param cells an array of cell names
#'
#' @returns `dropGenesCells()` returns a completely new `COTAN` object with the
#'   new raw data obtained after the indicated genes/cells were expunged. Only
#'   the meta-data for the data-set are kept, while the rest is dropped as no
#'   more relevant with the restricted matrix
#'
#' @export
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "dropGenesCells",
  "COTAN",
  function(objCOTAN, genes, cells) {
    if (!all(genes %in% getGenes(objCOTAN)) ||
        !all(cells %in% getCells(objCOTAN))) {
      stop("Asked to drop genes and/or cells",
           " that were not present in the 'COTAN' object")
    }

    genesPosToKeep <- which(!(getGenes(objCOTAN) %in% genes))
    cellsPosToKeep <- which(!(getCells(objCOTAN) %in% cells))

    if (length(genesPosToKeep) == getNumGenes((objCOTAN)) &&
        length(cellsPosToKeep) == getNumCells((objCOTAN))) {
      logThis("No genes/cells where dropped", logLevel = 2L)
    }

    # as all estimates would be wrong, a completely new object is created
    # with the same meta data for the data-set as the original
    output <- COTAN(objCOTAN@raw[genesPosToKeep, cellsPosToKeep])
    output@metaDataset <- getMetadataDataset(objCOTAN)

    return(output)
  }
)


# -------------- Clusterization handling ---------------

#' @details `addClusterization()` adds a *clusterization* to the current `COTAN`
#'   object, by adding a new column in the `metaCells` `data.frame` and adding a
#'   new element in the `clustersCoex` `list` using the passed in `COEX`
#'   `data.frame` or an empty `data.frame` if none were passed in.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName the name of the *clusterization* to be added. It cannot match
#'   an existing name unless `override = TRUE` is used
#' @param clusters a (factors) array of *cluster* **labels**
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#' @param override When `TRUE` silently allows overriding data for an existing
#'   *clusterization* name. Otherwise the default behavior will avoid potential
#'   data losses
#'
#' @returns `addClusterization()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "addClusterization",
  "COTAN",
  function(objCOTAN, clName, clusters,
           coexDF = data.frame(), override = FALSE) {
    internalName <- clName
    if (!startsWith(internalName, "CL_")) {
      internalName <- paste0("CL_", clName)
    }

    if (!override && (internalName %in% colnames(getMetadataCells(objCOTAN)))) {
      stop("A clusterization with name '", clName, "' already exists.")
    }

    if (length(clusters) != getNumCells(objCOTAN)) {
      stop("The passed clusterization has the wrong number of elements [",
           length(clusters), "] instead of the expected number of cells [",
           getNumCells(objCOTAN), "].")
    }
    if (!is_empty(coexDF) && !isa(coexDF, "data.frame")) {
      stop("'clusterCoex' is supposedly composed of data.frames.",
           " A '", class(coexDF), "' was given instead for clusterization '",
           clName, "'.")
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, clusters,
                                        internalName, getCells(objCOTAN))

    # this add a new entry in the list for the new name!
    objCOTAN@clustersCoex[[internalName]] <- coexDF

    validObject(objCOTAN)

    return(objCOTAN)
  }
)

#' @details `addClusterizationCoex()` adds a *clusterization* `COEX`
#'   `data.frame` to the current `COTAN` object. It requires the named
#'   *clusterization* to be already present.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName the name of an existing clusterization
#' @param coexDF a `data.frame` where each column indicates the `COEX` for all,
#'   or just some of, the clusters of the clusterization
#'
#' @returns `addClusterizationCoex()` returns the updated `COTAN` object
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "addClusterizationCoex",
  "COTAN",
  function(objCOTAN, clName, coexDF) {
    if (!isa(coexDF, "data.frame")) {
      stop("'clusterCoex' is supposedly composed of data.frames.",
           " A '", class(coexDF), "' was given instead for clusterization '",
           clName, "'.")
    }

    internalName <- clName
    if (!startsWith(internalName, "CL_")) {
      internalName <- paste0("CL_", clName)
    }

    if (!internalName %in% names(getClustersCoex(objCOTAN))) {
      stop("A clusterization with name '", clName, "' does not exists.")
    }

    # this should not add any new elements to the list!
    objCOTAN@clustersCoex[[internalName]] <- coexDF

    return(objCOTAN)
  }
)


#' @details `dropClusterization()` drops a *clusterization* from the current
#'   `COTAN` object, by removing the corresponding column in the `metaCells`
#'   `data.frame` and the corresponding `COEX` `data.frame` from the
#'   `clustersCoex` `list`.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName the name of an existing *clusterization*.
#'
#' @returns `dropClusterization()` returns the updated `COTAN` object
#'
#' @export
#'
#' @rdname HandlingClusterizations
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
      stop("A clusterization with name '", clName, "' does not exists.")
    }

    keptCols <- !colnames(objCOTAN@metaCells) %in% internalName
    objCOTAN@metaCells <- objCOTAN@metaCells[, keptCols, drop = FALSE]

    # assign NULL to drop elements from list
    objCOTAN@clustersCoex[[internalName]] <- NULL

    return(objCOTAN)
  }
)


# ------------ COEX erasers ---------------

#' @details `dropGenesCoex()` drops the `genesCoex` member from the given
#'   `COTAN` object
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `dropGenesCoex()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname CalculatingCOEX
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


#' @details `dropCellsCoex()` drops the `cellsCoex` member from the given
#'   `COTAN` object
#'
#' @param objCOTAN A `COTAN` object
#'
#' @returns `dropCellsCoex()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname CalculatingCOEX
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
