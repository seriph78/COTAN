
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


#' @aliases initializeMetaDataset
#'
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

#' @aliases addElementToMetaDataset
#'
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

#' @aliases findFullyExpressedGenes
#'
#' @details `findFullyExpressedGenes()` determines the fully-expressed genes
#'   inside the raw data
#'
#' @param objCOTAN a `COTAN` object
#' @param cellsThreshold any gene that is expressed in more cells than threshold
#'   times the total number of cells will be marked as **fully-expressed**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
#'
#' @returns `findFullyExpressedGenes()` returns the given `COTAN` object with
#'   updated **fully-expressed** genes' information
#'
#' @export
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "findFullyExpressedGenes",
  "COTAN",
  function(objCOTAN, cellsThreshold = 0.99) {

    threshold <- cellsThreshold * getNumCells(objCOTAN)
    feGenes <- getNumOfExpressingCells(objCOTAN) >= threshold

    objCOTAN@metaGenes <-
      setColumnInDF(objCOTAN@metaGenes, feGenes, "feGenes", getGenes(objCOTAN))

    return(objCOTAN)
  }
)


#' @aliases findFullyExpressingCells
#'
#' @details `findFullyExpressingCells()` determines the cells that are
#'   expressing all genes in the dataset
#'
#' @param objCOTAN a `COTAN` object
#' @param genesThreshold any cell that is expressing more genes than threshold
#'   times the total number of genes will be marked as **fully-expressing**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
#'
#' @returns `findFullyExpressingCells()` returns the given `COTAN` object  with
#'   updated **fully-expressing** cells' information
#'
#' @export
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "findFullyExpressingCells",
  "COTAN",
  function(objCOTAN, genesThreshold = 0.99) {

    threshold <- genesThreshold * getNumGenes(objCOTAN)
    feCells <- getNumExpressedGenes(objCOTAN) >= threshold

    objCOTAN@metaCells <-
      setColumnInDF(objCOTAN@metaCells, feCells, "feCells", getCells(objCOTAN))

    return(objCOTAN)
  }
)


#' @aliases storeGDI
#'
#' @details `storeGDI()` stored and already calculated genes' GDI `array` in a
#'   `COTAN` object. It can be retrieved using the method [getGDI()]
#'
#' @param objCOTAN a `COTAN` object
#' @param genesGDI the named genes' GDI `array` to store or the output
#'   `data.frame` of the function [calculateGDI()]
#'
#' @returns `storeGDI()` returns the given `COTAN` object  with updated
#'   **GDI** genes' information
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @rdname GenesStatistics
#'
setMethod(
  "storeGDI",
  "COTAN",
  function(objCOTAN, genesGDI) {
    if (isa(genesGDI, "data.frame")) {
      assert_that("GDI" %in% colnames(genesGDI),
                  msg = "Passed data.frame has no GDI column")
      genesGDI <- getColumnFromDF(genesGDI, colName = "GDI")
    }

    allGenes <- getGenes(objCOTAN)
    assert_that(setequal(names(genesGDI), allGenes),
                msg = paste("Passed GDI array must be named",
                            "and aligned to the object genes' list"))

    objCOTAN@metaGenes <-
      setColumnInDF(objCOTAN@metaGenes, genesGDI[allGenes], "GDI", allGenes)

    return(objCOTAN)
  }
)


#' @aliases dropGenesCells
#'
#' @details `dropGenesCells()` removes an array of genes and/or cells from the
#'   current `COTAN` object.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes an array of gene names
#' @param cells an array of cell names
#'
#' @returns `dropGenesCells()` returns a completely new `COTAN` object with the
#'   new raw data obtained after the indicated genes/cells were expunged. All
#'   remaining data is dropped too as no more relevant with the restricted
#'   matrix. Exceptions are:
#'   * the meta-data for the data-set that gets kept unchanged
#'   * the meta-data of genes/cells that gets restricted to the remaining
#'     elements. The columns calculated via `estimate` and `find` methods are
#'     dropped too
#'
#' @export
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom methods validObject
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "dropGenesCells",
  "COTAN",
  function(objCOTAN, genes, cells) {
    assert_that((all(genes %in% getGenes(objCOTAN)) &&
                 all(cells %in% getCells(objCOTAN))),
                msg = paste0("Asked to drop genes and/or cells",
                             " that were not present in the 'COTAN' object"))

    genesPosToKeep <- which(!(getGenes(objCOTAN) %in% genes))
    cellsPosToKeep <- which(!(getCells(objCOTAN) %in% cells))

    assert_that(length(genesPosToKeep) != 0L,
                length(cellsPosToKeep) != 0L,
                msg = "Asked to drop all genes and/or cells")

    if (length(genesPosToKeep) == getNumGenes(objCOTAN) &&
        length(cellsPosToKeep) == getNumCells(objCOTAN)) {
      logThis("Asked to drop no genes or cells", logLevel = 2L)
      # return original object
      return(objCOTAN)
    }

    logThis(paste("Asked to drop", length(genes), "genes and",
                  length(cells), "cells"), logLevel = 3L)

    # As all estimates would be wrong, a completely new object is created
    output <- COTAN(objCOTAN@raw[genesPosToKeep, cellsPosToKeep, drop = FALSE])

    # Copy the meta data for the data-set thus erasing the sync bits
    output@metaDataset <- getMetadataDataset(objCOTAN)
    output <- addElementToMetaDataset(output, tags[["gsync"]], FALSE)
    output <- addElementToMetaDataset(output, tags[["csync"]], FALSE)

    # Filter the meta data for genes keeping those not related to estimates
    colsToKeep <- !(names(getMetadataGenes(objCOTAN)) %in%
                      c("feGenes", "lambda", "dispersion", "GDI"))
    if (any(colsToKeep)) {
      output@metaGenes <-
        getMetadataGenes(objCOTAN)[genesPosToKeep, colsToKeep, drop = FALSE]
    }

    # Filter the meta data for cells keeping those not related to estimates
    colsToKeep <- !(names(getMetadataCells(objCOTAN)) %in%
                      c("feCells", "nu"))
    if (any(colsToKeep)) {
      output@metaCells <-
        getMetadataCells(objCOTAN)[cellsPosToKeep, colsToKeep, drop = FALSE]
    }

    # Drop all clusterizations' data.frames, but ensure object validity
    for (internalName in names(getClustersCoex(objCOTAN))) {
      output@clustersCoex[[internalName]] <- data.frame()
    }

    validObject(output)

    return(output)
  }
)


# -------------- Clusterization handling ---------------

#' @aliases addClusterization
#'
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
#' @importFrom assertthat assert_that
#'
#' @importFrom methods validObject
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

    clusters <- factor(clusters)

    assert_that(nchar(internalName) >= 4L,
                msg = "Given an empty name for the new clusterization.")

    assert_that(isTRUE(override) ||
                  !(internalName %in% colnames(getMetadataCells(objCOTAN))),
                msg = paste0("A clusterization with name '",
                             clName, "' already exists."))

    assert_that(length(clusters) == getNumCells(objCOTAN),
                msg = paste0("The passed clusterization has the ",
                             "wrong number of elements [", length(clusters),
                             "] instead of the expected number of cells [",
                             getNumCells(objCOTAN), "]."))

    assert_that(identical(names(clusters), getCells(objCOTAN)),
                msg = paste0("The passed clusterization must be named ",
                             "and aligned to the cells' list"))

    if (!is_empty(coexDF)) {
      assert_that(isa(coexDF, "data.frame"),
                  msg = paste0("'coexDF' is supposedly composed of ",
                               "data.frames. A '", class(coexDF),
                               "' was given  instead for clusterization '",
                               clName, "'."))

      assert_that(identical(rownames(coexDF), getGenes(objCOTAN)),
                  setequal(colnames(coexDF), levels(clusters)),
                  msg = "coex is not aligned to the given clusterization")
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, clusters,
                                        internalName, getCells(objCOTAN))

    # this add a new entry in the list for the new name!
    objCOTAN@clustersCoex[[internalName]] <- coexDF

    validObject(objCOTAN)

    return(objCOTAN)
  }
)

#' @aliases addClusterizationCoex
#'
#' @details `addClusterizationCoex()` adds a *clusterization* `COEX`
#'   `data.frame` to the current `COTAN` object. It requires the named
#'   *clusterization* to be already present.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName the name of an existing *clusterization*
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#'
#' @returns `addClusterizationCoex()` returns the updated `COTAN` object
#'
#' @export
#'
#' @importFrom assertthat assert_that
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "addClusterizationCoex",
  "COTAN",
  function(objCOTAN, clName, coexDF) {
    internalName <- getClusterizationName(objCOTAN, clName = clName,
                                          keepPrefix = TRUE)

    if (!is_empty(coexDF)) {
      assert_that(isa(coexDF, "data.frame"),
                  msg = paste0("'coexDF' is supposedly composed of ",
                               "data.frames. A '", class(coexDF),
                               "' was given  instead for clusterization '",
                               clName, "'."))

      assert_that(identical(rownames(coexDF), getGenes(objCOTAN)),
                  setequal(colnames(coexDF),
                           getClusters(objCOTAN, clName = internalName)),
                  msg = "coex is not aligned to the given clusterization")
    }

    # this should not add any new elements to the list!
    objCOTAN@clustersCoex[[internalName]] <- coexDF

    return(objCOTAN)
  }
)

#' @aliases dropClusterization
#'
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
#' @importFrom assertthat assert_that
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

    assert_that(internalName %in% colnames(getMetadataCells(objCOTAN)),
                msg = paste0("No clusterization with name '",
                             clName, "' is present."))

    keptCols <- !colnames(objCOTAN@metaCells) %in% internalName
    objCOTAN@metaCells <- objCOTAN@metaCells[, keptCols, drop = FALSE]

    # assign NULL to drop elements from list
    objCOTAN@clustersCoex[[internalName]] <- NULL

    return(objCOTAN)
  }
)


# -------------- Conditions handling ---------------

#' @aliases addCondition
#'
#' @details `addCondition()` adds a *condition* to the current `COTAN` object,
#'   by adding a new column in the `metaCells` `data.frame`
#'
#' @param objCOTAN a `COTAN` object
#' @param condName the name of the *condition* to be added. It cannot match an
#'   existing name unless `override = TRUE` is used
#' @param conditions a (factors) array of *condition* **labels**
#' @param override When `TRUE` silently allows overriding data for an existing
#'   *condition* name. Otherwise the default behavior will avoid potential
#'   data losses
#'
#' @returns `addCondition()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom methods validObject
#'
#' @export
#'
#' @rdname HandlingConditions
#'
setMethod(
  "addCondition",
  "COTAN",
  function(objCOTAN, condName, conditions, override = FALSE) {
    internalName <- condName
    if (!startsWith(internalName, "COND_")) {
      internalName <- paste0("COND_", condName)
    }

    assert_that(nchar(internalName) >= 6L,
                msg = "Given an empty name for the new condition")

    assert_that(isTRUE(override) ||
                  !(internalName %in% colnames(getMetadataCells(objCOTAN))),
                msg = paste0("A condition with name '",
                             condName, "' already exists."))

    assert_that(length(conditions) == getNumCells(objCOTAN),
                msg = paste0("The passed condition has the wrong ",
                             "number of elements [", length(conditions),
                             "] instead ofthe expected number of cells [",
                             getNumCells(objCOTAN), "]."))

    assert_that(identical(names(conditions), getCells(objCOTAN)),
                msg = paste0("The passed condition must be named ",
                             "and aligned to the cells' list"))

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, factor(conditions),
                                        internalName, getCells(objCOTAN))

    validObject(objCOTAN)

    return(objCOTAN)
  }
)


#' @aliases dropCondition
#'
#' @details `dropCondition()` drops a *condition* from the current `COTAN`
#'   object, by removing the corresponding column in the `metaCells`
#'   `data.frame`
#'
#' @param objCOTAN a `COTAN` object
#' @param condName the name of an existing *condition*.
#'
#' @returns `dropCondition()` returns the updated `COTAN` object
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @rdname HandlingConditions
#'
setMethod(
  "dropCondition",
  "COTAN",
  function(objCOTAN, condName) {
    internalName <- condName
    if (!startsWith(internalName, "COND_")) {
      internalName <- paste0("COND_", condName)
    }

    assert_that(internalName %in% colnames(getMetadataCells(objCOTAN)),
                msg = paste0("No condition with name '",
                             condName, "' is present."))

    keptCols <- !colnames(objCOTAN@metaCells) %in% internalName
    objCOTAN@metaCells <- objCOTAN@metaCells[, keptCols, drop = FALSE]

    return(objCOTAN)
  }
)


# ------------ COEX erasers ---------------

#' @aliases dropGenesCoex
#'
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


#' @aliases dropCellsCoex
#'
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
