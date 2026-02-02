# ---- Conversion to/from other classes ----

#' @title Data class conversions
#'
#' @description All functions to convert a [COTAN-class] object to/from other
#'   data classes used by the `BioConductor` analysis packages
#'
#' @name Conversions
#'
#' @examples
#'   data("test.dataset")
#'   obj <- COTAN(raw = test.dataset)
#'   obj <- proceedToCoex(obj, calcCoex = FALSE, saveObj = FALSE)
#'
#'   sce <- convertToSingleCellExperiment(objCOTAN = obj)
#'
#'   newObj <- convertFromSingleCellExperiment(sce)
#'
#'   stopifnot(identical(getDims(newObj)[-5:-6], getDims(obj)[-5:-6]),
#'             identical(unlist(getDims(newObj)[5:6]),
#'                       unlist(getDims(obj)[5:6]) + 1L))
#'
NULL

## ---- Convert to SingleCellExperiment ----

# convertToSingleCellExperiment

# @format `raw.dataset` is a data frame with \eqn{2000} genes and \eqn{815}
#   cells

#' @details `convertToSingleCellExperiment()` converts a [COTAN-class] object
#'   into a [SingleCellExperiment::SingleCellExperiment-class] object. Stores
#'   the raw counts in the `"counts"` [SummarizedExperiment::Assays-class], the
#'   metadata for genes and cells as `rowData` and `colData` slots respectively
#'   and finally the genes' and cells' `COEX` along the dataset metadata into
#'   the `metadata` slot.
#'
#'   The function performs the following steps:
#'   * Extracts the raw counts matrix, gene metadata, cell metadata, gene
#'     and cell *co-expression* matrix from the `COTAN`
#'     object; the `clustersCoex` slot is not converted
#'   * Identifies *clusterizations* and *conditions* in the cell metadata by the
#'     prefixes `"CL_"` and `"COND_"`
#'   * Renames *clusterization* columns with the prefix `"COTAN_clusters_"` and
#'     *condition* columns with the prefix `"COTAN_conditions_"`
#'   * Constructs a `SingleCellExperiment` object with the counts matrix, gene
#'     metadata, updated cell metadata, and stores the *co-expression* matrices
#'     in the `metadata` slot.
#'
#'   The resulting `SingleCellExperiment` object is compatible with downstream
#'   analysis packages and workflows within the Bioconductor ecosystem
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns A [SingleCellExperiment::SingleCellExperiment-class] object
#'   containing the data from the input [COTAN-class] object, with
#'   *clusterizations* and *conditions* appropriately prefixed and stored in the
#'   cell metadata.
#'
#' @importFrom methods is
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom tibble rownames_to_column
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
#'
#' @seealso [COTAN-class], [SingleCellExperiment::SingleCellExperiment]
#'
#' @rdname Conversions
#'

convertToSingleCellExperiment <- function(objCOTAN) {
  assert_that(is(objCOTAN, "COTAN"), validObject(objCOTAN),
              msg = "Input object should be of type `COTAN`.")

  # Identify clustering and condition columns
  genesMeta <- getMetadataGenes(objCOTAN)

  hasGenesNameCol <-
    !is_empty(grep(x = colnames(genesMeta), pattern = "Genes?[ -_]?Names?",
                   ignore.case = TRUE))
  if (!hasGenesNameCol) {
    genesMeta <- rownames_to_column(genesMeta, var = "GenesNames")
    rownames(genesMeta) <- getGenes(objCOTAN)
  }

  cellsMeta <- getMetadataCells(objCOTAN)

  hasCellsIDCol <-
    !is_empty(grep(colnames(cellsMeta), pattern = "Cells?[ -_]?IDs?",
                   ignore.case = TRUE))
  if (!hasCellsIDCol) {
    cellsMeta <- rownames_to_column(cellsMeta, var = "CellsIDs")
    rownames(cellsMeta) <- getCells(objCOTAN)
  }

  # Identify clustering and condition columns

  # Rename clustering columns replacing the 'CL_' prefix with 'COTAN_clusters_'
  colnames(cellsMeta) <-
    sub(x = colnames(cellsMeta), pattern = "^CL_",
        replacement = "COTAN_clusters_")

  # Rename conditions columns replacing the 'COND_' prefix with
  # 'COTAN_conditions_'
  colnames(cellsMeta) <-
    sub(x = colnames(cellsMeta), pattern = "^COND_",
        replacement = "COTAN_conditions_")

  # retrieves the genes' COEX if available
  genesCoex <- emptySymmetricMatrix()
  if (isCoexAvailable(objCOTAN, actOnCells = FALSE, ignoreSync = TRUE)) {
    genesCoex <- getGenesCoex(objCOTAN, zeroDiagonal = FALSE, ignoreSync = TRUE)
  }

  # retrieves the cells' COEX if available
  cellsCoex <- emptySymmetricMatrix()
  if (isCoexAvailable(objCOTAN, actOnCells = TRUE, ignoreSync = TRUE)) {
    cellsCoex <- getCellsCoex(objCOTAN, zeroDiagonal = FALSE, ignoreSync = TRUE)
  }

  # Create the SingleCellExperiment object
  objSCE <- SingleCellExperiment(
    assays = list(counts = getRawData(objCOTAN)),
    rowData = S4Vectors::DataFrame(genesMeta),
    colData = S4Vectors::DataFrame(cellsMeta),
    metadata = list(
      genesCoex = genesCoex,
      cellsCoex = cellsCoex,
      datasetMetadata = getMetadataDataset(objCOTAN)
    )
  )

  assert_that(identical(dim(objSCE),
                        c(getNumGenes(objCOTAN), getNumCells(objCOTAN))))

  return(objSCE)
}

## ---- Convert from SingleCellExperiment ----

#' convertFromSingleCellExperiment
#'
#' @details `convertFromSingleCellExperiment()` converts a
#'   [SingleCellExperiment::SingleCellExperiment-class] object back into a
#'   [COTAN-class] object. It supports `SCE` objects that were originally
#'   created from either a `COTAN` object or a `Seurat` object. The function
#'   extracts the `"counts"` matrix, genes' metadata, cells' metadata,
#'   *co-expression* matrices (if available), and reconstructs the `COTAN`
#'   object accordingly.

#' The function performs the following steps:
#'   * Extracts the raw matrix from the `"counts"`
#' [SummarizedExperiment::Assays-class]
#'   * Extracts gene metadata from `rowData`
#'   * Extracts cell metadata from `colData`, excluding any *clusterizations* or
#'     *conditions* present
#'   * Attempts to retrieve *co-expression* matrices from the `metadata` slot if
#' they exist
#'   * Constructs a `COTAN` object using the extracted data
#'   * Adds back the *clusterizations* and *conditions* using `COTAN` methods
#' If the COEX is not present (e.g., in `SCE` objects created from `Seurat`),
#' the `genesCoex` and `cellsCoex` slots in the resulting `COTAN` object will be
#' empty matrices
#'
#' @param objSCE A [SingleCellExperiment::SingleCellExperiment-class] object to
#'   be converted
#' @param assayName Name of the assay in `objSCE` to be used as raw counts.
#'   Default is `"counts"`. If not found, the function falls back to `"X"` and
#'   then to the first available assay, emitting a warning in both cases.
#' @param clNamesPattern A regular expression pattern used to identify the
#'   clustering columns in `colData`. Default supports `Seurat` conventions:
#'   `"^(COTAN_clusters_|seurat_clusters$|.*_snn_res\\..*|wsnn_res\\..*)"`
#' @param condNamesPattern A regular expression pattern used to identify the
#'   condition columns in `colData`. Default supports `Seurat` conventions:
#'   `"^(COTAN_conditions_|condition$|orig.ident$)"`
#' @param genesNamesPattern A regular expression pattern (case insensitive) used
#'   to identify the genes' names column in `rowData`. It used only if no names
#'   are available from the data matrix or the genes' data-set. Default supports
#'   is: `"^((Gene|Feature)s?[ -_]?(Name|ID)s?|IDs?|Symbols?)$"`
#' @param cellsIDsPattern A regular expression pattern (case insensitive) used
#'   to identify the cells' names column in `colData`. It used only if no names
#'   are available from the data matrix or the cells' data-set. Default supports
#'   is: `"^(Cells?[ -_]?IDs?|Barcodes?|Cells?)$"`
#'
#' @returns A [COTAN-class] object containing the data extracted from the input
#'   [SingleCellExperiment::SingleCellExperiment-class] object
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom methods new
#' @importFrom methods is
#' @importFrom methods validObject
#'
#' @importFrom rlang is_empty
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom Matrix dspMatrix
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment rowData
#' @importFrom SingleCellExperiment colData
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment altExp
#' @importFrom SingleCellExperiment altExpNames
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SingleCellExperiment reducedDimNames
#'
#' @export
#'
#' @seealso [COTAN-class], [SingleCellExperiment::SingleCellExperiment]
#'
#' @rdname Conversions
#'
convertFromSingleCellExperiment <- function(
    objSCE,
    assayName = "counts",
    clNamesPattern = "",
    condNamesPattern = "",
    genesNamesPattern = "",
    cellsIDsPattern = ""
  ) {
  assert_that(is(objSCE, "SingleCellExperiment"),
              msg = "Input object should be of type `SingleCellExperiment`.")

  # ---- Choose assay and extract raw matrix ----
  chooseAssay <- function(sceObj, assayName) {
    avail <- SummarizedExperiment::assayNames(sceObj)
    if (length(avail) < 1L) {
      stop("No assays present in the SingleCellExperiment object.")
    }
    if (!isEmptyName(assayName) && assayName %in% avail) {
      return(assayName)
    }
    # common fallbacks (e.g. zellkonverter uses 'X')
    fallback <- intersect(c("counts", "X"), avail)
    if (length(fallback) > 0L) {
        warning(sprintf("Assay '%s' not found; using '%s'.",
                        assayName, fallback[[1L]]))
        return(fallback[[1L]])
    }
    warning(sprintf("Assay '%s' not found; using first assay '%s'.",
                    assayName, avail[[1L]]))
    return(avail[[1L]])
  }

  rawSourceSCE <- objSCE
  chosenAssay <- chooseAssay(rawSourceSCE, assayName)

  # If requested, try altExp 'raw' when the desired assay is absent at top-level
  if (TRUE) {
    topAvail <- assayNames(objSCE)
    if (!(chosenAssay %in% topAvail)) {
      altNames <- tryCatch(altExpNames(objSCE),
                           error = function(e) character())
      if ("raw" %in% altNames) {
        rawSourceSCE <- altExp(objSCE, "raw")
        chosenAssay <- chooseAssay(rawSourceSCE, assayName)
      }
    }
  }

  rawMatrix <- assay(rawSourceSCE, chosenAssay)
  rawMatrix <- coerceToDgCMatrix(rawMatrix)

  validateRawCounts(rawMatrix)

  if (TRUE) {
    if (!identical(rawSourceSCE, objSCE)) {
      topAssays <- assayNames(objSCE)
      if (length(topAssays) > 0L) {
        warning("Dropping top-level assays: ",
                paste(topAssays, collapse = ", "))
      }
    }

    # assays
    keptAssay <- chosenAssay
    droppedAssays <- setdiff(assayNames(rawSourceSCE), keptAssay)
    if (length(droppedAssays) > 0L) {
      warning("Dropping assays: ", paste(droppedAssays, collapse = ", "))
    }

    # reducedDims (if available)
    redNames <- tryCatch(reducedDimNames(objSCE),
                         error = function(e) character())
    if (length(redNames) > 0L) {
      warning("Dropping reducedDims: ", paste(redNames, collapse = ", "))
    }

    # altExps
    altNames <- tryCatch(altExpNames(objSCE),
                         error = function(e) character())
    # don't warn about 'raw' if you actually used it
    if (identical(rawSourceSCE, objSCE)) {
      usedAlt <- character()
    } else {
      usedAlt <- "raw"
    }
    droppedAlt <- setdiff(altNames, usedAlt)
    if (length(droppedAlt) > 0L) {
      warning("Dropping altExps: ", paste(droppedAlt, collapse = ", "))
    }
  }


  # Extract gene metadata
  genesMetaSCE <- rowData(rawSourceSCE)

  ftCol <- intersect(c("feature_type", "FeatureType", "Type", "type"),
                     colnames(genesMetaSCE))
  if (length(ftCol) > 0L) {
    ft <- as.character(genesMetaSCE[[ftCol[[1L]]]])
    ft <- ft[!is.na(ft)]
    if (length(unique(ft)) > 1L) {
      stop(
        "Multiple feature types detected in rowData ('", ftCol[[1L]],
        "'). Subset the SingleCellExperiment to RNA features before conversion."
      )
    }
  }

  # Sanity check: same order and length
  assert_that(is.null(rownames(genesMetaSCE)) || is.null(rownames(rawMatrix)) ||
                identical(rownames(rawMatrix), rownames(genesMetaSCE)),
              msg = "Row names of counts and gene metadata differ.")

  if (!is_empty(rownames(rawMatrix))) {
    genesNames <- rownames(rawMatrix)
  } else {
    genesNames <- rownames(genesMetaSCE)
  }

  # retrieve gene names from appropriate meta column
  if (is_empty(genesNames)) {
    genesNamesPattern <-
      ifelse(
        isEmptyName(genesNamesPattern),
        # broadened default for common SCE conventions
        # (10X / scater / zellkonverter)
        "^((Gene|Feature)s?[ -_]?(Name|ID)s?|IDs?|Symbols?)$",
        genesNamesPattern)

    colToUse <- grep(x = colnames(genesMetaSCE),
                     pattern = genesNamesPattern,
                     ignore.case = TRUE)

    if (!is_empty(colToUse)) {
      # use first match
      genesNames <- genesMetaSCE[, colToUse[[1L]], drop = TRUE]
    }
  }
  if (is_empty(genesNames)) {
    warning("Could not retrieve genes' names from the input")
  } else {
    genesNames <- make.unique(as.character(genesNames))
    rownames(rawMatrix)    <- genesNames
    rownames(genesMetaSCE) <- genesNames
  }
  genesMeta <- genesMetaSCE

  # Extract cell metadata
  cellsMetaSCE <- colData(rawSourceSCE)

  # Sanity check: same order and length
  assert_that(is.null(rownames(cellsMetaSCE)) ||
                identical(colnames(rawMatrix), rownames(cellsMetaSCE)),
              msg = "Column names of counts and cell metadata differ.")

  if (!is_empty(colnames(rawMatrix))) {
    cellsIDs <- colnames(rawMatrix)
  } else {
    cellsIDs <- rownames(cellsMetaSCE)
  }

  # retrieve cells IDs from appropriate meta column
  if (is_empty(cellsIDs)) {
    cellsIDsPattern <-
      ifelse(
        isEmptyName(cellsIDsPattern),
        # broadened default for common SCE conventions
        "^(Cells?[ -_]?IDs?|barcodes?|cells?)$",
        cellsIDsPattern)

    colToUse <- grep(x = colnames(cellsMetaSCE),
                     pattern = cellsIDsPattern,
                     ignore.case = TRUE)

    if (!is_empty(colToUse)) {
      cellsIDs <- cellsMetaSCE[, colToUse[[1L]], drop = TRUE]
    }
  }
  if (is_empty(cellsIDs)) {
    warning("Could not retrieve cells' IDs from the input")
  } else {
    cellsIDs <- make.unique(as.character(cellsIDs))
    colnames(rawMatrix)    <- cellsIDs
    rownames(cellsMetaSCE) <- cellsIDs
  }

  # Identify clustering columns according to given pattern
  # (e.g., 'COTAN_clusters_', 'seurat_clusters', ...)

  if (isEmptyName(clNamesPattern)) {
    clNamesPattern <- paste0("^(COTAN_clusters_|seurat_clusters$|",
                             ".*_snn_res\\..*|wsnn_res\\..*|leiden_res\\..*)")
  }
  clCols <- grep(clNamesPattern, names(cellsMetaSCE), value = FALSE)

  # Identify condition columns
  # (e.g., 'COTAN_conditions_', 'condition', 'orig.ident', etc.)

  if (isEmptyName(condNamesPattern)) {
    condNamesPattern <- "^(COTAN_conditions_|condition$|orig.ident$)"
  }
  condCols <- grep(condNamesPattern, names(cellsMetaSCE), value = FALSE)

  colsToKeep <- !(seq_len(ncol(cellsMetaSCE)) %in% union(clCols, condCols))
  cellsMeta <- cellsMetaSCE[, colsToKeep, drop = FALSE]

  # Attempt to retrieve co-expression matrices from metadata
  sceMetadataList <- S4Vectors::metadata(objSCE)

  keptMd <- c("genesCoex", "cellsCoex", "datasetMetadata", "project.name")
  droppedMd <- setdiff(names(sceMetadataList), keptMd)
  if (length(droppedMd) > 0L) {
    warning("Dropping metadata entries: ", paste(droppedMd, collapse = ", "))
  }

  genesCoex <- emptySymmetricMatrix()
  if ("genesCoex" %in% names(sceMetadataList)) {
    genesCoex <- sceMetadataList[["genesCoex"]]
  }

  if (!is_empty(genesCoex)) {
    gd <- dim(genesCoex)
    if (length(gd) != 2L ||
        gd[[1L]] != nrow(rawMatrix) ||
        gd[[2L]] != nrow(rawMatrix)) {
      warning(paste("genesCoex dimensions do not match nrow(counts);",
              "dropping genesCoex."))
      genesCoex <- emptySymmetricMatrix()
    }
    genesCoex <-
      tryCatch(as(genesCoex, "dspMatrix"),
               error = function(e) {
                 warning("genesCoex cannot be made into a symmetric matrix")
                 emptySymmetricMatrix()
               })
  }

  cellsCoex <- emptySymmetricMatrix()
  if ("cellsCoex" %in% names(sceMetadataList)) {
    cellsCoex <- sceMetadataList[["cellsCoex"]]
  }

  if (!is_empty(cellsCoex)) {
    cd <- dim(cellsCoex)
    if (length(cd) != 2L ||
        cd[[1L]] != ncol(rawMatrix) ||
        cd[[2L]] != ncol(rawMatrix)) {
      warning(paste("cellsCoex dimensions do not match ncol(counts);",
                    "dropping cellsCoex."))
      cellsCoex <- emptySymmetricMatrix()
    }

    cellsCoex <-
      tryCatch(as(cellsCoex, "dspMatrix"),
               error = function(e) {
                 warning("cellsCoex cannot be made into a symmetric matrix")
                 emptySymmetricMatrix()
               })
  }

  # Extract dataset metadata if available
  datasetMeta <- data.frame()
  if ("datasetMetadata" %in% names(sceMetadataList)) {
    datasetMeta <- sceMetadataList[["datasetMetadata"]]
  } else if ("project.name" %in% names(sceMetadataList)) {
    # If the project name is stored directly in metadata
    datasetMeta <- data.frame(project.name = sceMetadataList[["project.name"]],
                              stringsAsFactors = FALSE)
  }

  # Create the COTAN object
  objCOTAN <- new(
    "COTAN",
    raw = rawMatrix,
    genesCoex = genesCoex,
    cellsCoex = cellsCoex,
    metaDataset = as.data.frame(datasetMeta),
    metaGenes = as.data.frame(genesMeta),
    metaCells = as.data.frame(cellsMeta),
    clustersCoex = list()
  )

  # Validate the COTAN object
  validObject(objCOTAN)

  if (!is_empty(clCols)) {
    logThis(paste("Clusterizations found:",
                  toString(names(cellsMetaSCE)[clCols])), logLevel = 2L)

    # Add clusterizations to the COTAN object
    for (col in names(cellsMetaSCE)[clCols]) {
      clusters <- getColumnFromDF(cellsMetaSCE, col)
      # Determine the clusterization name
      if (startsWith(col, "COTAN_clusters_")) {
        clName <- sub("^COTAN_clusters_", "", col)
      } else {
        clName <- col
      }
      # Add the clusterization using the method
      objCOTAN <- addClusterization(objCOTAN, clName = clName,
                                    clusters = clusters, override = FALSE)
    }
  }

  if (!is_empty(condCols)) {
    logThis(paste("Conditions found:",
                  toString(names(cellsMetaSCE)[condCols])), logLevel = 2L)

    # Add conditions to the COTAN object
    for (col in names(cellsMetaSCE)[condCols]) {
      conditions <- getColumnFromDF(cellsMetaSCE, col)
      # Determine the condition name
      if (startsWith(col, "COTAN_conditions_")) {
        condName <- sub("^COTAN_conditions_", "", col)
      } else {
        condName <- col
      }
      # Add the conditions using the method
      objCOTAN <- addCondition(objCOTAN, condName = condName,
                               conditions = conditions, override = FALSE)
    }
  }

  # Validate the COTAN object
  validObject(objCOTAN)

  return(objCOTAN)
}
