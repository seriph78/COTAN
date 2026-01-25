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
#'   stopifnot(identical(getDims(newObj), getDims(obj)))
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
    !is_empty(grep(x = colnames(genesMeta), pattern = "Genes?Names?",
                   ignore.case = TRUE))
  if (!hasGenesNameCol) {
    genesMeta <- rownames_to_column(genesMeta, var = "GenesNames")
    rownames(genesMeta) <- getGenes(objCOTAN)
  }

  cellsMeta <- getMetadataCells(objCOTAN)

  hasCellsIDCol <-
    !is_empty(grep(colnames(cellsMeta), pattern = "Cells?IDs?",
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
#' @param clNamesPattern A regular expression pattern used to identify the
#'   clustering columns in `colData`. Default supports `Seurat` conventions:
#'   `"^(COTAN_clusters_|seurat_clusters$|.*_snn_res\\..*|wsnn_res\\..*)"`
#' @param condNamesPattern A regular expression pattern used to identify the
#'   condition columns in `colData`. Default supports `Seurat` conventions:
#'   `"^(COTAN_conditions_|condition$|orig.ident$)"`
#' @param genesNamesPattern A regular expression pattern (case insentitive) used
#'   to identify the genes' names column in `rowData`. It used only if no names
#'   are available from the data matrix or the genes' data-set. Default supports
#'   is: `"Genes?Names?`
#' @param cellsIDsPattern A regular expression pattern (case insensitive) used
#'   to identify the cells' names column in `colData`. It used only if no names
#'   are available from the data matrix or the cells' data-set. Default supports
#'   is: `"Cells?IDs?"`
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
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment rowData
#' @importFrom SingleCellExperiment colData
#'
#' @export
#'
#' @seealso [COTAN-class], [SingleCellExperiment::SingleCellExperiment]
#'
#' @rdname Conversions
#'
convertFromSingleCellExperiment <- function(
    objSCE,
    clNamesPattern = "",
    condNamesPattern = "",
    genesNamesPattern = "",
    cellsIDsPattern = ""
  ) {
  assert_that(is(objSCE, "SingleCellExperiment"),
              msg = "Input object should be of type `SingleCellExperiment`.")

  # Extract counts matrix
  assert_that("counts" %in% SummarizedExperiment::assayNames(objSCE),
              msg = paste("The required 'counts' assay is missing",
                          "in the SingleCellExperiment object."))

  rawMatrix <- counts(objSCE)

  # Extract gene metadata
  genesMetaSCE <- rowData(objSCE)

  # Sanity check: same order and length
  assert_that(is.null(rownames(genesMetaSCE)) || is.null(rownames(rawMatrix)) ||
                identical(rownames(rawMatrix), rownames(genesMetaSCE)),
              msg = "Row names of counts and gene metadata differ.")

  genesNames <- ifelse(is_empty(rownames(rawMatrix)),
                       rownames(rawMatrix),
                       rownames(genesMetaSCE))

  # retrieve gene names from appropriate meta column
  if (is_empty(genesNames)) {
    genesNamesPattern <- ifelse(isEmptyName(genesNamesPattern),
                                "Genes?Names?",
                                genesNamesPattern)

    colToUse <- grep(x = colnames(genesMetaSCE),
                     pattern = genesNamesPattern,
                     ignore.case = TRUE)

    if (!is_empty(colToUse)) {
      # use first match
      genesNames <- genesMetaSCE[, colToUse[[1L]]]
    }
  }
  if (is_empty(genesNames)) {
    warning("Could not retrieve genes' names from the input")
  } else {
    rownames(rawMatrix)    <- make.unique(genesNames)
    rownames(genesMetaSCE) <- make.unique(genesNames)
  }

  # Extract cell metadata
  cellsMetaSCE <- colData(objSCE)

  # Sanity check: same order and length
  assert_that(is.null(rownames(cellsMetaSCE)) ||
                identical(colnames(rawMatrix), rownames(cellsMetaSCE)),
              msg = "Column names of counts and cell metadata differ.")

  cellsIDs <- ifelse(is_empty(colnames(rawMatrix)),
                     colnames(rawMatrix),
                     rownames(cellsMetaSCE))

  # retrieve cells IDs from appropriate meta column
  if (is_empty(cellsIDs)) {
    cellsIDsPattern <- ifelse(isEmptyName(cellsIDsPattern),
                              "CELLS?IDS?",
                              cellsIDsPattern)

    colToUse <- grep(x = colnames(cellsMetaSCE),
                     pattern = cellsIDsPattern,
                     ignore.case = TRUE)

        if (!is_empty(colToUse)) {
      cellsIDs <- cellsMetaSCE[, colToUse]
    }
  }
  if (is_empty(cellsIDs)) {
    warning("Could not retrieve cells' IDs from the input")
  } else {
    colnames(rawMatrix)    <- make.unique(cellsIDs)
    rownames(cellsMetaSCE) <- make.unique(cellsIDs)
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
  condCols <- grep("^(COTAN_conditions_|condition$|orig.ident$)",
                   names(cellsMetaSCE), value = FALSE)

  colsToKeep <- !(seq_len(ncol(cellsMetaSCE)) %in% union(clCols, condCols))
  cellsMeta <- cellsMetaSCE[, colsToKeep, drop = FALSE]

  # Attempt to retrieve co-expression matrices from metadata
  sceMetadataList <- S4Vectors::metadata(objSCE)

  genesCoex <- emptySymmetricMatrix()
  if ("genesCoex" %in% names(sceMetadataList)) {
    genesCoex <- sceMetadataList[["genesCoex"]]
  }

  cellsCoex <- emptySymmetricMatrix()
  if ("cellsCoex" %in% names(sceMetadataList)) {
    cellsCoex <- sceMetadataList[["cellsCoex"]]
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
    raw = as(as(rawMatrix, "Matrix"), "sparseMatrix"),
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
