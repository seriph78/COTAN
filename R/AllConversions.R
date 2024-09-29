# ---------- Conversion to/from other classes ----------

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
#'   identical(getDims(newObj), getDims(obj))
#'
NULL

# convertToSingleCellExperiment

# @format `raw.dataset` is a data frame with \eqn{2000} genes and \eqn{815}
#   cells

#' @details `convertToSingleCellExperiment()` converts a [COTAN-class] object
#'   into a [SingleCellExperiment-class] object. Stores the raw counts in the
#'   `"counts"` [Assays-class], the metadata for genes and cells as `rowData`
#'   and `colData` slots respectively and finally the genes' and cells' *coex*
#'   along the dataset metadata into the `metadata` slot.
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

#' @param objCOTAN a `COTAN` object
#'
#' @returns A [SingleCellExperiment-class] object containing the data from the
#'   input [COTAN-class] object, with clusterizations and conditions
#'   appropriately prefixed and stored in the cell metadata.
#'
#' @importFrom methods is
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
#' @seealso [COTAN-class], [SingleCellExperiment]
#'
#' @rdname Conversions
#'

convertToSingleCellExperiment <- function(objCOTAN) {
  assert_that(is(objCOTAN, "COTAN"), validObject(objCOTAN))

  # Identify clustering and condition columns
  cellsMeta <- getMetadataCells(objCOTAN)

  # Rename clustering columns replacing the "CL_" prefix with
  # 'COTAN_clusters_'
  clCols <- grep("^CL_", names(cellsMeta), value = FALSE)
  names(cellsMeta)[clCols] <-
    paste0("COTAN_clusters_", sub("^CL_", "", names(cellsMeta)[clCols]))

  # Rename conditions columns replacing the "COND_" prefix with
  # 'COTAN_conditions_'
  condCols <- grep("^COND_", names(cellsMeta), value = FALSE)
  names(cellsMeta)[condCols] <-
    paste0("COTAN_conditions_", sub("^COND_", "", names(cellsMeta)[condCols]))

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
    rowData = DataFrame(getMetadataGenes(objCOTAN)),
    colData = DataFrame(cellsMeta),
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


#' convertFromSingleCellExperiment
#'
#' @details `convertFromSingleCellExperiment()` converts a
#'   [SingleCellExperiment-class] object back into a [COTAN-class] object. It
#'   supports `SCE` objects that were originally created from either a `COTAN`
#'   object or a `Seurat` object. The function extracts the `"counts"` matrix,
#'   genes' metadata, cells' metadata, *co-expression* matrices (if available),
#'   and reconstructs the `COTAN` object accordingly.
#'
#'   The function performs the following steps:
#'   * Extracts the raw matrix from the `"counts"` [Assays-class]
#'   * Extracts gene metadata from `rowData`
#'   * Extracts cell metadata from `colData`, excluding any *clusterizations* or
#'     *conditions* present
#'   * Attempts to retrieve *co-expression* matrices from the `metadata` slot if
#'     they exist
#'   * Constructs a `COTAN` object using the extracted data
#'   * Adds back the *clusterizations* and *conditions* using `COTAN` methods
#'   If the COEX is not present (e.g., in `SCE` objects created from `Seurat`),
#'   the `genesCoex` and `cellsCoex` slots in the resulting `COTAN` object will
#'   be empty matrices
#'
#' @param objSCE A [SingleCellExperiment-class] object to be converted
#' @param clNamesPattern A regular expression pattern used to identify the
#'   clustering columns in `colData`. Default supports `Seurat` conventions:
#'   `"^(COTAN_clusters_|seurat_clusters$|.*_snn_res\\..*|wsnn_res\\..*)"`
#'
#' @returns A [COTAN-class] object containing the data extracted from the input
#'   [SingleCellExperiment-class] object
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom methods new
#' @importFrom methods is
#' @importFrom methods validObject
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment rowData
#' @importFrom SingleCellExperiment colData
#'
#' @importFrom SummarizedExperiment assayNames
#'
#' @importFrom S4Vectors metadata
#'
#' @export
#'
#' @seealso [COTAN-class], [SingleCellExperiment]
#'
#' @rdname Conversions
#'
convertFromSingleCellExperiment <- function(objSCE,
                                            clNamesPattern = "") {
  assert_that(is(objSCE, "SingleCellExperiment"),
              msg = "The input object must be of class 'SingleCellExperiment'.")

  # Extract counts matrix
  assert_that("counts" %in% assayNames(objSCE),
    msg = "The 'counts' assay is missing in the SingleCellExperiment object.")

  rawMatrix <- counts(objSCE)

  # Extract gene metadata
  genesMeta <- rowData(objSCE)

  # Extract cell metadata
  cellsMetaSCE <- colData(objSCE)

  # Identify clustering columns according to given pattern
  # (e.g., 'COTAN_clusters_', 'seurat_clusters', ...)

  if (isEmptyName(clNamesPattern)) {
    clNamesPattern <- "^(COTAN_clusters_|seurat_clusters$|.*_snn_res\\..*|wsnn_res\\..*|leiden_res\\..*)"
  }
  clCols <- grep(clNamesPattern, names(cellsMetaSCE), value = FALSE)

  # Identify condition columns
  # (e.g., 'COTAN_conditions_', 'condition', 'orig.ident', etc.)
  condCols <- grep("^(COTAN_conditions_|condition$|orig.ident$)",
                   names(cellsMetaSCE), value = FALSE)

  colsToKeep <- !(seq_len(ncol(cellsMetaSCE)) %in% union(clCols, condCols))
  cellsMeta <- cellsMetaSCE[, colsToKeep, drop = FALSE]

  # Attempt to retrieve co-expression matrices from metadata
  sceMetadataList <- metadata(objSCE)

  genesCoex <- emptySymmetricMatrix()
  if ("genesCoex" %in% names(sceMetadataList)) {
    genesCoex <- sceMetadataList$genesCoex
  }

  cellsCoex <- emptySymmetricMatrix()
  if ("cellsCoex" %in% names(sceMetadataList)) {
    cellsCoex <- sceMetadataList$cellsCoex
  }

  # Extract dataset metadata if available
  datasetMeta <- data.frame()
  if ("datasetMetadata" %in% names(sceMetadataList)) {
    datasetMeta <- sceMetadataList$datasetMetadata
  } else if ("project.name" %in% names(sceMetadataList)) {
    # If the project name is stored directly in metadata
    datasetMeta <- data.frame(project.name = sceMetadataList$project.name,
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
                  paste0(names(cellsMetaSCE)[clCols],
                  collapse = ", ")), logLevel = 2L)

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
                  paste0(names(cellsMetaSCE)[condCols],
                  collapse = ", ")), logLevel = 2L)

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

  return(objCOTAN)
}

