#' cellsHeatmapPlot
#'
#' @description Heatmap plot of the cells' coex matrix
#'
#' @param objCOTAN a `COTAN` object
#' @param cells Which cells to plot (all if no argument is given)
#' @param clusters Use this clusterization to select/reorder the cells to plot
#'
#' @returns the cells' coex `Heatmap` plot
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom ComplexHeatmap Heatmap
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- initializeMetaDataset(objCOTAN,
#'                                   GEO = "test_GEO",
#'                                   sequencingMethod = "test_method",
#'                                   sampleCondition = "test")
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#' hPlots <- cellsHeatmapPlot(objCOTAN)
#'
#' @rdname cellsHeatmapPlot
#'
cellsHeatmapPlot <- function(objCOTAN, cells = NULL, clusters = NULL) {
  coexMat <- as.matrix(getCellsCoex(objCOTAN))
  stopifnot("cells coex not found in the COTAN" <- !is_empty(coexMat))

  diag(coexMat) <- 0

  # if clustering is needed
  if (!is_empty(clusters)) {
    # identifier for each cluster
    clustersTags <- unique(clusters)

    # size of each cluster
    clustersList <- toClustersList(clusters)
    clustersSize <- sapply(clustersList, length)

    # cell names grouped by the identifier of the cluster to which they belong
    cellNames <- unlist(clustersList)
    coexMat <- coexMat[cellNames, cellNames]
    colnames(coexMat) <- cellNames
    rownames(coexMat) <- cellNames

    Heatmap(coexMat,
            border = TRUE,
            column_split = factor(rep(clustersTags, clustersSize),
                                  levels = clustersTags),
            row_split = factor(rep(clustersTags, clustersSize),
                               levels = clustersTags),
            cluster_rows = FALSE,
            cluster_columns = FALSE)
  } else {
    cells <- handleNamesSubsets(getCells(objCOTAN), cells)

    Heatmap(coexMat[cells, cells],
            cluster_rows = FALSE,
            cluster_columns = FALSE)
  }
}

