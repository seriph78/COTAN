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
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- initializeMetaDataset(GEO = "test_GEO",
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
  coexMat <- getCellsCoex(objCOTAN)
  stopifnot("cells coex not found in the COTAN" <- !is_empty(coexMat))

  diag(coexMat) <- 0

  # if clustering is needed
  if (!is_empty(clusters)) {
    # identifier for each cluster
    clustersIdentifier <- unique(sort(clusters))

    # size of each cluster
    clustersSize <- sapply(clustersIdentifier, function(x) sum(clusters == x))

    # cell names ordered by the identifier of the cluster to which they belong
    nameSort <- names(sort(clusters))
    coexMat <- coexMat[nameSort, nameSort]

    colnames(coexMat) <- clusters[nameSort]
    rownames(coexMat) <- clusters[nameSort]

    Heatmap(coexMat,
            border = TRUE,
            column_split = factor(rep(clustersIdentifier, clustersSize),
                                  levels = title),
            row_split = factor(rep(clustersIdentifier, clustersSize),
                               levels = title),
            cluster_rows = FALSE,
            cluster_columns = FALSE)
  } else {
    cells <- handleNamesSubsets(getCells(objCOTAN), cells)

    Heatmap(coexMat[cells, cells],
            cluster_rows = FALSE,
            cluster_columns = FALSE)
  }
}

