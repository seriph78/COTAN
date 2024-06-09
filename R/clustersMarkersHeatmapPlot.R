
#' @details `clustersMarkersHeatmapPlot()` returns the heatmap plot of a summary
#'   score for each *cluster* and each gene marker list in the given
#'   *clusterization*. It also returns the numerosity and percentage of each
#'   *cluster* on the right and a gene *clusterization* dendogram on the left
#'   (as returned by the function [geneSetEnrichment()]) that allows to estimate
#'   which markers groups are more or less expressed in each *cluster* so it is
#'   easier to derive the *clusters*' cell types.
#'
#' @param objCOTAN a `COTAN` object
#' @param groupMarkers a named `list` with an element for each group comprised
#'   of one or more marker genes
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName` that, in such a case, will only indicate
#'   the relevant column name in the returned `data.frame`
#' @param kCuts the number of estimated *cluster* (this defines the height for
#'   the tree cut and the associated colors)
#' @param condNameList a `list` of *conditions*' names to be used for additional
#'   columns in the final plot. When none are given no new columns will be added
#'   using data extracted via the function [clustersSummaryData()]
#' @param conditionsList a `list` of *conditions* to use. If given they will
#'   take precedence on the ones indicated by `condNameList`
#'
#' @returns `clustersMarkersHeatmapPlot()` returns a list with:
#'  * `"heatmapPlot"` the complete heatmap plot
#'  * `"dataScore"` the `data.frame` with the score values
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom dendextend set
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr %>%
#'
#' @importFrom grid unit
#' @importFrom grid gpar
#' @importFrom grid grid.text
#'
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom ComplexHeatmap anno_barplot
#' @importFrom ComplexHeatmap anno_numeric
#' @importFrom ComplexHeatmap anno_text
#' @importFrom ComplexHeatmap Legend
#' @importFrom ComplexHeatmap Heatmap
#'
#' @importFrom circlize colorRamp2
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
clustersMarkersHeatmapPlot <- function(objCOTAN, groupMarkers,
                                       clName = "", clusters = NULL,
                                       kCuts = 3L,
                                       condNameList = NULL,
                                       conditionsList = NULL) {
  assert_that(is_empty(conditionsList) ||
                length(conditionsList) == length(condNameList),
              msg = "Explicitly given conditions must have corresponding names")

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  expressionCl <- clustersDeltaExpression(objCOTAN, clusters = clusters,
                                          clName = clName)
  scoreDF <- geneSetEnrichment(groupMarkers = groupMarkers,
                               clustersCoex = expressionCl)

  scoreDFT <- t(scoreDF[, 1L:(ncol(scoreDF) - 2L), drop = FALSE])
  dend <- clustersTreePlot(objCOTAN, kCuts = kCuts, clName = clName)[["dend"]]

  dend <- set(dend = dend, "branches_lwd", 2L)

  hbList <- NULL
  lgdList <- list()
  if (!is_empty(condNameList)) {
    allColors <- getColorsVector()

    for (i in seq_along(condNameList)) {
      # gets a dummy condition if none was given
      c(condName, conditions) %<-%
        normalizeNameAndLabels(objCOTAN, name = condNameList[[i]],
                               labels = conditionsList[[i]], isCond = TRUE)
      condDF <- clustersSummaryData(objCOTAN, clName = clName,
                                    clusters = clusters,
                                    condName = condName,
                                    conditions = conditions)

      cNames <- c(clName, condName, "CellNumber")

      condDF <- condDF[, cNames] %>% pivot_wider(names_from = condName,
                                                 values_from = cNames[[3L]])

      # Here we have a column per condition alternative plus 1
      nCols <- ncol(condDF)
      assert_that(nCols > 1L, msg = "Internal error: provided empty conditions")
      condDF[, 2L:nCols] <- condDF[, 2L:nCols] / rowSums(condDF[, 2L:nCols])

      condDF <- as.data.frame(condDF)
      rownames(condDF) <- condDF[[1L]]
      condDF <- condDF[rownames(scoreDFT), ]
      condDF[is.na(condDF)] <- 0L

      conds <- colnames(condDF)[2L:nCols]

      colrs <- set_names(head(allColors, nCols - 1L), conds)
      allColors <- tail(allColors, length(allColors) - nCols + 1L)

      hb <- HeatmapAnnotation(
        which = "row",
        condName = anno_barplot(condDF[, 2L:nCols],
                                width = unit(3.0, "cm"),
                                gp = gpar(fill = colrs, col = "black"),
                                align_to = "right",
                                labels_gp = gpar(fontsize = 12L)),
        annotation_label = condName,
        annotation_name_rot = 0L)
      hbList <- append(hbList, hb)

      lgd <- Legend(labels = conds, title = condName,
                    legend_gp = gpar(fill = colrs))
      lgdList <- append(lgdList, lgd)
    }
  }

  clustDF <- clustersSummaryData(objCOTAN, clName = clName, clusters = clusters)
  rownames(clustDF) <- clustDF[[1L]]
  clustDF <- clustDF[rownames(scoreDFT), ]

  freq1 <- set_names(clustDF[["CellNumber"]], rownames(clustDF))

  freq2 <- set_names(paste0(clustDF[["CellPercentage"]], "%"),
                     rownames(clustDF))

  haList <- c(
    ha1 = rowAnnotation(
      cellNumber = anno_numeric(freq1,
                                rg = c(0L, max(freq1)),
                                bg_gp = gpar(fill = "orange", col = "black"),
                                labels_gp = gpar(fontsize = 10.0)),
      annotation_label = "CellNumber",
      annotation_name_rot = 0.0),

    ha2 = rowAnnotation(
      cellPerc = anno_text(freq2,
                           gp = gpar(fontsize = 10.0)),
      annotation_label = "CellPerc",
      annotation_name_rot = 0.0)
  )

  colorFunc <- colorRamp2(c(0.0, 1.0), c("lightblue", "red"))
  cellFunc <- function(j, i, x, y, width, height, fill) {
    grid.text(formatC(scoreDFT[i, j], digits = 1L, format = "f"),
              x, y, gp = gpar(fontsize = 9L))
  }

  heatmap <- Heatmap(scoreDFT,
                     rect_gp = gpar(col = "white", lwd = 1L),
                     cluster_rows = dend,
                     cluster_columns = FALSE,
                     col = colorFunc,
#                     row_dend_width = unit(8.0, "cm"),
#                     width = unit(10, "cm"),
#                     height = unit(10, "cm")
                     column_names_gp = gpar(fontsize = 11L),
                     row_names_gp = gpar(fontsize = 11L),
                     cell_fun = cellFunc,
                     right_annotation = haList,
                     left_annotation = hbList,
                     heatmap_legend_param = list(title = "score"))

  # TODO: avoid using the draw command here...
  heatmap <- draw(heatmap, annotation_legend_list = lgdList)

  return(list("heatmapPlot" = heatmap, "dataScore" = scoreDF))
}
