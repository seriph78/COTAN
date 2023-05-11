
#' @details `clustersMarkersHeatmapPlot()` returns the heatmap plot of a summary
#'   score for each *cluster* and each gene marker list in the given
#'   *clusterization*. It also returns the numerosity and percentage of each
#'   *cluster* on the right and a gene clusterization dendogram on the left (as
#'   returned by the function [geneSetEnrichment()]) that allows to estimate
#'   which markers groups are more or less expressed in each *cluster* so it is
#'   easier to derive the *clusters*' cell types.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the clusterization. If not given the last available
#'   clusterization will be used, as it is probably the most significant!
#' @param groupMarkers a named `list` with an element for each group of one or
#'   more marker genes for each group.
#' @param kCuts the number of estimated *cluster* (this defines the height for
#'   the tree cut and the associated colors)
#' @param conditionsList a list of `data.frames` coming from the
#'   [clustersSummaryPlot()] function
#'
#' @returns `clustersMarkersHeatmapPlot()` returns a list with:
#'  * "heatmapPlot" the complete heatmap plot
#'  * "dataScore" the `data.frame` with the score values
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
#' @importFrom tidyr `%>%`
#'
#' @importFrom grid unit
#' @importFrom grid gpar
#' @importFrom grid grid.text
#'
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom ComplexHeatmap anno_barplot
#' @importFrom ComplexHeatmap anno_numeric
#' @importFrom ComplexHeatmap anno_text
#' @importFrom ComplexHeatmap Legend
#' @importFrom ComplexHeatmap Heatmap
#'
#' @importFrom circlize colorRamp2
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
clustersMarkersHeatmapPlot <- function(objCOTAN, groupMarkers, clName = NULL,
                                       kCuts = 3, conditionsList = NULL) {
  clName <- getClusterizationName(objCOTAN, clName = clName)

  expressionCl <- clustersDeltaExpression(objCOTAN, clName = clName)
  scoreDF <- geneSetEnrichment(groupMarkers = groupMarkers,
                               clustersCoex = expressionCl)

  scoreDFT <- t(scoreDF[, 1L:(ncol(scoreDF) - 2L)])
  dend <- clustersTreePlot(objCOTAN, kCuts = kCuts)[["dend"]]

  dend <- set(dend = dend, "branches_lwd", 2)

  {
    numDigits <- floor(log10(nrow(scoreDFT))) + 1L
    rownames(scoreDFT) <- formatC(as.numeric(rownames(scoreDFT)),
                                  width = numDigits, flag = "0")
  }

  if (!is_empty(conditionsList)) {
    hbList <- list()
    lgdList <- list()
    allColors <- getColorsVector()

    for (condDF in conditionsList) {
      # a data frame coming from clustersSummaryPlot()
      assert_that(ncol(condDF) >= 3 &&
                    colnames(condDF)[[1L]] == "Cluster" &&
                    colnames(condDF)[[3L]] == "CellNumber",
                  msg = "Passed condition data.frame must have 3 columns")

      assert_that(all.equal(levels(ordered(condDF[[1L]])),
                            levels(ordered(clusters))),
                  msg = "Passed conditions do not agree with given clusters")

      condName <- colnames(condDF)[[2L]]
      cNames <- c("Cluster", condName, "CellNumber")

      if (is.numeric(condDF[[1L]])) {
        numDigits <- floor(log10(nrows(condDF))) + 1L
        condDF[[1L]] <-
          factor(formatC(condDF[[1L]], width = numDigits, flag = "0"))
      }

      condDF <- condDF[, cNames] %>%
        pivot_wider(names_from = condName,
                    values_from = cNames[[3L]])

      # Here we have a column per condition alternative plus 1
      nCols <- ncol(condDF)
      assert_that(nCols > 1L, msg = "Internal error: provided empty conditions")
      condDF[, 2:nCols] <- condDF[, 2:nCols] / rowSums(condDF[, 2:nCols])

      condDF <- as.data.frame(condDF)
      rownames(condDF) <- condDF[["Cluster"]]
      condDF <- condDF[rownames(scoreDFT), ]
      condDF[is.na(condDF)] <- 0L

      conds <- colnames(condDF)[2:nCols]

      colrs <- set_names(head(allColors, nCols - 1L), conds)
      allColors <- tail(allColors, length(allColors) - nCols + 1L)

      hb <- rowAnnotation(
        condName = anno_barplot(condDF[, 2:nCols],
                                width = unit(3.0, "cm"),
                                gp = gpar(fill = colrs, col = "black"),
                                align_to = "right",
                                labels_gp = gpar(fontsize = 12L)),
        annotation_name_rot = 0L)
      hbList <- append(hbList, hb)

      lgd <-  Legend(labels = conds, title = condName,
                     legend_gp = gpar(fill = colrs))
      lgdList <- append(lgdList, lgd)
    }
  } else {
    hb <- c()
    lgdList <- list()
  }

  clsInfo <- clustersSummaryPlot(objCOTAN, clName = clName)[["data"]]
  {
    numDigits <- floor(log10(nrow(clsInfo))) + 1L
    clsInfo[["Cluster"]] <-
      formatC(clsInfo[["Cluster"]], width = numDigits, flag = "0")
  }

  rownames(clsInfo) <- clsInfo[["Cluster"]]
  clsInfo <- clsInfo[rownames(scoreDFT), ]

  freq1 <- set_names(clsInfo[["CellNumber"]], rownames(clsInfo))

  freq2 <- set_names(paste0(clsInfo[["CellPercentage"]], "%"),
                     rownames(clsInfo))

  ha <- c(
    ha1 = rowAnnotation(
      cell.number = anno_numeric(freq1,
                                 bg_gp = gpar(fill = "orange", col = "black"),
                                 labels_gp = gpar(fontsize = 10)),
      annotation_name_rot = 0),

    ha2 = rowAnnotation(
      cell.perc = anno_text(freq2,
                            gp = gpar(fontsize = 10)),
      annotation_name_rot = 0)
  )

  colorFunc <- colorRamp2(c(0, 1), c("lightblue", "red"))
  cellFunc <- function(j, i, x, y, width, height, fill) {
    grid.text(formatC(scoreDFT[i, j], digits = 1L, format = "f"),
              x, y, gp = gpar(fontsize = 9L))
  }

  finalHeatmap <- Heatmap(scoreDFT, rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = FALSE,
                          col = colorFunc,
                          width = unit(28.0, "cm"),
                          row_dend_width = unit(8.0, "cm"),
                          #height = unit(6.0, "cm"),
                          column_names_gp = gpar(fontsize = 11L),
                          row_names_gp = gpar(fontsize = 11L),
                          cell_fun = cellFunc,
                          right_annotation = ha,
                          left_annotation = hb
             #, column_title = paste0("Gene marker set expression in ", sample)
                         )

  finalHeatmap <- draw(finalHeatmap, annotation_legend_list = lgdList)

  return(list("heatmapPlot" = finalHeatmap, "dataScore" = scoreDF))
}
