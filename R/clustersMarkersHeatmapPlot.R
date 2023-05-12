
#' @details `clustersMarkersHeatmapPlot()` returns the heatmap plot of a summary
#'   score for each *cluster* and each gene marker list in the given
#'   *clusterization*. It also returns the numerosity and percentage of each
#'   *cluster* on the right and a gene clusterization dendogram on the left (as
#'   returned by the function [geneSetEnrichment()]) that allows to estimate
#'   which markers groups are more or less expressed in each *cluster* so it is
#'   easier to derive the *clusters*' cell types.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
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
    stop("Uasge of condition list is not supported yet")
    # # a data frame coming from clustersSummaryPlot()
    # cond1 <- conditionsList[[1L]]
    #
    # if (is.numeric(cond1[["Cluster"]])) {
    #   numDigits <- floor(log10(length(cond1[["Cluster"]]))) + 1L
    #   cond1[["Cluster"]] <-
    #     formatC(cond1[["Cluster"]], width = numDigits, flag = "0")
    # }
    #
    # cond1 <- cond1[, c("Cluster", "cond1" ,"CellNumber")] %>%
    #   pivot_wider(names_from = cond1, values_from = CellNumber)
    #
    # cond1[, 2L:3L] <- cond1[, 2L:3L] / rowSums(cond1[, 2L:3L])
    # cond1 <- as.data.frame(cond1)
    # rownames(cond1) <- cond1[["Cluster"]]
    # cond1 <- cond1[rownames(scoreDFT), ]
    # cond1[is.na(cond1)] <- 0
    #
    # # TODO: here instead of F and M we need a flexible number of conditions
    # cond1Col <- c("F" = "deeppink1", "M" = "darkturquoise")
    # cond2Col <- c("F"="red4", "R"="sienna3", "H"="seagreen")
    # hb <- list(
    #   hb1 = rowAnnotation(
    #     cond1 = anno_barplot(cond1[, 2L:3L],
    #                          width = unit(3.0, "cm"),
    #                          gp = gpar(fill = cond1Col, col = "black"),
    #                          align_to = "right",
    #                          labels_gp = gpar(fontsize = 12L)),
    #     annotation_name_rot = 0L),
    #   hb2 = rowAnnotation(condition = anno_barplot(cond2[, 2L:3L],
    #                                     width = unit(3.0, "cm"),
    #                                     gp = gpar(fill = cond2Col,
    #                                               col = "black"),
    #                                     align_to = "right",
    #                                     labels_gp = gpar(fontsize = 12L)),
    #                       annotation_name_rot = 0L)
    # )
    # lgdList <- list(
    #   lb1 = Legend(labels = c("Female", "Male"), title = "cond1",
    #                legend_gp = gpar(fill = cond1Col))#,
    #   lb2 = Legend(labels = c("Flare", "Remission", "Healthy"), title = "cond2",
    #                legend_gp = gpar(fill = cond2col))
    # )
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

  freq <- set_names(clsInfo[["CellNumber"]], rownames(clsInfo))

  freq2 <- set_names(paste0(clsInfo[["CellPercentage"]], "%"),
                     rownames(clsInfo))

  ha <- c(
    ha1 = rowAnnotation(cell.number = anno_numeric(freq,
                                                   bg_gp = gpar(fill = "orange",
                                                                col = "black")
                                                #labels_gp = gpar(fontsize = 10)
                                                  ),
                        annotation_name_rot = 0),

    ha2 = rowAnnotation(cell.perc = anno_text(freq2,
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
