
#' @details `clustersMarkersHeatmapPlot()` returns the heatmap plot of a summary
#'   score for each *cluster* and each gene marker in the given
#'   *clusterization*. It also returns the size and percentage of each
#'   *cluster* on the right and a *clusterization* `dendogram` on the left, as
#'   returned by the function [clustersTreePlot()]. The heatmap cells' colors
#'   express the **DEA**, that is whether a gene is enriched or depleted in the
#'   cluster, while the stars are aligned to the corresponding adjusted
#'   \eqn{p-}value: `***` for \eqn{p < 0.1\%}, `**` for \eqn{p < 1\%}, `*` for
#'   \eqn{p < 5\%}, `.` for \eqn{p < 10\%}
#'
#' @param objCOTAN a `COTAN` object
#' @param groupMarkers an optional named `list` with an element for each group
#'   comprised of one or more marker genes
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName` that, in such a case, will only indicate
#'   the relevant column name in the returned `data.frame`
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#' @param kCuts the number of estimated *cluster* (this defines the height for
#'   the tree cut and the associated colors)
#' @param adjustmentMethod *p-value* multi-test adjustment method. Defaults to
#'   `"bonferroni"`; use `"none"` for no adjustment
#' @param condNameList a `list` of *conditions*' names to be used for additional
#'   columns in the final plot. When none are given no new columns will be added
#'   using data extracted via the function [clustersSummaryData()]
#' @param conditionsList a `list` of *conditions* to use. If given they will
#'   take precedence on the ones indicated by `condNameList`
#'
#' @returns `clustersMarkersHeatmapPlot()` returns a list with:
#'  * `"heatmapPlot"` the complete heatmap plot
#'  * `"dataScore"` the `data.frame` with the score values
#'  * `"pValueDF"`  the `data.frame` with the corresponding adjusted
#'   \eqn{p-}values
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom stats symnum
#'
#' @importFrom dendextend set
#'
#' @importFrom tidyr pivot_wider
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
clustersMarkersHeatmapPlot <- function(objCOTAN, groupMarkers = list(),
                                       clName = "", clusters = NULL,
                                       coexDF = NULL, kCuts = 3L,
                                       adjustmentMethod = "bonferroni",
                                       condNameList = NULL,
                                       conditionsList = NULL) {
  assert_that(is_empty(conditionsList) ||
                length(conditionsList) == length(condNameList),
              msg = "Explicitly given conditions must have corresponding names")

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  if (is_empty(coexDF)) {
    if (clName %in% getClusterizations(objCOTAN)) {
      coexDF <- getClusterizationData(objCOTAN, clName = clName)[["coex"]]
    }
    if (is_empty(coexDF)) {
      coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)
    }
  }

  genesGroupMap <- NULL
  if (!is_empty(groupMarkers)) {
    genesGroupMap <-
      unlist(lapply(
        names(groupMarkers),
        \(group) {
          genes <- groupMarkers[[group]]
          set_names(rep_len(group, length(genes)), genes)
        }
      ))

    genesGroupMap <- factor(genesGroupMap, levels = names(groupMarkers))
  }

  scoreDF <- coexDF[names(genesGroupMap), , drop = FALSE]
  rownames(scoreDF) <- names(genesGroupMap)

  pValueDF <- pValueFromDEA(coexDF = coexDF, numCells = getNumCells(objCOTAN),
                            adjustmentMethod = adjustmentMethod)
  pValueDF <- pValueDF[names(genesGroupMap), , drop = FALSE]
  rownames(pValueDF) <- names(genesGroupMap)


  dend <- clustersTreePlot(objCOTAN, kCuts = kCuts, clName = clName)[["dend"]]
  dend <- set(dend, "branches_lwd", 2L)

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

      condDF <- pivot_wider(condDF[, cNames],
                            names_from = condName,
                            values_from = cNames[[3L]])

      # Here we have a column per condition alternative plus 1
      nCols <- ncol(condDF)
      assert_that(nCols > 1L, msg = "Internal error: provided empty conditions")
      condDF[, 2L:nCols] <- condDF[, 2L:nCols] / rowSums(condDF[, 2L:nCols])

      condDF <- as.data.frame(condDF)
      rownames(condDF) <- condDF[[1L]]
      condDF <- condDF[colnames(scoreDF), ]
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
  clustDF <- clustDF[colnames(scoreDF), ]

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

  # Define the color function for the heatmap
  colExtr <- suppressWarnings(max(-min(scoreDF, na.rm = TRUE),
                                   max(scoreDF, na.rm = TRUE)))
  colorFunc <- colorRamp2(c(-colExtr, 0.0, colExtr),
                          c("#E64B35FF", "gray93", "#3C5488FF"))

  # nolint start: spaces_inside_linter
  cutPoints <- c(0.0,        0.001,       0.01,      0.05,      0.1,     1.0)
  marks     <- c(     "***",        "**",       "*",       ".",      " "    )
  # nolint end

  # Define the function to annotate each cell with its significance
  cellFunc <- function(j, i, x, y, width, height, fill) {
    pValue <- pValueDF[j, i]

    # Use symnum() to get the significance mark
    mark <- symnum(pValue, corr = FALSE, na = FALSE,
                   cutpoints = cutPoints, symbols = marks)

    # Display the mark in the cell
    grid.text(mark, x, y, gp = gpar(fontsize = 12L, fontface = "bold"))
  }

  heatMapPl <- Heatmap(
    t(scoreDF),
    name = "Score",
    rect_gp = gpar(col = "white", lwd = 1L),
    cluster_rows = dend,
    row_dend_width = unit(2.0, "cm"),
    cluster_columns = FALSE,  # Columns are not clustered but grouped
    col = colorFunc,
    column_split = genesGroupMap,      # Group columns by gene groups
    column_gap = unit(1.0, "mm"),        # Gap between gene groups
    column_names_gp = gpar(fontsize = 10L, angle = 45.0, hjust = 1.0),
    row_names_gp = gpar(fontsize = 10L),
    cell_fun = cellFunc,
    top_annotation = NULL,
    right_annotation = haList,
    left_annotation = hbList,
    heatmap_legend_param = list(title = "Score")
  )

  # If there are additional legends, add them
  if (!is.null(lgdList)) {
    heatMapPl <- heatMapPl + lgdList
  }

  return(list("heatmapPlot" = heatMapPl,
              "dataScore" = scoreDF,
              "pValues" = pValueDF))
}
