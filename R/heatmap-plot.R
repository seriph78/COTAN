
#' Heatmap Plots
#'
#' @description These functions create heatmap `COEX` plots.
#'
#' @name HeatmapPlots
NULL

#' @details `heatmapPlot()` creates the `heatmap` of one or more `COTAN` objects
#'
#' @param genesLists A `list` of genes' `array`s. The first `array` defines the
#'   genes in the columns
#' @param sets A numeric array indicating which fields in the previous `list`
#'   should be used
#' @param conditions An array of prefixes indicating the different files
#' @param pValueThreshold The p-value threshold. Default is 0.01
#' @param dir The directory in which are all `COTAN` files (corresponding to the
#'   previous prefixes)
#'
#' @returns `heatmapPlot()` returns a `ggplot2` object
#'
#' @importFrom Matrix forceSymmetric
#'
#' @importFrom tidyr pivot_longer
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_fill_gradient2
#'
#' @importFrom scales squish
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#'
#' ## Save the `COTAN` object to file
#' data_dir <- tempdir()
#' saveRDS(objCOTAN, file = file.path(data_dir, "test.dataset.cotan.RDS"))
#'
#' ## some genes
#' primaryMarkers <- c("g-000010", "g-000020", "g-000030")
#'
#' ## an example of named list of different gene set
#' groupMarkers <- list(G1 = primaryMarkers,
#'                      G2 = c("g-000300", "g-000330"),
#'                      G3 = c("g-000510", "g-000530", "g-000550",
#'                             "g-000570", "g-000590"))
#'
#' hPlot <- heatmapPlot(genesLists = groupMarkers, sets = c(2, 3),
#'                      pValueThreshold = 0.05, conditions = c("test.dataset"),
#'                      dir = paste0(data_dir, "/"))
#'
#' @rdname HeatmapPlots
#'
heatmapPlot <- function(genesLists, sets, conditions, dir,
                        pValueThreshold = 0.01) {
  logThis("heatmap plot: START", logLevel = 2L)

  colGenes <- genesLists[[1L]]
  allGenes <- unique(sort(unlist(genesLists[sets])))

  dfToPrint <- data.frame()
  for (cond in conditions) {
    logThis(paste0("Loading condition '", cond, "'"), logLevel = 3L)
    obj <- readRDS(file.path(dir, paste0(cond, ".cotan.RDS")))

    assert_that(any(colGenes %in% getGenes(obj)),
                msg = "primary markers all absent")

    pValue <- calculatePValue(obj, geneSubsetCol = colGenes,
                              geneSubsetRow = allGenes)
    pValue <- as.data.frame(pValue)

    pValue[["g2"]] <- as.vector(rownames(pValue))
    dfTempVal <- pivot_longer(pValue, cols = seq_along(colnames(pValue)) - 1L,
                              names_to = "g1", values_to = "pValue")

    #---------------------------------------------------------
    coex <- getGenesCoex(obj)
    coex <- coex[getGenes(obj) %in% allGenes, getGenes(obj) %in% colGenes]
    coex <- as.data.frame(coex)

    coex[["g2"]] <- as.vector(rownames(coex))
    dfTempCoex <- pivot_longer(coex, cols = seq_along(colnames(pValue)) - 1L,
                               names_to = "g1", values_to = "coex")

    dfTemp <- merge(dfTempCoex, dfTempVal)
    dfTemp[["time"]] <- cond
    dfTemp[["type"]] <- NA
    dfTemp[["absent"]] <- NA
    dfTemp2 <- data.frame()
    for (type in names(genesLists)[sets]) {
      for (g1 in colGenes) {
        tt <- dfTemp[dfTemp[["g2"]] %in% genesLists[[type]] &
                       dfTemp[["g1"]] == g1, ]
        # control if the subset is smaller than the number of wanted genes
        if (dim(tt)[[1L]] < length(genesLists[[type]])) {
          nRow <- length(genesLists[[type]]) - dim(tt)[[1L]]
          tRows <- as.data.frame(matrix(nrow = nRow, ncol = 7L))
          colnames(tRows) <- colnames(tt)
          tRows[, "g1"] <- g1
          tRows[, "time"] <- cond
          tRows[, "absent"] <- "yes"
          tRows[, "pValue"] <- 1.0
          tRows[, "g2"] <-
            genesLists[[type]][!genesLists[[type]] %in% tt[["g2"]]]
          tt <- rbind(tt, tRows)
        }
        tt[["type"]] <- type
        dfTemp2 <- rbind(dfTemp2, tt)
      }
      logThis(paste("Done type ", type), logLevel = 3L)
    }
    dfTemp <- dfTemp2
    dfTemp[["t_hk"]] <-
      ifelse((dfTemp[["g2"]] %in% getFullyExpressedGenes(obj)) |
             (dfTemp[["g1"]] %in% getFullyExpressedGenes(obj)),
             "hk", "n")
    dfTemp[dfTemp[["pValue"]] > pValueThreshold, "coex"] <- 0L
    dfToPrint <- rbind(dfToPrint, dfTemp)
  }
  logThis(paste("min coex:", min(dfToPrint[["coex"]], na.rm = TRUE),
                "max coex:", max(dfToPrint[["coex"]], na.rm = TRUE)),
          logLevel = 2L)

  heatmap <- ggplot(data = subset(dfToPrint, type %in% names(genesLists)[sets]),
                    aes(time, factor(g2, levels = rev(levels(factor(g2)))))) +
             geom_tile(aes(fill = coex), colour = "black", show.legend = TRUE) +
             facet_grid(type ~ g1, scales = "free", space = "free") +
             scale_fill_gradient2(low = "#E64B35FF", mid = "gray93",
                                  high = "#3C5488FF", midpoint = 0L,
                                  na.value = "grey80", space = "Lab",
                                  guide = "colourbar", aesthetics = "fill",
                                  oob = squish) +
             plotTheme("heatmap", textSize = 9L)

  return(heatmap)
}



#' @details `genesHeatmapPlot()` is used to plot an *heatmap* made using only
#'   some genes, as markers, and collecting all other genes correlated with
#'   these markers with a p-value smaller than the set threshold. Than all
#'   relations are plotted. Primary markers will be plotted as groups of rows.
#'   Markers list will be plotted as columns.
#'
#' @param objCOTAN a `COTAN` object
#' @param primaryMarkers A set of genes plotted as rows
#' @param secondaryMarkers A set of genes plotted as columns
#' @param pValueThreshold The p-value threshold. Default is 0.01
#' @param symmetric A Boolean: default `TRUE`. When `TRUE` the union of
#'   `primaryMarkers` and `secondaryMarkers` is used for both rows and column
#'   genes
#'
#' @returns `genesHeatmapPlot()` returns a `ggplot2` object
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap Legend
#' @importFrom ComplexHeatmap draw
#'
#' @importFrom grid gpar
#'
#' @importFrom stats quantile
#'
#' @importFrom circlize colorRamp2
#'
#' @importFrom Matrix forceSymmetric
#'
#' @importFrom rlang is_empty
#'
#' @examples
#' ghPlot <- genesHeatmapPlot(objCOTAN, primaryMarkers = primaryMarkers,
#'                            secondaryMarkers = groupMarkers,
#'                            pValueThreshold = 0.05, symmetric = FALSE)
#'
#' @rdname HeatmapPlots
#'
genesHeatmapPlot <-
  function(objCOTAN, primaryMarkers, secondaryMarkers = c(),
           pValueThreshold = 0.01, symmetric = TRUE) {
    if (isTRUE(symmetric)) {
      secondaryMarkers <- as.list(c(unlist(primaryMarkers),
                                    unlist(secondaryMarkers)))
    }

    if (is_empty(secondaryMarkers)) {
      secondaryMarkers <- as.list(primaryMarkers)
    } else {
      secondaryMarkers <- as.list(secondaryMarkers)
    }

    allMarkers <- unique(c(unlist(secondaryMarkers), unlist(primaryMarkers)))
    unexpressedGenes <- allMarkers[!allMarkers %in% getGenes(objCOTAN)]

    if (!is_empty(unexpressedGenes)) {
      warning("The markers ", toString(unexpressedGenes),
              " are not present in the COTAN object!")
    }

    pValue <- calculatePValue(objCOTAN)[, allMarkers]

    pvalFloored <- apply(pValue, 1.0, FUN = min)
    rowGenes <- names(pvalFloored)[pvalFloored < pValueThreshold]

    rowGenes <- unique(c(allMarkers, rowGenes))
    pValue <- as.data.frame(pValue)

    coex <- getGenesCoex(objCOTAN)[getGenes(objCOTAN) %in% rowGenes, ]
    if (symmetric == TRUE) {
      coex <- coex[, getGenes(objCOTAN) %in% rowGenes]
    }

    listRows <- c()
    for (m in unlist(secondaryMarkers)) {
      genes <- rownames(pValue[pValue[, m] < pValueThreshold, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0.0, ])]
      listRows[[m]] <- genes
    }

    clGenesRows <- c()
    for (g in names(listRows)) {
      tmp <- data.frame("genes" = listRows[[g]],
                        "cl" = rep(g, length(listRows[[g]])))
      clGenesRows <- rbind(clGenesRows, tmp)
    }

    reorderIdxRow <- match(clGenesRows[["genes"]], rownames(coex))

    listCols <- c()
    for (m in primaryMarkers) {
      genes <- rownames(pValue[pValue[, m] < pValueThreshold, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0.0, ])]
      listCols[[m]] <- genes
    }

    if (symmetric == TRUE) {
      clGenesCols <- clGenesRows
    } else {
      clGenesCols <- data.frame()
      for (g in names(listCols)) {
        tmp <- data.frame("genes" = listCols[[g]],
                          "cl" = rep(g, length(listCols[[g]])))
        clGenesCols <- rbind(clGenesCols, tmp)
      }
    }

    clGenesCols <- clGenesCols[clGenesCols[["genes"]] %in% colnames(coex), ]

    reorderIdxCol <- match(clGenesCols[["genes"]], colnames(coex))

    coex <- coex[reorderIdxRow, reorderIdxCol]

    colFun <- colorRamp2(
      c(round(quantile(as.matrix(coex), probs = 0.001), digits = 3L),
        0.0,
        round(quantile(as.matrix(coex), probs = 0.999), digits = 3L)),
      c("#E64B35FF", "gray93", "#3C5488FF")
    )

    # The next line is to set the columns and raws order
    # need to be implemented
    part1 <- Heatmap(as.matrix(coex),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_split = clGenesRows[["cl"]],
      column_split = clGenesCols[["cl"]],
      col = colFun,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title_gp = gpar(
        fill = "#8491B44C", font = 3L,
        col = "#3C5488FF"
      ),
      row_title_gp = gpar(fill = "#8491B44C", font = 3L, col = "#3C5488FF")
    )
    lgd <- Legend(
      col_fun = colFun, title = "coex", grid_width = unit(0.3, "cm"),
      direction = "horizontal", title_position = "topcenter",
      title_gp = gpar(fontsize = 10L, fontface = "bold", col = "#3C5488FF"),
      labels_gp = gpar(col = "#3C5488FF", font = 3L)
    )
    draw(part1,
      show_heatmap_legend = FALSE,
      annotation_legend_list = lgd, annotation_legend_side = "bottom"
    )
  }


#' @details `cellsHeatmapPlot()` creates the heatmap plot of the cells' `COEX`
#'   matrix
#'
#' @param objCOTAN a `COTAN` object
#' @param cells Which cells to plot (all if no argument is given)
#' @param clusters Use this clusterization to select/reorder the cells to plot
#'
#' @returns `cellsHeatmapPlot()` returns the cells' `COEX` *heatmap* plot
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom ComplexHeatmap Heatmap
#'
#' @examples
#' clusters <- c(rep_len("1", getNumCells(objCOTAN)/2),
#'               rep_len("2", getNumCells(objCOTAN)/2))
#' names(clusters) <- getCells(objCOTAN)
#'
#' chPlot <- cellsHeatmapPlot(objCOTAN, clusters = clusters)
#'
#' @rdname HeatmapPlots
#'
cellsHeatmapPlot <- function(objCOTAN, cells = NULL, clusters = NULL) {
  coexMat <- as.matrix(getCellsCoex(objCOTAN))
  assert_that(!is_empty(coexMat), msg = "cells coex not found in the COTAN")

  # if clustering is needed
  if (!is_empty(clusters)) {
    # identifier for each cluster
    clustersTags <- unique(clusters)

    # size of each cluster
    clustersList <- toClustersList(clusters)
    clustersSize <- vapply(clustersList, length, integer(1))

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
