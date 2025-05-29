
#' Heatmap Plots
#'
#' @description These functions create heatmap `COEX` plots.
#'
#' @name HeatmapPlots
NULL

#' @details `singleHeatmapDF()` creates the `heatmap` `data.frame` of one
#'   `COTAN` object
#'
#' @param objCOTAN a `COTAN` object
#' @param genesLists A `list` of genes' `array`s. The first `array` defines the
#'   genes in the columns
#' @param sets A numeric array indicating which fields in the previous `list`
#'   should be used. Defaults to all fields
#' @param pValueThreshold The p-value threshold. Default is 0.01
#'
#' @returns `singleHeatmapDF()` returns a `data.frame`
#'
#' @importFrom tidyr pivot_longer
#'
#' @importFrom rlang is_empty
#' @importFrom rlang is_integer
#'
#' @importFrom assertthat assert_that
#'
#' @rdname HeatmapPlots
#'

singleHeatmapDF <- function(objCOTAN,
                            genesLists, sets,
                            pValueThreshold = 0.01) {
  assert_that(!is_empty(genesLists), !is_empty(names(genesLists)),
              msg = "genesLists must be a named list of genes arrays")

  assert_that(!is_empty(sets), is_integer(sets),
              max(sets) <= length(genesLists),
              msg = "sets must be positions in the genes list")

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  logThis(paste("Hangling COTAN object with condition:", cond), logLevel = 2L)

  allGenes <- getGenes(objCOTAN)

  colGenes <- genesLists[[1L]]
  genesLists <- genesLists[sets]

  askedGenes <- unique(unlist(genesLists))

  assert_that(any(colGenes %in% allGenes),
              msg = "none of primary markers are available in the COTAN object")

  #---------------------------------------------------------
  pValue <-
    calculatePValue(objCOTAN,
                    geneSubsetCol = allGenes[allGenes %in% colGenes],
                    geneSubsetRow = allGenes[allGenes %in% askedGenes])
  pValue <- as.data.frame(pValue)

  pValue <- setColumnInDF(pValue, colToSet = rownames(pValue), colName = "g2")

  dfValue <- pivot_longer(pValue, cols = seq_len(ncol(pValue)) - 1L,
                          names_to = "g1", values_to = "pValue")
  rm(pValue)

  #---------------------------------------------------------
  coex <- getGenesCoex(objCOTAN, genes = allGenes[allGenes %in% colGenes])
  coex <- coex[allGenes %in% askedGenes, , drop = FALSE]
  coex <- as.data.frame(coex)

  coex <- setColumnInDF(coex, colToSet = rownames(coex), colName = "g2")
  dfCoex <- pivot_longer(coex, cols = seq_len(ncol(coex)) - 1L,
                         names_to = "g1", values_to = "coex")
  rm(coex)

  #---------------------------------------------------------
  assert_that(identical(colnames(dfValue)[-3L], colnames(dfCoex)[-3L]))

  dfMerge <- merge(dfCoex, dfValue)

  dfMerge[["cond"]] <- cond
  dfMerge[["type"]] <- NA
  dfMerge[["absent"]] <- FALSE

  rm(dfCoex, dfValue)
  #---------------------------------------------------------

  dfRet <- data.frame()
  for (genesType in names(genesLists)) {
    logThis(paste("Handling genes type:", genesType), logLevel = 3L)
    rowGenes <- genesLists[[genesType]]
    for (g1 in colGenes) {
      genesToUse <- (dfMerge[["g2"]] %in% rowGenes) & (dfMerge[["g1"]] == g1)
      dfSlice <- dfMerge[genesToUse, , drop = FALSE]

      # control if the subset is smaller than the number of wanted genes
      numMissingRows <- length(rowGenes) - nrow(dfSlice)
      if (numMissingRows > 0L) {
        extraRows <- as.data.frame(matrix(nrow = numMissingRows,
                                          ncol = ncol(dfSlice)))
        colnames(extraRows) <- colnames(dfSlice)
        extraRows[, "g1"] <- g1
        extraRows[, "cond"] <- cond
        extraRows[, "absent"] <- TRUE
        extraRows[, "pValue"] <- 1.0
        extraRows[, "g2"] <- rowGenes[!(rowGenes %in% dfSlice[["g2"]])]
        dfSlice <- rbind(dfSlice, extraRows)
      }
      # now dfSlice has the right number of rows
      assert_that(nrow(dfSlice) == length(rowGenes))

      dfSlice[["type"]] <- genesType
      dfRet <- rbind(dfRet, dfSlice)
    }
  }

  dfRet[["fe_genes"]] <- (dfRet[["g2"]] %in% getFullyExpressedGenes(objCOTAN)) |
                         (dfRet[["g1"]] %in% getFullyExpressedGenes(objCOTAN))

  dfRet[dfRet[["pValue"]] > pValueThreshold, "coex"] <- 0.0

  return(dfRet)
}


#' @details `heatmapPlot()` creates the `heatmap` of one or more `COTAN` objects
#'
#' @param objCOTAN a `COTAN` object
#' @param genesLists A `list` of genes' `array`s. The first `array` defines the
#'   genes in the columns
#' @param sets A numeric array indicating which fields in the previous `list`
#'   should be used. Defaults to all fields
#' @param pValueThreshold The \eqn{p}-value threshold. Default is \eqn{0.01}
#' @param conditions An `array` of prefixes indicating the different files
#' @param dir The directory in which are all `COTAN` files (corresponding to the
#'   previous prefixes)
#'
#' @returns `heatmapPlot()` returns a `ggplot2` object
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
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateLambdaLinear(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 6L)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
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
#' hPlot <- heatmapPlot(objCOTAN, pValueThreshold = 0.05,
#'                      genesLists = groupMarkers, sets = 2L:3L)
#' plot(hPlot)
#'
#' @rdname HeatmapPlots
#'
heatmapPlot <- function(objCOTAN = NULL,
                        genesLists,
                        sets = NULL,
                        pValueThreshold = 0.01,
                        conditions = NULL,
                        dir = ".") {
  logThis("heatmap plot: START", logLevel = 2L)

  assert_that(is.null(objCOTAN) != is_empty(conditions),
              msg = paste("Please pass either a COTAN object",
                          "or an array of conditions"))

  if (is_empty(sets)) {
    sets <- seq_along(genesLists)
  }

  dfToPrint <- data.frame()
  if (!is.null(objCOTAN)) {
    dfToPrint <- singleHeatmapDF(objCOTAN, genesLists = genesLists,
                                 sets = sets, pValueThreshold = pValueThreshold)
  } else {
    for (cond in conditions) {
      logThis(paste0("Loading condition '", cond, "'"), logLevel = 3L)
      tempObj <- tryCatch({
          readRDS(file.path(dir, paste0(cond, ".cotan.RDS")))
        },
        error = function(err) {
          logThis(paste("While loading the COTAN object with file name",
                        paste0(cond, ".cotan.RDS"), "got an error:", err),
                  logLevel = 0L)
          return(NULL)
        })
      if (is.null(tempObj)) {
        next
      }
      dfTemp <- singleHeatmapDF(tempObj, genesLists = genesLists,
                                sets = sets, pValueThreshold = pValueThreshold)
      rm(tempObj)
      dfToPrint <- rbind(dfToPrint, dfTemp)
    }
  }

  logThis(paste("min coex:", min(dfToPrint[["coex"]], na.rm = TRUE),
                "max coex:", max(dfToPrint[["coex"]], na.rm = TRUE)),
          logLevel = 2L)

  heatmap <- ggplot(data = subset(dfToPrint, type %in% names(genesLists)[sets]),
                    aes(cond, factor(g2, levels = rev(levels(factor(g2)))))) +
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
#' @importFrom rlang is_empty
#'
#' @examples
#' ghPlot <- genesHeatmapPlot(objCOTAN, primaryMarkers = primaryMarkers,
#'                            secondaryMarkers = groupMarkers,
#'                            pValueThreshold = 0.05, symmetric = FALSE)
#' plot(ghPlot)
#'
#' @rdname HeatmapPlots
#'
genesHeatmapPlot <-
  function(objCOTAN, primaryMarkers,
           secondaryMarkers = vector(mode = "character"),
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
    if (isTRUE(symmetric)) {
      coex <- coex[, getGenes(objCOTAN) %in% rowGenes]
    }

    listRows <- list()
    for (m in unlist(secondaryMarkers)) {
      genes <- rownames(pValue[pValue[, m] < pValueThreshold, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0.0, ])]
      listRows[[m]] <- genes
    }

    clGenesRows <- data.frame()
    for (g in names(listRows)) {
      tmp <- data.frame("genes" = listRows[[g]],
                        "cl" = rep(g, length(listRows[[g]])))
      clGenesRows <- rbind(clGenesRows, tmp)
    }

    reorderIdxRow <- match(clGenesRows[["genes"]], rownames(coex))

    listCols <- list()
    for (m in primaryMarkers) {
      genes <- rownames(pValue[pValue[, m] < pValueThreshold, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0.0, ])]
      listCols[[m]] <- genes
    }

    if (isTRUE(symmetric)) {
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
#' ## plot(chPlot)
#'
#' @rdname HeatmapPlots
#'
cellsHeatmapPlot <- function(objCOTAN, cells = NULL, clusters = NULL) {
  coexMat <- as.matrix(getCellsCoex(objCOTAN))
  assert_that(!is_empty(coexMat), msg = "cells coex not found in the COTAN")

  # if clustering is needed
  if (!is_empty(clusters)) {
    # identifier for each cluster
    clustersTags <- levels(factor(clusters))

    # size of each cluster
    clustersList <- toClustersList(clusters)
    clustersSize <- lengths(clustersList)

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
