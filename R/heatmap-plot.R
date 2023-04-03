
#' heatmapPlot
#'
#' @description This function creates the `heatmap` of one or more `COTAN`
#'   objects.
#'
#' @param genesLists A `list` of genes' `array`s. The first `array` defines the
#'   genes in the columns.
#' @param sets A numeric array indicating which fields in the
#'   previous `list` should be used.
#' @param conditions An array of prefixes indicating the different files.
#' @param pValueThreshold The p-value threshold. Default is 0.05.
#' @param dir The directory in which are all `COTAN` files (corresponding to the
#'   previous prefixes)
#'
#' @returns a `ggplot` object
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
#' @export
#'
#' @examples
#' data("ERCCraw")
#' rownames(ERCCraw) = ERCCraw$V1
#' ERCCraw = ERCCraw[,2:ncol(ERCCraw)]
#' objCOTAN <- COTAN(raw = ERCCraw)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)
#' data_dir <- tempdir()
#' saveRDS(objCOTAN, file = file.path(data_dir, "ERCC.cotan.RDS"))
#' # some genes
#' primary.markers <- c("ERCC-00154", "ERCC-00156", "ERCC-00164")
#' # a example of named list of different gene set
#' gene.sets.list <- list(
#'   "primary.markers" = primary.markers,
#'   "2.R" = c("ERCC-00170", "ERCC-00158"),
#'   "3.S" = c("ERCC-00160", "ERCC-00162")
#' )
#' heatmapPlot(
#'   pValueThreshold = 0.05, genesLists = gene.sets.list,
#'   sets = c(2, 3), conditions = c("ERCC"), dir = paste0(data_dir, "/")
#' )
#'
#' @rdname heatmapPlot
#'
heatmapPlot <- function(genesLists, sets, conditions,
                        pValueThreshold = 0.05, dir) {
  #time <- g2 <- NULL
  logThis("heatmap plot: START", logLevel = 2)

  colGenes <- genesLists[[1]]
  allGenes <- unique(sort(unlist(genesLists[sets])))

  df.to.print <- data.frame()
  for (cond in conditions) {
    logThis(paste0("Loading condition '", cond, "'"), logLevel = 3)
    obj <- readRDS(file.path(dir, paste0(cond, ".cotan.RDS")))

    stopifnot("primary markers all absent" <- any(colGenes %in% getGenes(obj)))

    pValue <- calculatePValue(obj, geneSubsetCol = colGenes, geneSubsetRow = allGenes)
    pValue <- as.data.frame(pValue)

    pValue$g2 <- as.vector(rownames(pValue))
    df.temp.pval <- pivot_longer(pValue, cols = seq_along(colnames(pValue)) - 1,
                                 names_to = "g1", values_to = "pValue")

    #---------------------------------------------------------
    coex <- getGenesCoex(obj)
    coex <- coex[getGenes(obj) %in% allGenes, getGenes(obj) %in% colGenes]
    coex <- as.data.frame(coex)

    coex$g2 <- as.vector(rownames(coex))
    df.temp.coex <- pivot_longer(coex, cols = seq_along(colnames(pValue)) - 1,
                                 names_to = "g1", values_to = "coex")

    df.temp <- merge(df.temp.coex, df.temp.pval)
    df.temp$time <- cond
    df.temp$type <- NA
    df.temp$absent <- NA
    df.temp2 <- data.frame()
    for (type in names(genesLists)[sets]) {
      for (g1 in colGenes) {
        tt <- df.temp[df.temp$g2 %in% genesLists[[type]] & df.temp$g1 == g1, ]
        # control if the subset is smaller than the number of wanted genes
        if (dim(tt)[1] < length(genesLists[[type]])) {
          n.row <- length(genesLists[[type]]) - dim(tt)[1]
          t.rows <- as.data.frame(matrix(nrow = n.row, ncol = 7))
          colnames(t.rows) <- colnames(tt)
          t.rows[, "g1"] <- g1
          t.rows[, "time"] <- cond
          t.rows[, "absent"] <- "yes"
          t.rows[, "pValue"] <- 1
          t.rows[, "g2"] <- genesLists[[type]][!genesLists[[type]] %in% tt$g2]
          tt <- rbind(tt, t.rows)
        }
        tt$type <- type
        df.temp2 <- rbind(df.temp2, tt)
      }
      print(type)
    }
    df.temp <- df.temp2
    df.temp$t_hk <- ifelse((df.temp$g2 %in% getHousekeepingGenes(obj)) |
                             (df.temp$g1 %in% getHousekeepingGenes(obj)), "hk", "n")
    df.temp[df.temp$pValue > pValueThreshold, ]$coex <- 0
    df.to.print <- rbind(df.to.print, df.temp)
  }
  logThis(paste("min coex:", min(df.to.print$coex, na.rm = TRUE),
                "max coex",  max(df.to.print$coex, na.rm = TRUE)), logLevel = 2)

  heatmap <- ggplot(data = subset(df.to.print, type %in% names(genesLists)[sets]),
                    aes(time, factor(g2, levels = rev(levels(factor(g2)))))) +
             geom_tile(aes(fill = coex), colour = "black", show.legend = TRUE) +
             facet_grid(type ~ g1, scales = "free", space = "free") +
             scale_fill_gradient2(low = "#E64B35FF", mid = "gray93", high = "#3C5488FF",
                                  midpoint = 0,na.value = "grey80", space = "Lab",
                                  guide = "colourbar", aesthetics = "fill", oob = squish) +
             plotTheme("heatmap", textSize = 9)

  return(heatmap)
}



#' genesHeatmapPlot
#'
#' @description This function is used to plot an heatmap made using only some
#'   genes, as markers, and collecting all other genes correlated with these
#'   markers with a p-value smaller than the set threshold. Than all relations
#'   are plotted. Primary markers will be plotted as groups of rows. Markers
#'   list will be plotted as columns.
#'
#' @param objCOTAN A `COTAN` object.
#' @param primaryMarkers A set of genes plotted as rows.
#' @param secondaryMarkers A set of genes plotted as columns.
#' @param pValue The p-value threshold
#' @param symmetric A Boolean: default `TRUE`. When `TRUE` the union of
#'   `primaryMarkers` and `secondaryMarkers` is used for both rows and column
#'   genes
#'
#' @returns A ggplot2 object
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
#' @rdname genesHeatmapPlot
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)
#' # some genes
#' primaryMarkers <- c("g-000010", "g-000020", "g-000030")
#' # a example of named list of different gene set
#' groupMarkers <- list(G1 = primaryMarkers,
#'                      G2 = c("g-000300", "g-000330"),
#'                      G3 = c("g-000510", "g-000530", "g-000550", "g-000570", "g-000590"))
#' genesHeatmapPlot(objCOTAN,
#'                  primaryMarkers = primary.markers,
#'                  secondaryMarkers = groupMarkers,
#'                  pValue = 0.05, symmetric = FALSE)
genesHeatmapPlot <-
  function(objCOTAN, primaryMarkers, secondaryMarkers = c(),
           pValue = 0.001, symmetric = TRUE) {
    if (isTRUE(symmetric)) {
      secondaryMarkers <- as.list(c(unlist(primaryMarkers), unlist(secondaryMarkers)))
    }

    if (is_empty(secondaryMarkers)) {
      secondaryMarkers <- as.list(primaryMarkers)
    } else {
      secondaryMarkers <- as.list(secondaryMarkers)
    }

    allMarkers <- unique(c(unlist(secondaryMarkers), unlist(primaryMarkers)))
    unexpressedGenes <- allMarkers[!allMarkers %in% getGenes(objCOTAN)]

    if (!is_empty(unexpressedGenes)) {
      warning(paste("The markers", paste(unexpressedGenes, collapse = ", "),
                    "are not present in the COTAN object!"))
    }

    pval <- calculatePValue(objCOTAN)[, allMarkers]

    pvalFloored <- apply(pval, 1, FUN = min)
    rowGenes <- names(pvalFloored)[pvalFloored < pValue]

    rowGenes <- unique(c(allMarkers, rowGenes))
    pval <- as.data.frame(pval)

    coex <- getGenesCoex(objCOTAN)[getGenes(objCOTAN) %in% rowGenes, ]
    if (symmetric == TRUE) {
      coex <- coex[ , getGenes(objCOTAN) %in% rowGenes]
    }

    listRows <- c()
    for (m in unlist(secondaryMarkers)) {
      genes <- rownames(pval[pval[, m] < pValue, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0, ])]
      listRows[[m]] <- genes
    }

    clGenesRows <- c()
    for (g in names(listRows)) {
      tmp <- data.frame("genes" = listRows[[g]], "cl" = rep(g, length(listRows[[g]])))
      clGenesRows <- rbind(clGenesRows, tmp)
    }

    reorder_idx_row <- match(clGenesRows[["genes"]], rownames(coex))

    listCols <- c()
    for (m in primaryMarkers) {
      genes <- rownames(pval[pval[, m] < pValue, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0, ])]
      listCols[[m]] <- genes
    }

    if (symmetric == TRUE) {
      clGenesCols <- clGenesRows
    } else {
      clGenesCols <- data.frame()
      for (g in names(listCols)) {
        tmp <- data.frame("genes" = listCols[[g]], "cl" = rep(g, length(listCols[[g]])))
        clGenesCols <- rbind(clGenesCols, tmp)
      }
    }

    clGenesCols <- clGenesCols[clGenesCols[["genes"]] %in% colnames(coex), ]

    reorder_idx_col <- match(clGenesCols[["genes"]], colnames(coex))

    coex <- coex[reorder_idx_row, reorder_idx_col]

    col_fun <- colorRamp2(
      c(
        round(quantile(as.matrix(coex), probs = 0.001),
          digits = 3
        ), 0,
        round(quantile(as.matrix(coex), probs = 0.999),
          digits = 3
        )
      ),
      c("#E64B35FF", "gray93", "#3C5488FF")
    )

    # The next line is to set the columns and raws order
    # need to be implemented
    part1 <- Heatmap(as.matrix(coex),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_split = clGenesRows[["cl"]],
      column_split = clGenesCols[["cl"]],
      col = col_fun,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title_gp = gpar(
        fill = "#8491B44C", font = 3,
        col = "#3C5488FF"
      ),
      row_title_gp = gpar(fill = "#8491B44C", font = 3, col = "#3C5488FF")
    )
    lgd <- Legend(
      col_fun = col_fun, title = "coex", grid_width =
        unit(0.3, "cm"),
      direction = "horizontal", title_position = "topcenter",
      title_gp = gpar(fontsize = 10, fontface = "bold", col = "#3C5488FF"),
      labels_gp = gpar(col = "#3C5488FF", font = 3)
    )
    draw(part1,
      show_heatmap_legend = FALSE,
      annotation_legend_list = lgd, annotation_legend_side = "bottom"
    )
  }


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
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
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

