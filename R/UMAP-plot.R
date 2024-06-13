
#' @details `UMAPPlot()` plots the given `data.frame` containing genes
#'   information related to clusters after applying the [umap::umap()]
#'   transformation
#'
#' @param df The `data.frame` to plot. It must have a row names containing the
#'   given elements
#' @param clusters The **clusterization**. Must be a named `array` aligned to
#'   the rows in the `data.frame`.
#' @param elements a named `list` of elements to label. Each array in the list
#'   will be shown with a different color
#' @param title a string giving the plot title. Will default to UMAP Plot if not
#'   specified
#' @param colors an `array` of colors to use in the plot. If not sufficient
#'   colors are given it will complete the list using colors from
#'   [getColorsVector()]
#' @param numNeighbors Overrides the `n_neighbors` value from [umap.defaults]
#' @param minPointsDist Overrides the `min_dist` value from [umap.defaults]
#'
#'
#' @returns `UMAPPlot()` returns a `ggplot2` object
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggrepel geom_label_repel
#'
#' @importFrom stats quantile
#'
#' @importFrom umap umap
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
UMAPPlot <- function(df, clusters = NULL, elements = NULL, title = "",
                     colors = NULL, numNeighbors = 0L, minPointsDist = NA) {
  logThis("UMAP plot", logLevel = 2L)

  assert_that(is_empty(clusters) || length(clusters) == nrow(df),
              msg = paste("UMAPPlot - clusters vector must have size equal to",
                          "the number of rows in the data.frame"))

  # empty title
  if (isEmptyName(title)) {
    title <- "UMAP Plot"
  }

  if (!is_empty(clusters)) {
    clusters <- factor(clusters)
  }

  entryType <- rep_len("none", nrow(df))

  # assign a different color to each list of elements
  for (nm in names(elements)) {
    selec <- rownames(df) %in% elements[[nm]]
    if (any(selec)) {
      entryType[selec] <- nm
    } else {
      logThis(paste("UMAPPlot - none of the elements in group", nm,
                    "is present: will be ignored"), logLevel = 1L)
    }
  }

  labelled <- entryType != "none"

  # assign a different color to each cluster
  for (cl in levels(clusters)) {
    selec <- !labelled & clusters == cl
    min
    if (any(selec)) {
      entryType[selec] <- cl
    } else {
      logThis(paste("UMAPPlot - none of the elements of the cluster", cl,
                    "is present unlabelled: cluster will be ignored"),
              logLevel = 1L)
    }
    rm(selec)
  }

  clustered <- !labelled & entryType != "none"

  umapConfig <- umap.defaults
  if (numNeighbors != 0L) {
    umapConfig$n_neighbors <- numNeighbors
  }
  if (!is.na(minPointsDist)) {
    umapConfig$min_dist <- minPointsDist
  }
  umapConfig$verbose <- TRUE

  umap <- umap(df, config = umapConfig)

  plotDF <- data.frame(x = umap[["layout"]][, 1L],
                       y = umap[["layout"]][, 2L],
                       types = entryType)

  # add the centroids to the data.frame
  centroids <- rep(FALSE, times = nrow(plotDF))

  if (!is_empty(clusters)) {
    clList <- toClustersList(clusters)

    numericDF <- plotDF[, -3L]
    for (clName in names(clList)) {
      subsetDF <- as.matrix(numericDF[clList[[clName]], , drop = FALSE])
      plotDF <- rbind(plotDF, c(colMeans(subsetDF), clName))
      rownames(plotDF)[nrow(plotDF)] <- clName
    }

    numCl <- length(clList)
    centroids <- c(centroids, rep(TRUE, numCl))
    entryType <- c(as.character(entryType), rep("centroid", numCl))
    labelled <- c(labelled,  rep(FALSE, numCl))
    clustered <- c(clustered,  rep(FALSE, numCl))
    rm(numericDF, subsetDF, clList, clName, numCl)
  }

  # Ensure x and y columns are numeric
  plotDF[["x"]] <- as.numeric(plotDF[["x"]])
  plotDF[["y"]] <- as.numeric(plotDF[["y"]])
  plotDF[["types"]] <- factor(plotDF[["types"]])

  generic <- (!labelled & !clustered & !centroids)

  allTypes  <- c(unique(setdiff(entryType, c("none", "centroid"))),
                        "none", "centroid")
  myColours <- colors
  if (length(colors) < length(allTypes)) {
    if (!is_empty(colors)) {
      warning("UMAPPlot - not enough colors passed in")
    }
    numMissing <- length(allTypes) - length(myColours)
    myColours <- c(myColours, getColorsVector(numMissing))
    names(myColours) <- allTypes
    if (numMissing > 1L){
      myColours[["none"]] <- "#8491B4B2"
    }
    myColours[["centroid"]] <- "#000000"
  } else {
    myColours <- myColours[1L:length(allTypes)]
    names(myColours) <- allTypes
  }

  assert_that(setequal(c(entryType, "none", "centroid"), names(myColours)))

  pointSize <- min(max(1.0, 200000.0/dim(plotDF)[1L]), 3.0)

  plot <- ggplot() +
    geom_point(data = plotDF[generic, , drop = FALSE],
               aes(x, y, colour = types),
               size = pointSize, alpha = 0.3) +
    geom_point(data = plotDF[clustered, , drop = FALSE],
               aes(x, y, colour = types),
               size = pointSize, alpha = 0.5) +
    geom_point(data = plotDF[labelled, , drop = FALSE],
               aes(x, y, colour = types),
               size = 1.2 * pointSize, alpha = 0.8) +
    scale_color_manual("Status", values = myColours) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    ggtitle(title) +
    plotTheme("UMAP", textSize = 10L)

  plot <- plot +
    geom_label_repel(data = plotDF[labelled, , drop = FALSE],
                     aes(x, y, fill = types,
                         label = rownames(plotDF)[labelled]),
                     label.size = NA, show.legend = FALSE, force = 2.0,
                     box.padding = 0.25,
                     max.overlaps = 40L, alpha = 0.8,
                     direction = "both", na.rm = TRUE, seed = 1234L) +
    geom_text(data = plotDF[centroids, , drop = FALSE],
              aes(x, y, colour = "centroid"),
              label = rownames(plotDF)[centroids],
              show.legend = FALSE, alpha = 0.8,
              fontface = "bold", size = 1.5 * pointSize) +
    scale_fill_manual("Status", values = myColours)

  return(plot)
}



#' @details `cellsUMAPPlot()` returns a `ggplot2` plot where the given
#'   *clusters* are placed on the base of their relative distance. Also if
#'   needed calculates and stores the `DEA` of the relevant *clusterization*.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be returned, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#' @param method selects the method to use to create the `data.frame` to pass to
#'   the [UMAPPlot()]. The following methods calculate, for each cell, a
#'   statistic for each gene based on available data/model. The following
#'   methods are supported:
#'   * `"Likelyhood"` uses the likelyhood of observed presence/absence of each
#'   gene. The default method
#'   * `"LogNorm"` uses the log-normalized counts
#' @param geneSel Decides whether and how to perform gene-selection. It can be a
#'   string indicating one of the following methods or a straight list of genes.
#'   The following methods are supported:
#'   * `"HGDI"` Will pick-up the genes with highest **GDI**. The default method.
#'     Since it requires an available `COEX` matrix it will fall-back to `"HVG"`
#'     when not available
#'   * `"HVG"` Will pick-up the genes with the highest variability
#'   * `"None"` Will use all genes
#' @param numGenes the number of genes to select using the above method. Will be
#'   ignored when no selection have been asked or when an explicit list of genes
#'   was passed in
#' @param colors an `array` of colors to use in the plot. If not sufficient
#'   colors are given it will complete the list using colors from
#'   [getColorsVector()]
#'
#' @returns `cellsUMAPPlot()` returns a list with 2 objects:
#'  * `"plot"` a `ggplot2` object representing the `umap` plot
#'  * `"cellsPCA"` the `data.frame` PCA used to create the plot
#'
#' @importFrom rlang is_empty
#'
#' @importFrom stringr str_equal
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom PCAtools pca
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
cellsUMAPPlot <- function(objCOTAN,
                          clName = "",
                          clusters = NULL,
                          method = "Likelyhood",
                          geneSel = NULL,
                          numGenes = 2000L,
                          colors = NULL) {
  # pick last if no name was given
  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  if (is.null(clusters)) {
    logThis("No clusters given: using all dataset as one")
    clName <- "AllCells"
    clusters <- factor(set_names(rep("1", getNumCells(objCOTAN)),
                                 getCells(objCOTAN)))
  }

  assert_that(inherits(clusters, "factor"),
              msg = "Internal error - clusters must be factors")



  if (isEmptyName(method)) {
    method = "Likelyhood"
  }

  cellsMatrix <- NULL
  if (str_equal(method, "Likelyhood", ignore_case = TRUE)) {
    cellsMatrix <- calculateLikelihoodOfObserved(objCOTAN)
  } else if (str_equal(method, "LogNorm", ignore_case = TRUE)) {
    cellsMatrix <- getNormalizedData(objCOTAN, retLog = TRUE)
  } else {
    stop("Unrecognised `method` passed in: ", method)
  }

  if (is_empty(genesSel)) {
    genesSel <- "HGDI"
  }

  genesPos <- rep(TRUE, getNumGenes(objCOTAN))
  if (length(genesSel) > 1L) {
    genesPos <- getGenes(objCOTAN) %in% genesSel
  } else if (str_equal(genesSel, "HGDI", ignore_case = TRUE)) {
    gdi <- calculateGDI(objCOTAN, statType = "S", rowsFraction = 0.05)[["GDI"]]

    if(sum(gdi >= 1.5) > numGenes) {
      genesPos <- order(gdi, decreasing = TRUE)[1:numGenes]
    } else {
      genesPos <- gdi >= 1.4
    }
    rm(gdi)
  } else if (str_equal(genesSel, "HVG", ignore_case = TRUE)) {
    stop("Unsupported `geneSel`")
  } else {
    stop("Unrecognised `geneSel` passed in: ", genesSel)
  }

  cellsMatrix <- cellsMatrix[genesPos, ]

  pcaCells <- pca(mat = cellsMatrix, rank = 50L,
                  transposed = FALSE, BSPARAM = IrlbaParam())[["rotated"]]

  gc()

  pcaCells <- as.data.frame(pcaCells)

  umapPlot <- UMAPPlot(pcaCells, clusters = clusters, colors = colors,
                       title = paste("UMAP of clusterization", clName))

  return(list("plot" = umapPlot, "cellsPCA" = cellsPCA))
}
