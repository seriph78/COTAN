
#' @details `UMAPPlot()` plots the given `data.frame` containing genes
#'   information related to clusters after applying the umap transformation via
#'   [Seurat::RunUMAP()]
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
#' @param numNeighbors Overrides the `n_neighbors` value from
#'   [umap::umap.defaults]
#' @param minPointsDist Overrides the `min_dist` value from
#'   [umap::umap.defaults]
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
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat Embeddings
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
UMAPPlot <- function(df,
                     clusters = NULL,
                     elements = NULL,
                     title = "",
                     colors = NULL,
                     numNeighbors = 0L,
                     minPointsDist = NaN) {
  logThis("UMAP plot", logLevel = 2L)

  assert_that(!is_empty(rownames(df)),
              msg = "UMAPPlot - data.frame must have proper row-names")

  assert_that(is_empty(clusters) || identical(names(clusters), rownames(df)),
              msg = paste("UMAPPlot - clusters' names must be the same",
                          "as the row-names of the data.frame"))

  # empty title
  if (isEmptyName(title)) {
    title <- "UMAP Plot"
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
  clusters <- factor(clusters)
  for (cl in levels(clusters)) {
    selec <- !labelled & clusters == cl
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

  if (numNeighbors == 0L) {
    numNeighbors <- 30L
  }
  if (!is.finite(minPointsDist)) {
    minPointsDist <- 0.3
  }

  logThis("Calculating UMAP: START", logLevel = 3L)

  umap <- suppressWarnings(
    Embeddings(RunUMAP(df,
                       assay = "Generic",
                       reduction.key = "UMAP_",
                       umap.method = "uwot",
                       n.neighbors = numNeighbors, #30L
                       n.components = 2L,
                       metric = "cosine",
                       learning.rate = 1,
                       min.dist = minPointsDist, # 0.3
                       uwot.sgd = FALSE,
                       seed.use = 42,
                       verbose = TRUE)))

  plotDF <- data.frame(x = umap[, 1L],
                       y = umap[, 2L],
                       types = entryType)

  logThis("Calculating UMAP: DONE", logLevel = 3L)

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
    if (numMissing > 1L) {
      myColours[["none"]] <- "#8491B4B2"
    }
    myColours[["centroid"]] <- "#000000"
  } else {
    myColours <- myColours[seq_along(allTypes)]
    names(myColours) <- allTypes
  }

  assert_that(setequal(c(entryType, "none", "centroid"), names(myColours)))

  pointSize <- min(max(0.33, 5000.0 / nrow(plotDF)), 3.0)

  plot <- ggplot() +
    scale_color_manual("Status", values = myColours) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    ggtitle(title) +
    plotTheme("UMAP", textSize = 10L)

  if (any(generic)) {
    plot <- plot +
      geom_point(data = plotDF[generic, , drop = FALSE],
                 aes(x, y, colour = types),
                 size = pointSize, alpha = 0.3)
  }

  if (any(clustered | centroids)) {
    plot <- plot +
      geom_point(data = plotDF[clustered, , drop = FALSE],
                 aes(x, y, colour = types),
                 size = pointSize, alpha = 0.5) +
      geom_text_repel(data = plotDF[centroids, , drop = FALSE],
                      aes(x, y, colour = "centroid"),
                      label = rownames(plotDF)[centroids],
                      max.overlaps = 40L, fontface = "bold",
                      show.legend = FALSE, alpha = 0.8)
  }
  if (any(labelled | clustered | centroids)) {
    plot <- plot +
      scale_fill_manual("Status", values = myColours)
  }

  if (any(labelled)) {
    plot <- plot +
      geom_point(data = plotDF[labelled, , drop = FALSE],
                 aes(x, y, colour = types),
                 size = 1.2 * pointSize, alpha = 0.8) +
      geom_label_repel(data = plotDF[labelled, , drop = FALSE],
                       aes(x, y, fill = types,
                           label = rownames(plotDF)[labelled]),
                       label.size = NA, show.legend = FALSE, force = 2.0,
                       box.padding = 0.25,
                       max.overlaps = 40L, alpha = 0.8,
                       direction = "both", na.rm = TRUE, seed = 1234L)
  } else {
    plot <- plot +
      theme(legend.position = "none")
  }

  return(plot)
}

#' @title Get High Variable Genes running the `Seurat` package
#'
#' @description The function uses the [Seurat::Seurat] package to extract the
#'   high variable genes given the counts raw data
#'
#' @param rawData the raw counts
#' @param hvgMethod the HVG method
#' @param cond the sample condition
#' @param numFeatures the number of genes to select. It defaults to 2000
#'
#' @returns a subset of the genes in the dataset
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat VariableFeatures
#'
#' @importFrom withr with_options
#'
#' @noRd
#'
seuratHVG <- function(rawData, hvgMethod, cond, numFeatures = 2000L) {
  tryCatch({
    logThis("Creating Seurat object: START", logLevel = 2L)

    srat <- CreateSeuratObject(counts = as.data.frame(rawData),
                               project = paste0(cond, "_UMAP"),
                               min.cells = 3L, min.features = 4L)
    srat <- NormalizeData(srat)
    srat <- FindVariableFeatures(srat, selection.method = hvgMethod,
                                 nfeatures = numFeatures)

    hvg <- VariableFeatures(object = srat)

    rm(srat)
    logThis("Creating Seurat object: DONE", logLevel = 2L)

    return(hvg)
  },
  error = function(e) {
    logThis(msg = paste("Seurat HVG failed with", nrow(rawData),
                        "genes with the following error:"), logLevel = 1L)
    logThis(msg = conditionMessage(e), logLevel = 1L)
    return(NULL)
  })
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
#' @param dataMethod selects the method to use to create the `data.frame` to
#'   pass to the [UMAPPlot()]. To calculate, for each cell, a statistic for each
#'   gene based on available data/model, the following methods are supported:
#'   * `"NuNorm"` uses the \eqn{\nu}*-normalized* counts
#'   * `"LogNormalized"` uses the *log-normalized* counts. The default method
#'   * `"Likelihood"` uses the likelihood of observed presence/absence of each
#'     gene
#'   * `"LogLikelihood"` uses the likelihood of observed presence/absence of
#'     each gene
#'   * `"Binarized"` uses the binarized data matrix
#'   * `"AdjBinarized"` uses the binarized data matrix where ones and zeros
#'     are replaced by the per-gene estimated probability of zero and its
#'     complement respectively
#' @param genesSel Decides whether and how to perform gene-selection. It can be
#'   a straight list of genes or a string indicating one of the following
#'   selection methods:
#'   * `"HGDI"` Will pick-up the genes with highest **GDI**. Since it requires
#'     an available `COEX` matrix it will fall-back to `"HVG_Seurat"` when the
#'     matrix is not available
#'   * `"HVG_Seurat"` Will pick-up the genes with the highest variability
#'     via the \pkg{Seurat} package (the default method)
#'   * `"HVG_Scanpy"` Will pick-up the genes with the highest variability
#'     according to the `Scanpy` package (using the \pkg{Seurat} implementation)
#' @param numGenes the number of genes to select using the above method. Will be
#'   ignored when no selection have been asked or when an explicit list of genes
#'   was passed in
#' @param colors an `array` of colors to use in the plot. If not sufficient
#'   colors are given it will complete the list using colors from
#'   [getColorsVector()]
#' @param numNeighbors Overrides the `n_neighbors` value from
#'   [umap::umap.defaults]
#' @param minPointsDist Overrides the `min_dist` value from
#'   [umap::umap.defaults]
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
#' @importFrom BiocSingular runPCA
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
cellsUMAPPlot <- function(objCOTAN,
                          clName = "",
                          clusters = NULL,
                          dataMethod = "",
                          genesSel = "HVG_Seurat",
                          numGenes = 2000L,
                          colors = NULL,
                          numNeighbors = 0L,
                          minPointsDist = NA) {
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

  if (isEmptyName(dataMethod)) {
    dataMethod <- "LogNormalized"
  }

  cellsMatrix <- NULL
  if (str_equal(dataMethod, "Binarized", ignore_case = TRUE)) {
    cellsMatrix <- getZeroOneProj(objCOTAN)
  } else if (str_equal(dataMethod, "AdjBinarized", ignore_case = TRUE)) {
    zeroOne <- getZeroOneProj(objCOTAN)
    rwMns <- rowMeans(zeroOne)
    cellsMatrix <- (zeroOne * (1.0 - rwMns) + (1.0 - zeroOne) * rwMns)
  } else if (str_equal(dataMethod, "Likelihood", ignore_case = TRUE)) {
    cellsMatrix <- calculateLikelihoodOfObserved(objCOTAN)
  } else if (str_equal(dataMethod, "LogLikelihood", ignore_case = TRUE)) {
    cellsMatrix <- log(calculateLikelihoodOfObserved(objCOTAN))
  } else if (str_equal(dataMethod, "NuNorm", ignore_case = TRUE) ||
             str_equal(dataMethod, "Normalized", ignore_case = TRUE)) {
    cellsMatrix <- getNuNormData(objCOTAN)
  } else if (str_equal(dataMethod, "LogNorm", ignore_case = TRUE) ||
             str_equal(dataMethod, "LogNormalized", ignore_case = TRUE)) {
    cellsMatrix <- getLogNormData(objCOTAN)
  } else {
    stop("Unrecognised `dataMethod` passed in: ", dataMethod)
  }

  genesPos <- rep(TRUE, getNumGenes(objCOTAN))
  if (length(genesSel) > 1L) {
    genesPos <- getGenes(objCOTAN) %in% genesSel
    logThis(paste("Given", sum(genesPos), "genes as input"), logLevel = 2L)
  } else {
    if (str_equal(genesSel, "HGDI", ignore_case = TRUE) &&
        !isCoexAvailable(objCOTAN)) {
      logThis(paste("The COEX matrix is not available: falling back",
                    "to HVG_Seurat for genes' selection"), logLevel = 1L)
      genesSel <- "HVG_Seurat"
    }

    if (str_equal(genesSel, "HGDI", ignore_case = TRUE)) {
      gdi <- getGDI(objCOTAN)
      if (is_empty(gdi)) {
        gdi <- getColumnFromDF(calculateGDI(objCOTAN, statType = "S",
                                            rowsFraction = 0.05), "GDI")
      }
      if (sum(gdi >= 1.5) > numGenes) {
        genesPos <-
          seq_along(gdi) %in% order(gdi, decreasing = TRUE)[1L:numGenes]
      } else {
        genesPos <- gdi >= 1.4
      }
      rm(gdi)
    } else if (str_equal(genesSel, "HVG_Seurat", ignore_case = TRUE)) {
      condition <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])
      genesPos <-
        getGenes(objCOTAN) %in% seuratHVG(getRawData(objCOTAN),
                                          hvgMethod = "vst",
                                          cond = condition,
                                          numFeatures = numGenes)
      rm(condition)
    } else if (str_equal(genesSel, "HVG_Scanpy", ignore_case = TRUE)) {
      condition <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])
      genesPos <-
        getGenes(objCOTAN) %in% seuratHVG(getRawData(objCOTAN),
                                          hvgMethod = "mean.var.plot",
                                          cond = condition,
                                          numFeatures = numGenes)
      rm(condition)
    } else {
      stop("Unrecognised `genesSel` passed in: ", genesSel)
    }
    logThis(paste("Selected", sum(genesPos), "genes using",
                  genesSel, "selector"), logLevel = 2L)
  }

  cellsMatrix <- cellsMatrix[genesPos, ]

  logThis("Elaborating PCA - START", logLevel = 3L)
  cellsPCA <- runPCA(x = t(cellsMatrix), rank = 50L,
                     BSPARAM = IrlbaParam(), get.rotation = FALSE)[["x"]]

  gc()

  cellsPCA <- as.data.frame(cellsPCA)
  logThis("Elaborating PCA - END", logLevel = 3L)

  umapTitle <- paste("UMAP of clusterization", clName, "using", dataMethod,
                     "matrix with", genesSel, "genes selector")
  umapPlot <- UMAPPlot(cellsPCA,
                       clusters = clusters,
                       colors = colors,
                       numNeighbors = numNeighbors,
                       minPointsDist = minPointsDist,
                       title = umapTitle)

  return(list("plot" = umapPlot, "cellsPCA" = cellsPCA))
}
