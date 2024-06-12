
#' @details `UMAPPlot()` plots the given `data.frame` containing genes
#'   information related to clusters after applying the UMAP transformation.
#'
#' @param df The `data.frame` to plot. It must have a row names containing the
#'   given elements
#' @param clusters The **clusterization**. Must be aligned to the rows in the
#'   `data.frame`
#' @param elements a named `list` of elements to label. Each array in the list
#'   will have different color. When not given it will label the clusters using
#'   their centroids
#' @param title a string giving the plot title. Will default to UMAP Plot if not
#'   specified
#'
#'
#' @returns `UMAPPlot()` returns a `ggplot2` object
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
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
UMAPPlot <- function(df, clusters = NULL, elements = NULL, title = "") {
  logThis("UMAP plot", logLevel = 2L)

  assert_that(is_empty(clusters) || length(clusters) == nrow(df),
              msg = paste("Clusters vector must have size equal to",
                          "the number of rows in the data.frame"))

  # empty title
  if (isEmptyName(title)) {
    title <- "UMAP Plot"
  }

  if (!is_empty(clusters)) {
    clusters <- factor(clusters)
  }

  colors <- rep_len("none", nrow(df))

  # assign a different color to each list of elements
  for (nm in names(elements)) {
    selec <- rownames(df) %in% elements[[nm]]
    if (any(selec)) {
      colors[selec] <- nm
    } else {
      logThis(paste("UMAPPlot - none of the elements in group", nm,
                    "is present: will be ignored"), logLevel = 1L)
    }
  }

  labelled <- colors != "none"

  # assign a different color to each cluster
  for (cl in levels(clusters)) {
    selec <- !labelled & clusters == cl
    if (any(selec)) {
      colors[selec] <- cl
    } else {
      logThis(paste("UMAPPlot - none of the elements of the cluster", cl,
                    "is present unlabelled: cluster will be ignored"),
              logLevel = 1L)
    }
  }

  clustered <- !labelled & colors != "none"

  umap <- umap(df)

  plotDF <- data.frame(x = umap[["layout"]][, 1L],
                       y = umap[["layout"]][, 2L],
                       colors = colors)

  centroids <- NULL
  if (is_empty(elements) && !is_empty(clusters)) {
    centroids <- data.frame()

    clList <- toClustersList(clusters)
    for (clName in names(clList)) {
      row <- colMeans(plotDF[clList[[clName]], 1:2, drop = FALSE])
      centroids <- rbind(centroids, row)
    }
    colnames(centroids) <- c("x", "y")

    centroids <- setColumnInDF(centroids, colName = "colors",
                               colToSet = names(clList))
  }

  allTypes <- setdiff(unique(colors), "none")
  myColours <- set_names(getColorsVector(length(allTypes)), allTypes)

  pointSize <- min(max(1.0, 200000.0/dim(plotDF)[1]), 5)

  plot <- ggplot(subset(plotDF, (!labelled & !clustered))) +
    geom_point(aes(x, y, colour = "#8491B4B2"),
               size = pointSize, alpha = 0.3) +
    geom_point(data = subset(plotDF, clustered),
               aes(x, y, colour = colors),
               size = pointSize, alpha = 0.5) +
    geom_point(data = subset(plotDF, labelled),
               aes(x, y, colour = colors),
               size = 1.5*pointSize, alpha = 0.8) +
    scale_color_manual("Status", values = myColours) +
    scale_fill_manual( "Status", values = myColours) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    geom_label_repel(data = subset(plotDF, labelled),
                     aes(x, y, fill = colors),
                     label = rownames(plotDF)[labelled],
                     label.size = NA, max.overlaps = 40L, alpha = 0.8,
                     direction = "both", na.rm = TRUE, seed = 1234L) +
    ggtitle(title) +
    plotTheme("UMAP", textSize = 10L)

  if (!is.null(centroids)) {
    plot <- plot +
      geom_text_repel(data = centroids,
                      aes(x, y, label = colors, colour = colors),
                      fontface = "bold",
                      size = 0.1*pointSize)
  }

  return(plot)
}
