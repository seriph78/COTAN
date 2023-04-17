#' @details `UMAPPlot()` plots the given `data.frame` containing genes
#'   information related to cleusters after applying the UMAP transformation.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes a named `list` of genes to label. Each array will have different
#'   color.
#' @param df the `data.frame` to plot. It must have a row for each gene.
#' @param cond a string corresponding to the condition/sample (it is used only
#'   for the title).
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
  emptyTitle <- !(length(title) && any(nzchar(title)))
  if (emptyTitle) {
    title <- "UMAP Plot"
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
  for (cl in unique(clusters)) {
    selec <- !labelled & clusters == cl
    if (any(selec)) {
      colors[selec] <- cl
    } else {
      logThis(paste("UMAPPlot - none of the elements of the cluster", cl,
                    "is present: will be ignored"), logLevel = 1L)
    }
  }

  clustered <- !labelled & colors != "none"

  umap <- umap(df)

  plotDF <- data.frame(x = umap[["layout"]][,1],
                       y = umap[["layout"]][,2],
                       colors = colors)

  allTypes <- setdiff(unique(colors), c("none"))
  myColours <- set_names(getColorsVector(length(allTypes)), allTypes)

  plot <- ggplot(subset(plotDF, (!labelled & !clustered))) +
    geom_point(aes(x, y, colour = "#8491B4B2"),
               size = 1.5, alpha = 0.3) +
    geom_point(data = subset(plotDF, clustered),
               aes(x, y, colour = colors),
               size = 1.5, alpha = 0.5) +
    geom_point(data = subset(plotDF, labelled),
               aes(x, y, colour = colors),
               size = 2.0, alpha = 0.8) +
    scale_color_manual("Status", values = myColours) +
    scale_fill_manual( "Status", values = myColours) +
    xlab("") + ylab("") +
    geom_label_repel(data = subset(plotDF, labelled),
                     aes(x, y, fill = colors),
                     label = rownames(plotDF)[labelled],
                     label.size = NA, max.overlaps = 40L, alpha = 0.8,
                     direction = "both", na.rm = TRUE, seed = 1234L) +
    ggtitle(title) +
    plotTheme("UMAP", textSize = 10L)

  return(plot)
}
