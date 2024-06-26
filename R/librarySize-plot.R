
#' @importFrom ggplot2 ggproto
#' @importFrom ggplot2 Geom
#' @importFrom ggplot2 GeomPolygon
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 resolution
#' @importFrom ggplot2 draw_key_polygon
#' @importFrom ggplot2 layer
#'
#' @importFrom plyr arrange
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#'
#' @noRd
#'
geom_flat_violin <- function(
    mapping = NULL,
    data = NULL,
    stat = "ydensity",
    position = "dodge",
    trim = TRUE,
    scale = "area",
    show.legend = NA,
    inherit.aes = TRUE,
    ...) {
  GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,
     setup_data = function(data, params) {
       `%||%` <- function(a, b) {
         if (!is.null(a)) a else b
       }

       data$width <- data$width %||%
         params$width %||% (resolution(data$x, FALSE) * 0.9)

       # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
       data %>%
         group_by(group) %>%
         mutate(
           ymin = min(y),
           ymax = max(y),
           xmin = x,
           xmax = x + width / 2.0
         )
     },

     draw_group = function(data, panel_scales, coord) {
       # Find the points for the line to go all the way around
       data <- transform(data,
                         xminv = x,
                         xmaxv = x + violinwidth * (xmax - x)
       )

       # Make sure it's sorted properly to draw the outline
       newdata <- rbind(
         arrange(transform(data, x = xminv), y),
         arrange(transform(data, x = xmaxv), -y)
       )

       # Close the polygon: set first and last point the same
       # Needed for coord_polar and such
       newdata <- rbind(newdata, newdata[1L, ])

       ggplot2:::ggname("geom_flat_violin",
                        GeomPolygon$draw_panel(newdata, panel_scales, coord))
     },

     draw_key = draw_key_polygon,

     default_aes = aes(
       weight = 1.0, colour = "grey20", fill = "white", linewidth = 0.5,
       alpha = NA, linetype = "solid"
     ),

     required_aes = c("x", "y")
  )

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}



#' @details `cellSizePlot()` plots the raw library size for each cell and
#'   sample.
#'
#' @param objCOTAN a `COTAN` object
#' @param splitPattern Pattern used to extract, from the column names, the
#'   sample field (default " ")
#' @param numCol Once the column names are split by splitPattern, the column
#'   number with the sample name (default 2)
#'
#' @returns `cellSizePlot()` returns the `violin-boxplot` plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 position_jitter
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 ylim
#'
#' @importFrom stringr str_split
#'
#' @export
#'
#' @examples
#' lsPlot <- cellSizePlot(objCOTAN)
#' plot(lsPlot)
#'
#' @rdname RawDataCleaning
#'
cellSizePlot <- function(objCOTAN, splitPattern = " ", numCol = 2L) {
  sizes <- sort(getCellsSize(objCOTAN))
  sizes <- as.data.frame(sizes)
  sizes <- setColumnInDF(sizes, seq_len(nrow(sizes)), colName = "n")
  if (TRUE) {
    splitNames <- str_split(rownames(sizes),
                            pattern = splitPattern, simplify = TRUE)
    if (ncol(splitNames) < numCol) {
      # no splits found take all as a single group
      sampleCol <- rep("1", nrow(splitNames))
    } else {
      sampleCol <- splitNames[, numCol]
    }
    sizes <- setColumnInDF(sizes, sampleCol, colName = "sample")
  }

  plot <-
    sizes %>%
    ggplot(aes(x = sample, y = sizes, fill = sample)) +
    #geom_point(size = 0.5) +
    geom_flat_violin(position = position_nudge(x = 0.15, y = 0.0),
                     adjust = 2.0, alpha = 0.5) +
    geom_point(position = position_jitter(width = 0.1),
               size = 0.4, color = "black", alpha = 0.5) +
    geom_boxplot(aes(x = sample, y = sizes, fill = sample),
                 outlier.shape = NA, alpha = 0.8, width = 0.15,
                 colour = "gray65", size = 0.6) +
    labs(title = "Cell library size",
         y = "Size (read number)",
         x = "") +
    # scale_y_continuous(expand = c(0L, 0L)) +
    ylim(0.0, max(sizes[["sizes"]])) +
    plotTheme("size-plot")

  return(plot)
}


#' @details `genesSizePlot()` plots the raw gene number (reads > 0) for each
#'   cell and sample
#'
#' @param objCOTAN a `COTAN` object
#' @param splitPattern Pattern used to extract, from the column names, the
#'   sample field (default " ")
#' @param numCol Once the column names are split by splitPattern, the column
#'   number with the sample name (default 2)
#'
#' @returns `genesSizePlot()` returns the `violin-boxplot` plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 position_nudge
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 position_jitter
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 ylim
#'
#' @importFrom stringr str_split
#'
#' @export
#'
#' @examples
#' gsPlot <- genesSizePlot(objCOTAN)
#' plot(gsPlot)
#'
#' @rdname RawDataCleaning
#'
genesSizePlot <- function(objCOTAN, splitPattern = " ", numCol = 2L) {
  sizes <- sort(getNumExpressedGenes(objCOTAN))
  sizes <- as.data.frame(sizes)
  sizes <- setColumnInDF(sizes, seq_len(nrow(sizes)), colName = "n")
  if (TRUE) {
    splitNames <- str_split(rownames(sizes),
                            pattern = splitPattern, simplify = TRUE)
    if (ncol(splitNames) < numCol) {
      # no splits found take all as a single group
      sampleCol <- rep("1", nrow(splitNames))
    } else {
      sampleCol <- splitNames[, numCol]
    }
    sizes <- setColumnInDF(sizes, sampleCol, colName = "sample")
  }

  plot <-
    sizes %>%
    ggplot(aes(x = sample, y = sizes, fill = sample)) +
    #geom_point(size = 0.5) +
    geom_flat_violin(position = position_nudge(x = 0.15, y = 0.0),
                     adjust = 2.0, alpha = 0.5) +
    geom_point(position = position_jitter(width = 0.1),
               size = 0.4, color = "black", alpha = 0.5) +
    geom_boxplot(aes(x = sample, y = sizes, fill = sample),
                 outlier.shape = NA, alpha = 0.8, width = 0.15,
                 colour = "gray65", size = 0.6) +
    labs(title = "Detected gene number",
         y = "Size (number of genes)",
         x = "") +
    # scale_y_continuous(expand = c(0L, 0L)) +
    ylim(0L, max(sizes[["sizes"]])) +
    plotTheme("size-plot")

  return(plot)
}
