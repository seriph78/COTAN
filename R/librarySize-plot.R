
#' @importFrom ggplot2 ggproto
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 draw_key_polygon
#' @importFrom ggplot2 layer
#'
#' @importFrom plyr arrange
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
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
  GeomFlatViolin <- ggplot2::ggproto("GeomFlatViolin", Geom,
     setup_data = function(data, params) {
       "%||%" <- function(a, b) {
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
           xmax = x + width / 2
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
         plyr::arrange(transform(data, x = xminv), y),
         plyr::arrange(transform(data, x = xmaxv), -y)
       )

       # Close the polygon: set first and last point the same
       # Needed for coord_polar and such
       newdata <- rbind(newdata, newdata[1, ])

       ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
     },

     draw_key = ggplot2::draw_key_polygon,

       default_aes = ggplot2::aes(
         weight = 1, colour = "grey20", fill = "white", linewidth = 0.5,
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



#' librarySizePlot
#'
#' @description Function that plots the row library size for each cell and
#'   sample (if there are more samples).
#'
#' @param objCOTAN a `COTAN` object
#' @param splitPattern Pattern used to extract, from the column names, the
#'   sample field (default " ")
#' @param numCols Once the column names are splitted by splitPattern, the column
#'   number with the sample name (default 2)
#'
#' @returns the violin-boxplot plot
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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' lsPlot <- librarySizePlot(objCOTAN)
#' plot(lsPlot)
#'
#' @rdname librarySizePlot
#'
librarySizePlot <- function(objCOTAN, splitPattern = " ", numCols = 2) {
  sizes <- getCellsSize(objCOTAN)
  sizes <- sort(sizes)
  sizes <- as.data.frame(sizes)
  sizes$n <- c(1:dim(sizes)[1])
  sizes$sample <- str_split(rownames(sizes), pattern = splitPattern, simplify = T)[, numCols]

  plot <-
    sizes %>%
    ggplot(aes(x = sample, y = sizes, fill = sample)) +
    #geom_point(size = 0.5) +
    geom_flat_violin(position = position_nudge(x = .15, y = 0),adjust =2, alpha = 0.5) +
    geom_point(position = position_jitter(width = .1), size = .4, color="black", alpha= 0.5) +
    geom_boxplot(aes(x = sample, y = sizes, fill = sample),
                 outlier.shape = NA, alpha = .8, width = .15, colour = "gray65", size = 0.6)+
    labs(title = "Cell library size",
         y = "Size (read number)",
         x = "") +
    scale_y_continuous(expand = c(0, 0)) +
    ylim(0,max(sizes$sizes)) +
    plotTheme("size-plot")

  return(plot)
}



