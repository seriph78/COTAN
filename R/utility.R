#' plotTheme
#'
#' This function returns the appropriate theme for the selected plot kind
#'
#' @param plotKind a string indicating the plot kind; supported kinds are:
#' "common", "pca", "genes", "UDE", "heatmap", "GDI", "size-plot"
#' @param textSize axes and strip text size (default=14)
#' @return a ggplot2 theme object
#
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 unit
#'
#' @importFrom ggthemes theme_tufte
#'
#' @export
#' @rdname plotTheme
plotTheme <- function(plotKind = "common", textSize = 14) {
  myDarkBlue <- "#3C5488FF"
  ts <- textSize

  basic_theme <- theme(
    axis.text.x  = element_text(size = ts, angle = 0,  hjust = .5, vjust = .5, face = "plain", colour = myDarkBlue),
    axis.text.y  = element_text(size = ts, angle = 0,  hjust = 0,  vjust = .5, face = "plain", colour = myDarkBlue),
    axis.title.x = element_text(size = ts, angle = 0,  hjust = .5, vjust = 0,  face = "plain", colour = myDarkBlue),
    axis.title.y = element_text(size = ts, angle = 90, hjust = .5, vjust = .5, face = "plain", colour = myDarkBlue) )

  if (plotKind == "common") {
    return(basic_theme)
  }

  if (plotKind == "pca") {
    return( basic_theme +
            theme(legend.title = element_blank(),
                  legend.text = element_text(size = 12, face = "italic", color = myDarkBlue),
                  legend.position = "bottom") )
  }

  if (plotKind == "genes") {
    return( basic_theme +
            theme(plot.title = element_text(size = 20, hjust = 0.02, vjust = -10,
                                            face = "italic", color = myDarkBlue),
                  plot.subtitle = element_text(vjust = -15, hjust = 0.01,
                                               color = "darkred")) )
  }

  if (plotKind == "UDE") {
    return( basic_theme +
            theme(plot.title   = element_text(size = 20, color = myDarkBlue),
                  legend.title = element_text(size = 14, color = myDarkBlue, face = "italic"),
                  legend.text  = element_text(size = 11, color = myDarkBlue),
                  legend.key.width = unit(2, "mm"),
                  legend.position  ="right") )
  }

  if (plotKind == "heatmap") {
    return( basic_theme +
            theme( axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   panel.spacing = unit(0, "lines"),
                   strip.background = element_rect(fill = "#8491B44C"),
                   strip.text.y = element_text(size = ts, colour = myDarkBlue),
                   strip.text.x = element_text(size = ts, angle = 90, colour = myDarkBlue),
                   legend.text = element_text(color = myDarkBlue, face = "italic"),
                   legend.position = "bottom",
                   legend.title = element_blank(),
                   legend.key.height = unit(2, "mm") ) )
  }

  if (plotKind == "GDI") {
    return( basic_theme +
            theme( legend.title = element_blank(),
                   plot.title = element_text(size = ts + 2, face = "bold.italic", color = myDarkBlue),
                   legend.text = element_text(color = myDarkBlue, face = "italic"),
                   legend.position = "bottom") )
  }

  if (plotKind == "size-plot") {
    return( ggthemes::theme_tufte() +
            theme(legend.position = "none") )
                  # axis.text.x  = element_blank(),
                  # axis.ticks.x = element_blank()) )
  }

  warning(paste("plotTheme: no match found in listed themes for:", plotKind))
  return(basic_theme)
}

