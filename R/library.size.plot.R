### This script creates an R function to generate raincloud plots, then simulates
### data for plots. If using for your own data, you only need lines 1-80.
### It relies largely on code previously written by David Robinson
### (https://gist.github.com/dgrtwo/eb7750e74997891d7c20)
### and the package ggplot2 by Hadley Wickham

# Check if required packages are installed ----
#packages <- c("cowplot", "readr", "ggplot2", "dplyr", "lavaan", "Hmisc")
#if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#  install.packages(setdiff(packages, rownames(installed.packages())))
#}

# Load packages ----

# Defining the geom_flat_violin function ----
# Note: the below code modifies the
# existing github page by removing a parenthesis in line 50

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
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

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
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
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )



#' library.size.plot
#' 
#' Function that plots the row library size for each cell and sample (if there 
#' are more samples).
#'
#' @param obj COTAN object
#' @param split.pattern pattern used to extract, from the column names, the 
#' sample field (default " ")
#' @param n.col ones the column names are splitted by split.pattern, the column 
#' number with the sample name (default 2)
#' 
#'@importFrom Matrix colSums
#'@importFrom ggthemes theme_tufte
#'@import ggplot2
#'@importFrom stringr str_split
#' 
#' @return the violin-boxplot plot
#' @export
#'
#' @examples
library.size.plot <- function(obj, split.pattern = " ", n.col=2){
  sizes <- Matrix::colSums(obj@raw)
  sizes <- sort(sizes)
  sizes <- as.data.frame(sizes)
  sizes$n <- c(1:dim(sizes)[1])
  sizes$sample <- str_split(rownames(sizes),pattern = split.pattern,simplify = T)[,n.col]

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
    ylim(0,max(sizes$sizes))+
    ggthemes::theme_tufte()+
    theme(legend.position = "none")#,
  #axis.text.x=element_blank(),
  #axis.ticks.x=element_blank())
  
  return(plot)
   
}


  
