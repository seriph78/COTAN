funProbZero <- function(disp, mu) {
  ad <- abs(disp)
  (disp <= 0) * (exp(-(1 + ad) * mu)) +
    (disp >  0) * (1 + ad * mu)^(-1 / ad)
}

#' dispersionBisection
#'
#' private function invoked by 'estimateDispersion' for the estimation
#' of 'dispersion' field of a COTAN object via a bisection solver
#'
#' the goal is to find dispersion value that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param gene name of the relevant gene
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param muEstimator a matrix estimator of vector mu
#' @param threshold minimal solution precision
#'
#' @return r, data.frame(a, u)
#'
#' @rdname dispersionBisection
dispersionBisection <- function(gene,
                                zeroOneMatrix,
                                muEstimator,
                                housekeepingGenes,
                                threshold = 0.001) {
  if(gene %in% housekeepingGenes) {
    return(NA)
  }

  sumZeros <- sum(zeroOneMatrix[gene, ] == 0)
  muEstimator <- muEstimator[gene, ]

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  disp1 <- 0
  diff1 <- sum(funProbZero(disp1, muEstimator)) - sumZeros
  if (abs(diff1) <= threshold) {
    return(disp1)
  }

  disp2 <- -1 * sign(diff1) # we assume error is an increasing function of disp
  repeat {
    diff2 <- sum(funProbZero(disp2, muEstimator)) - sumZeros

    if (diff2 * diff1 < 0) {
      break
    }

    disp1 <- disp2 # disp is closer to producing 0
    diff1 <- diff2

    disp2 <- 2 * disp2 # we double at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  repeat {
    disp <- (disp1 + disp2) / 2
    diff <- sum(funProbZero(disp, muEstimator)) - sumZeros

    if (abs(diff) <= threshold) {
      return(disp)
    }

    # drop same sign diff point
    if (diff * diff2 > 0) {
      disp2 <- disp
      diff2 <- diff
    } else {
      disp1 <- disp
      diff1 <- diff
    }
  }
}


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

