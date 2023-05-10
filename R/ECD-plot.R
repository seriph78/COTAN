

#' @details `ECDPlot()` plots the empirical distribution function of library
#'   sizes (UMI number). It helps to define where to drop "cells" that are
#'   simple background signal.
#'
#' @param objCOTAN a `COTAN` object
#' @param yCut y threshold of library size to drop
#'
#' @returns `ECDPlot()` returns an ECD plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#'
#' @export
#'
#' @examples
#' ## These plots might help to identify genes/cells that need to be dropped
#' ecdPlot <- ECDPlot(objCOTAN, yCut = 100)
#' plot(ecdPlot)
#'
#' @rdname RawDataCleaning
#'
ECDPlot <- function(objCOTAN, yCut) {
  libSize <- sort(getCellsSize(objCOTAN), decreasing = TRUE)
  df <- setColumnInDF(data.frame(), colToSet = libSize, colName = "libSize")
  n <- seq_along(libSize)
  df <- setColumnInDF(df, colToSet = n, colName = "n")
  ggplot(df, aes(y = log(libSize), x = log(n))) +
    geom_point() +
    geom_hline(yintercept = log(yCut), linetype = "dashed", color = "red")
}
