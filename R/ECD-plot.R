

#' @details `ECDPlot()` plots the *Empirical Cumulative Distribution* function
#'   of library sizes (`UMI` number). It helps to define where to drop "cells"
#'   that are simple background signal.
#'
#' @param objCOTAN a `COTAN` object
#' @param yCut y threshold of library size to drop. Default is `NaN`
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#'
#' @returns `ECDPlot()` returns an `ECD` plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#'
#' @export
#'
#' @examples
#' ## These plots might help to identify genes/cells that need to be dropped
#' ecdPlot <- ECDPlot(objCOTAN, yCut = 100.0)
#' plot(ecdPlot)
#'
#' @rdname RawDataCleaning
#'
ECDPlot <- function(objCOTAN, yCut = NaN,
                    condName = "", conditions = NULL) {
  libSize <- sort(getCellsSize(objCOTAN), decreasing = TRUE)
  sizes <- setColumnInDF(data.frame(), colToSet = libSize, colName = "libSize")

  n <- seq_along(libSize)
  sizes <- setColumnInDF(sizes, colToSet = n, colName = "n")

  c(., conditions) %<-%
    normalizeNameAndLabels(objCOTAN, name = condName,
                           labels = conditions, isCond = TRUE)
  assert_that(!is_empty(conditions))

  sizes <- setColumnInDF(sizes, conditions[rownames(sizes)], colName = "sample")

  plot <- ggplot(sizes, aes(x = log(n), y = log(libSize),
                            fill = sample, colour = sample)) +
    geom_point()

  if (!is.nan(yCut)) {
    plot <- plot +
      geom_hline(yintercept = log(yCut), linetype = "dashed", color = "red")
  }
  return(plot)
}
