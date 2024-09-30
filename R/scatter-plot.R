
#' scatterPlot
#'
#' @details `scatterPlot()` creates a plot that check the relation between the
#'   library size and the number of genes detected.
#'
#' @param objCOTAN a `COTAN` object
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#' @param splitSamples Boolean. Whether to plot each sample in a different panel
#'   (default `FALSE`)
#'
#' @returns `scatterPlot()` returns the scatter plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 scale_y_log10
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 vars
#'
#' @importFrom scales trans_breaks
#' @importFrom scales math_format
#' @importFrom scales trans_format
#'
#' @importFrom stringr str_split
#'
#' @export
#'
#' @examples
#' scPlot <- scatterPlot(objCOTAN)
#' plot(scPlot)
#'
#' @rdname RawDataCleaning
#'
scatterPlot <-
  function(objCOTAN, condName = "", conditions = NULL, splitSamples = TRUE) {
  cellsSize <- getCellsSize(objCOTAN)
  genesSize <- getNumExpressedGenes(objCOTAN)

  sizes <- cbind(cellsSize, genesSize)
  sizes <- as.data.frame(sizes)

  c(., conditions) %<-%
    normalizeNameAndLabels(objCOTAN, name = condName,
                           labels = conditions, isCond = TRUE)
  assert_that(!is_empty(conditions))

  sizes <- setColumnInDF(sizes, conditions[rownames(sizes)], colName = "sample")

  plot <- ggplot(sizes, aes(x = cellsSize, y = genesSize, color = sample)) +
    geom_point(size = 0.5, alpha = 0.8) +
    labs(title = "Scatter plot of library size VS gene detected for each cell",
         y = "Gene number",
         x = "Library size (UMI)") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    plotTheme("size-plot")

  if (isTRUE(splitSamples)) {
    plot <- plot + facet_grid(cols = vars(sample))
  }

  return(plot)
}
