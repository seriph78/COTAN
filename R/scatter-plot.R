
#' scatterPlot
#'
#' @details `scatterPlot()` creates a plot that check the relation between the
#'   library size and the number of genes detected.
#'
#' @param objCOTAN a `COTAN` object
#' @param splitPattern Pattern used to extract, from the column names, the
#'   sample field (default " ")
#' @param numCol Once the column names are split by splitPattern, the column
#'   number with the sample name (default 2)
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
scatterPlot <- function(objCOTAN, splitPattern = " ",
                        numCol = 2L, splitSamples = FALSE) {
  cellsSize <- getCellsSize(objCOTAN)
  genesSize <- getNumExpressedGenes(objCOTAN)

  toPlot <- cbind(cellsSize, genesSize)
  toPlot <- as.data.frame(toPlot)

  if (TRUE) {
    splitNames <- str_split(rownames(toPlot),
                            pattern = splitPattern, simplify = TRUE)
    if (ncol(splitNames) < numCol) {
      # no splits found take all as a single group
      sampleCol <- rep("1", nrow(splitNames))
    } else {
      sampleCol <- splitNames[, numCol]
    }
    toPlot <- setColumnInDF(toPlot, sampleCol, colName = "sample")
  }

  plot <- ggplot(toPlot, aes(x = cellsSize, y = genesSize, color = sample)) +
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
