
#' scatterPlot
#'
#' @description Plot that check the relation between the library size and the
#'   number of genes detected.
#'
#' @param objCOTAN a `COTAN` object
#' @param splitPattern Pattern used to extract, from the column names, the
#'   sample field (default " ")
#' @param numCols Once the column names are splitted by splitPattern, the column
#'   number with the sample name (default 2)
#' @param split.samples Boolean. Whether to plot each sample in a different
#'   panel (default `FALSE`)
#'
#' @returns The scatter plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 scale_y_log10
#'
#' @importFrom scales trans_breaks
#' @importFrom scales math_format
#' @importFrom scales trans_format
#'
#' @importFrom Matrix colSums
#'
#' @importFrom stringr str_split
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' scPlot <- scatterPlot(objCOTAN)
#' plot(scPlot)
#'
#' @rdname scatterPlot
#'
scatterPlot <- function(obj, splitPattern = " ",
                        numCols=2, split.samples = FALSE) {
  lib.size <- getCellsSize(objCOTAN)
  gene.size <- Matrix::colSums(getZeroOneProj(objCOTAN))

  to.plot <- cbind(lib.size, gene.size)
  to.plot <- as.data.frame(to.plot)
  to.plot$sample <- stringr::str_split(rownames(to.plot),pattern = splitPattern,simplify = T)[,numCols]

  plot <- ggplot(to.plot, aes(x=lib.size,y=gene.size, color = sample)) +
    geom_point(size = 0.5, alpha= 0.8) +
    labs(title = "Scatter plot of library size VS gene detected for each cell",
         y = "Gene number",
         x = "Library size (UMI)") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    plotTheme("size-plot")

  if (split.samples == T){
    plot <- plot + facet_grid(cols = vars(sample))
  }

  return(plot)
}
