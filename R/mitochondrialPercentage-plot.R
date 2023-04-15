
#' @details `mitochondrialPercentagePlot()` plots the raw library size for each
#'   cell and sample.
#'
#' @param objCOTAN a `COTAN` object
#' @param splitPattern a pattern used to extract, from the column names, the
#'   sample field (default " ")
#' @param numCol Once the column names are split by splitPattern, the column
#'   number with the sample name (default 2)
#' @param genePrefix Prefix for the mitochondrial genes (default "^MT-" for
#'   Human, mouse "^mt-")
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 position_jitter
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 ylim
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#'
#' @returns `mitochondrialPercentagePlot()` returns a `list` with:
#'   * "plot" a `violin-boxplot` object
#'   * "sizes" a sizes `data.frame`
#'
#' @export
#'
#' @examples
#' mitPercPlot <- mitochondrialPercentagePlot(objCOTAN)[["plot"]]
#'
#' @rdname RawDataCleaning
#'
mitochondrialPercentagePlot <- function(objCOTAN, splitPattern = " ",
                                        numCol = 2L, genePrefix = "^MT-") {
  sizes <- getCellsSize(objCOTAN)
  sizes <- as.data.frame(sizes)
  sizes <- setColumnInDF(sizes, seq_len(nrow(sizes)), colName = "n")
  {
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

  mitGenes <- getGenes(objCOTAN)[str_detect(getGenes(objCOTAN),
                                            pattern = genePrefix)]
  mitGenesData <- getRawData(objCOTAN)[getGenes(objCOTAN) %in% mitGenes, ]
  if (!identical(colnames(mitGenesData), rownames(sizes))) {
    warning("Problem with cells' order!")
  }
  mitGenesData <- t(mitGenesData)
  sizes <- setColumnInDF(sizes, rowSums(mitGenesData), colName = "sum.mit")
  sizes <- setColumnInDF(sizes, colName = "mit.percentage",
                         colToSet = round(100.0 * sizes[["sum.mit"]] /
                                            sizes[["sizes"]], digits = 2L))

  plot <-
    sizes %>%
    ggplot(aes(x = sample, y = mit.percentage, fill = sample)) +
    #geom_point(size = 0.5) +
    geom_flat_violin(position = position_nudge(x = .15, y = 0.0),
                     adjust = 2.0, alpha = 0.5) +
    geom_point(position = position_jitter(width = .1),
               size = .4, color = "black", alpha = 0.5) +
    geom_boxplot(aes(x = sample, y = mit.percentage, fill = sample),
                 outlier.shape = NA, alpha = .8,
                 width = .15, colour = "gray65", size = 0.6) +
    labs(title = "Mitochondrial percentage of reads",
         y = "% (mit. reads / tot reads * 100)",
         x = "") +
    scale_y_continuous(expand = c(0.0, 0.0)) +
    #ylim(0,max(sizes$sizes)) +
    plotTheme("size-plot")

  return(list("plot" = plot, "sizes" = sizes))
}
