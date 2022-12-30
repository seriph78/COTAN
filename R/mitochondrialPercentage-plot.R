
#' mitochondrialPercentagePlot
#'
#' @description Function that plots the raw library size for each cell and
#'   sample.
#'
#' @param objCOTAN a `COTAN` object
#' @param splitPattern Pattern used to extract, from the column names, the
#'   sample field (default " ")
#' @param numCols Once the column names are split by splitPattern, the column
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
#'
#' @returns a list with the violin-boxplot plot and a sizes df
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' mpPlot <- mitochondrialPercentagePlot(objCOTAN)[["plot"]]
#' plot(mpPlot)
#'
#' @rdname mitochondrialPercentagePlot
#'
mitochondrialPercentagePlot <- function(objCOTAN, splitPattern = " ",
                                        numCols = 2, genePrefix = "^MT-"){
  sizes <- getCellsSize(objCOTAN)
  sizes <- as.data.frame(sizes)
  sizes$n <- c(1:dim(sizes)[1])
  sizes$sample <- stringr::str_split(rownames(sizes),
                                     pattern = splitPattern,simplify = T)[, numCols]

  mitGenes <- getGenes(obj)[stringr::str_detect(getGenes(obj), pattern = genePrefix)]
  mitGenesData <- as.data.frame(getRawData(obj)[getGenes(obj) %in% mitGenes, ])
  if(!identical(colnames(mitGenesData), rownames(sizes))){
    warning("Problem with cells' order!")
  }
  mitGenesData <- t(mitGenesData)
  sizes$sum.mit <- Matrix::rowSums(mitGenesData)
  sizes$mit.percentage <- round(sizes$sum.mit / sizes$sizes * 100, digits = 2)

  plot <-
    sizes %>%
    ggplot(aes(x = sample, y = mit.percentage, fill = sample)) +
    #geom_point(size = 0.5) +
    geom_flat_violin(position = position_nudge(x = .15, y = 0),adjust =2, alpha = 0.5) +
    geom_point(position = position_jitter(width = .1), size = .4, color="black", alpha= 0.5) +
    geom_boxplot(aes(x = sample, y = mit.percentage, fill = sample),
                 outlier.shape = NA, alpha = .8, width = .15, colour = "gray65", size = 0.6)+
    labs(title = "Mitochondrial percentage of reads",
         y = "% (mit. reads / tot reads * 100)",
         x = "") +
    scale_y_continuous(expand = c(0, 0)) +
    #ylim(0,max(sizes$sizes))+
    plotTheme("size-plot")

  return(list("plot" = plot, "sizes" = sizes))
}



