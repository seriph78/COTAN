
## ----- Mitochondrial percentage plot -----

#' @details `mitochondrialPercentagePlot()` plots the raw library size for each
#'   cell and sample.
#'
#' @param objCOTAN a `COTAN` object
#' @param genePrefix Prefix for the mitochondrial genes (default "^MT-" for
#'   Human, mouse "^mt-")
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`

#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 position_jitter
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 after_stat
#'
#' @importFrom ggdist stat_slabinterval
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom rlang is_empty
#'
#' @importFrom stringr str_detect
#'
#' @returns `mitochondrialPercentagePlot()` returns a `list` with:
#'   * `"plot"` a `half-violin-boxplot` object
#'   * `"sizes"` a sizes `data.frame`
#'
#' @export
#'
#' @examples
#' mitPercPlot <-
#'   mitochondrialPercentagePlot(objCOTAN, genePrefix = "g-0000")[["plot"]]
#' plot(mitPercPlot)
#'
#' @rdname RawDataCleaning
#'
mitochondrialPercentagePlot <- function(objCOTAN, genePrefix = "^MT-",
                                        condName = "", conditions = NULL) {
  sizes <- getCellsSize(objCOTAN)
  sizes <- as.data.frame(sizes)

  sizes <- setColumnInDF(sizes, seq_len(nrow(sizes)), colName = "n")

  c(., conditions) %<-%
    normalizeNameAndLabels(objCOTAN, name = condName,
                           labels = conditions, isCond = TRUE)
  assert_that(!is_empty(conditions))

  sizes <- setColumnInDF(sizes, conditions[rownames(sizes)], colName = "sample")

  mitGenes <- getGenes(objCOTAN)[str_detect(getGenes(objCOTAN),
                                            pattern = genePrefix)]
  if (is_empty(mitGenes)) {
    stop("gene prefix resulted in no matches")
  }

  mitGenesData <-
    getRawData(objCOTAN)[getGenes(objCOTAN) %in% mitGenes, , drop = FALSE]
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
    ggplot(aes(x = .data$sample, y = .data$mit.percentage,
               fill = .data$sample)) +
    ggdist::stat_slabinterval(
      aes(thickness = after_stat(pdf)),
      side = "right",                 # half-violin to the right
      show_point = FALSE,             # slab only
      show_interval = FALSE,          # no intervals
      position = position_nudge(x = 0.15),
      normalize = "groups",
      adjust = 2.0,
      trim = TRUE,
      alpha = 0.5,
      slab_colour = "black",      # outline color
      slab_linewidth = 0.6,       # outline width (ggplot2 ≥ 3.4)
      slab_alpha = 0.5            # fill alpha only; leaves stroke opaque
    ) +
    geom_point(
      position = position_jitter(width = 0.1),
      size = 0.4, color = "black", alpha = 0.5
    ) +
    geom_boxplot(
      aes(x = .data$sample, y = .data$mit.percentage, fill = .data$sample),
      outlier.shape = NA, alpha = 0.8,
      width = 0.15, colour = "gray65", size = 0.6
    ) +
    labs(
      title = "Mitochondrial percentage of reads",
      y = "% (mit. reads / tot reads * 100)",
      x = ""
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.0, max(sizes[["mit.percentage"]]))) +
    plotTheme("size-plot")

  return(list("plot" = plot, "sizes" = sizes))
}
