
## ----- Cells size plot -----

#' @details `cellSizePlot()` plots the raw library size for each cell and
#'   sample.
#'
#' @param objCOTAN a `COTAN` object
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#'
#' @returns `cellSizePlot()` returns a `half-violin-boxplot` object
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
#' @export
#'
#' @examples
#' lsPlot <- cellSizePlot(objCOTAN)
#' plot(lsPlot)
#'
#' @rdname RawDataCleaning
#'
cellSizePlot <- function(objCOTAN, condName = "", conditions = NULL) {
  sizes <- sort(getCellsSize(objCOTAN))
  sizes <- as.data.frame(sizes)

  sizes <- setColumnInDF(sizes, colToSet = seq_len(nrow(sizes)), colName = "n")

  c(., conditions) %<-%
    normalizeNameAndLabels(objCOTAN, name = condName,
                           labels = conditions, isCond = TRUE)
  assert_that(!is_empty(conditions))

  sizes <- setColumnInDF(sizes, conditions[rownames(sizes)], colName = "sample")

  plot <-
    sizes %>%
    ggplot(aes(x = sample, y = sizes, fill = sample)) +
    stat_slabinterval(
      aes(thickness = after_stat(pdf)),
      side = "right",                   # draw only the right half
      show_point = FALSE,               # slab only (no dot)
      show_interval = FALSE,            # slab only (no interval)
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
      aes(x = sample, y = sizes, fill = sample),
      outlier.shape = NA, alpha = 0.8, width = 0.15,
      colour = "gray65", size = 0.6
    ) +
    labs(
      title = "Cell library size",
      y = "Size (read number)",
      x = ""
    ) +
    coord_cartesian(ylim = c(0.0, max(sizes[["sizes"]]))) +
    plotTheme("size-plot")

  return(plot)
}


## ----- Genes size plot -----

#' @details `genesSizePlot()` plots the raw gene number (reads > 0) for each
#'   cell and sample
#'
#' @param objCOTAN a `COTAN` object
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#'
#' @returns `genesSizePlot()` returns a `half-violin-boxplot` object
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 position_nudge
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 position_jitter
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 after_stat
#'
#' @importFrom ggdist stat_slabinterval
#'
#' @export
#'
#' @examples
#' gsPlot <- genesSizePlot(objCOTAN)
#' plot(gsPlot)
#'
#' @rdname RawDataCleaning
#'
genesSizePlot <- function(objCOTAN, condName = "", conditions = NULL) {
  sizes <- sort(getNumExpressedGenes(objCOTAN))
  sizes <- as.data.frame(sizes)

  sizes <- setColumnInDF(sizes, seq_len(nrow(sizes)), colName = "n")

  c(., conditions) %<-%
    normalizeNameAndLabels(objCOTAN, name = condName,
                           labels = conditions, isCond = TRUE)
  assert_that(!is_empty(conditions))

  sizes <- setColumnInDF(sizes, conditions[rownames(sizes)], colName = "sample")

  plot <-
    sizes %>%
    ggplot(aes(x = sample, y = sizes, fill = sample)) +
    stat_slabinterval(
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
      aes(x = sample, y = sizes, fill = sample),
      outlier.shape = NA, alpha = 0.8, width = 0.15,
      colour = "gray65", size = 0.6
    ) +
    labs(
      title = "Detected gene number",
      y = "Size (number of genes)",
      x = ""
    ) +
    coord_cartesian(ylim = c(0.0, max(sizes[["sizes"]]))) +
    plotTheme("size-plot")

  return(plot)
}
