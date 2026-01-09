
## ----- Genes percentage plot -----

#' @details `genesPercentagePlot()` plots the percentage of expression of the
#'   passed-in genes against the overall cell's library size
#'
#' @param objCOTAN a `COTAN` object
#' @param genes an array of gene names
#' @param title The title of the plot; it defaults to `"Genes' percentage of
#'   reads"`
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
#' @returns `genesPercentagePlot()` returns a `list` with:
#'   * `"plot"` a `half-violin-boxplot` object
#'   * `"sizes"` a sizes `data.frame`
#'
#' @export
#'
#' @examples
#' genes <- getGenes(objCOTAN)[1:30]
#' genesPercPlot <-
#'   genesPercentagePlot(objCOTAN, genes = genes)[["plot"]]
#' plot(genesPercPlot)
#'
#' @rdname RawDataCleaning
#'
genesPercentagePlot <- function(objCOTAN,
                                genes,
                                title = "Genes' percentage of reads",
                                condName = "",
                                conditions = NULL) {
  df <- data.frame()
  df <- setColumnInDF(df, colToSet = getCellsSize(objCOTAN), colName = "sizes")
  df <- setColumnInDF(df, colToSet = seq_len(nrow(df)),      colName = "n")
  assert_that(identical(rownames(df), getCells(objCOTAN)))

  c(., conditions) %<-%
    normalizeNameAndLabels(objCOTAN, name = condName,
                           labels = conditions, isCond = TRUE)
  assert_that(!is_empty(conditions),
              identical(rownames(df), names(conditions)))

  df <- setColumnInDF(df, colToSet = conditions, colName = "sample")

  genesInObj <- genes %in% getGenes(objCOTAN)
  genes <- genes[genesInObj]
  if (is_empty(genes)) {
    stop("passed genes are not part of the overall genes' list")
  }
  if (!all(genesInObj)) {
    warning("some of passed genes are not part of the overall genes' list")
  }

  genesData <-
    getRawData(objCOTAN)[getGenes(objCOTAN) %in% genes, , drop = FALSE]
  if (!identical(colnames(genesData), rownames(df))) {
    warning("Problem with cells' order!")
  }
  df <- setColumnInDF(df, colToSet = colSums(genesData), colName = "sum")

  perc <- round(100.0 * df[["sum"]] / df[["sizes"]], digits = 2L)
  df <- setColumnInDF(df, colToSet = perc, colName = "percentage")

  plot <-
    df %>%
    ggplot(
      aes(x = .data$sample, y = .data$percentage, fill = .data$sample)
    ) +
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
      aes(x = .data$sample, y = .data$percentage, fill = .data$sample),
      outlier.shape = NA, alpha = 0.8,
      width = 0.15, colour = "gray65", size = 0.6
    ) +
    labs(
      title = title,
      y     = "% (genes reads / tot reads * 100)",
      x     = ""
    ) +
    scale_y_continuous(expand = c(0.0, 0.0)) +
    coord_cartesian(ylim = c(0.0, max(df[["percentage"]]))) +
    plotTheme("size-plot")

  return(list("plot" = plot, "sizes" = df))
}



## ----- Mitochondrial percentage plot -----

#' @details `mitochondrialPercentagePlot()` plots the percentage of expression
#'   of the mitochondrial genes against the overall cell's library size
#'
#' @param objCOTAN a `COTAN` object
#' @param genePrefix Prefix for the mitochondrial genes; use `"^MT-"` for
#'   Human (default) and `"^mt-"` for Mouse
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#'
#' @importFrom assertthat assert_that
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
mitochondrialPercentagePlot <- function(objCOTAN,
                                        genePrefix = "^MT-",
                                        condName = "",
                                        conditions = NULL) {
  isMitoch <- str_detect(getGenes(objCOTAN), pattern = genePrefix)
  mitGenes <- getGenes(objCOTAN)[isMitoch]
  if (is_empty(mitGenes)) {
    stop("gene prefix resulted in no matches")
  }

  title <- "Mitochondrial percentage of reads"

  res <- genesPercentagePlot(objCOTAN, genes = mitGenes, title = title,
                             condName = condName, conditions = conditions)
  assert_that(identical(colnames(res[["sizes"]]),
                        c("sizes", "n", "sample", "sum", "percentage")))

  # ensures backward compatibility
  colnames(res[["sizes"]]) <-
    c("sizes", "n", "sample", "sum.mit", "mit.percentage")

  return(res)
}
