# --------------------- Uniform Clusters ----------------------

#' @name UniformClusters
#'
#' @details `GDIPlot()` directly evaluates and plots the `GDI` for a sample.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes a named `list` of genes to label. Each array will have different
#'   color.
#' @param condition a string corresponding to the condition/sample (it is used
#'   only for the title).
#' @param statType type of statistic to be used. Default is "S": Pearson's
#'   chi-squared test statistics. "G" is G-test statistics
#' @param GDIThreshold the threshold level that discriminates uniform clusters.
#'   It defaults to \eqn{1.43}
#' @param GDIIn when the `GDI` data frame was already calculated, it can be put
#'   here to speed up the process (default is `NULL`)
#'
#' @returns `GDIPlot()` returns a `ggplot2` object
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggrepel geom_label_repel
#'
#' @importFrom stats quantile
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname UniformClusters
#'
GDIPlot <- function(objCOTAN, genes, condition = "",
                    statType = "S", GDIThreshold = 1.43,
                    GDIIn = NULL) {
  logThis("GDI plot", logLevel = 2L)

  if (!is_empty(GDIIn)) {
    GDIDf <- GDIIn
  } else {
    gdi <- getGDI(objCOTAN)
    if (!is_empty(gdi) && statType == "S") {
      # complete the GDIDf
      GDIDf <- setColumnInDF(GDIDf, colToSet = gdi, colName = "GDI",
                             rowNames = getGenes(objCOTAN))
      sumRawNorm <- log(rowSums(getNormalizedData(objCOTAN, retLog = FALSE)))
      GDIDf <- setColumnInDF(GDIDf, colToSet = sumRawNorm,
                             colName = "sum.raw.norm")
    } else {
      GDIDf <- calculateGDI(objCOTAN, statType = statType)
    }
  }

  GDIDf[["colors"]] <- "none"
  for (n in names(genes)) {
    genesSelec <- rownames(GDIDf) %in% genes[[n]]
    if (any(genesSelec)) {
      GDIDf[genesSelec, "colors"] <- n
    } else {
      logThis(paste("GDIPlot - none of the genes in group", n,
                    "is present: will be ignored"), logLevel = 1L)
    }
  }

  # drop fully-expressed genes i.e. those that have GDI <= -5.0
  {
    genesToKeep <- (GDIDf[["GDI"]] > -5.0)
    logThis(paste("Removed", sum(!genesToKeep), "low GDI genes",
                  "(such as the fully-expressed) in GDI plot"), logLevel = 1L)
    GDIDf <- GDIDf[genesToKeep, ]
  }

  myColours <- set_names(getColorsVector(length(names(genes))), names(genes))

  labelledGenes <- GDIDf[["colors"]] != "none"

  if (isEmptyName(condition)) {
    condition <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])
  }
  title <- paste0("GDI plot - ", condition)

  plot <- ggplot(subset(GDIDf, colors == "none"),
                 aes(x = sum.raw.norm, y = GDI)) +
          geom_point(alpha = 0.3, color = "#8491B4B2", size = 2.5) +
          geom_point(data = subset(GDIDf, colors != "none"),
                     aes(x = sum.raw.norm, y = GDI, colour = colors),
                     size = 2.5, alpha = 0.8) +
          geom_hline(yintercept = quantile(GDIDf[["GDI"]])[[4L]],
                     linetype = "dashed", color = "darkblue") +
          geom_hline(yintercept = quantile(GDIDf[["GDI"]])[[3L]],
                     linetype = "dashed", color = "darkblue") +
          geom_hline(yintercept = GDIThreshold, linetype = "dotted",
                     color = "red", linewidth = 0.5) +
          scale_color_manual("Status", values = myColours) +
          scale_fill_manual( "Status", values = myColours) +
          xlab("log normalized counts") +
          ylab("GDI") +
          geom_label_repel(data = GDIDf[labelledGenes, ],
                           aes(x = sum.raw.norm, y = GDI, fill = colors),
                           label = rownames(GDIDf)[labelledGenes],
                           label.size = NA, max.overlaps = 40L, alpha = 0.8,
                           direction = "both", na.rm = TRUE, seed = 1234L) +
          ggtitle(title) +
          plotTheme("GDI", textSize = 10L)

  return(plot)
}
