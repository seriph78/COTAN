#' GDIPlot
#'
#' This function directly evaluate and plot the GDI for a sample.
#'
#' @param objCOTAN a `COTAN` object
#' @param genes a named `list` of genes to label. Each array will have different
#'   color.
#' @param cond a string corresponding to the condition/sample (it is used only
#'   for the title).
#' @param statType type of statistic to be used. Default is "S": Pearson's
#'   chi-squared test statistics. "G" is G-test statistics
#' @param GDIThreshold the threshold level that discriminates uniform clusters.
#'   It defaults to \eqn{1.5}
#' @param GDIIn when the GDI data frame was already calculated, it can be put
#'   here to speed up the process (default is `NULL`)
#'
#' @returns A `ggplot2` object
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
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#'
#' @importFrom stats quantile
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)
#' genes <- c("g-000010", "g-000020", "g-000030", "g-000300", "g-000330",
#'            "g-000510", "g-000530", "g-000550", "g-000570", "g-000590")
#' gdiPlot <- GDIPlot(objCOTAN, genes = genes, cond = "raw")
#' plot(gdiPlot)
#'
#' @rdname GDIPlot
#'
GDIPlot <- function(objCOTAN, genes, cond = "",
                    statType = "S", GDIThreshold = 1.5,
                    GDIIn = NULL) {
  logThis("GDI plot", logLevel = 2L)

  if (is_empty(GDIIn)) {
    GDIDf <- calculateGDI(objCOTAN, statType = statType)
  } else {
    GDIDf <- GDIIn
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

  # drop housekeeping genes i.e. those that has GDI <= -5
  {
    genesToKeep <- (GDIDf[["GDI"]] > -5.0)
    logThis(paste("Removed", sum(!genesToKeep), "low GDI genes",
                  "(such as the housekeeping) in GDI plot"), logLevel = 1L)
    GDIDf <- GDIDf[genesToKeep, ]
  }

  qualColPals <- brewer.pal.info[brewer.pal.info[["category"]] == "qual", ]
  colVector <- unlist(mapply(brewer.pal, qualColPals[["maxcolors"]],
                             rownames(qualColPals)))

  myColours <- set_names(colVector[seq_along(names(genes))], names(genes))

  labelledGenes <- GDIDf[["colors"]] != "none"

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
          ggtitle(paste("GDI plot ", cond)) +
          plotTheme("GDI", textSize = 10L)

  return(plot)
}
