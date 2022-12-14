#' GDIPlot
#'
#' This function directly evaluate and plot the GDI for a sample.
#'
#' @param objCOTAN A COTAN object
#' @param genes A named list of genes to label. Each array will have different
#'   color.
#' @param cond A string corresponding to the condition/sample (it is used only
#'   for the title).
#' @param statType Type of statistic to be used. Default is "S": Pearson's
#'   chi-squared test statistics. "G" is G-test statistics
#' @param GDI.df When the GDI data frame was already calculated, it can be put
#'   here to speed up the process. Default is NULL.
#'
#' @returns A ggplot2 object
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
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggrepel geom_label_repel
#'
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' gdiPlot <- GDIPlot(objCOTAN, cond = "raw")
#' plot(gdiPlot)
#'
#' @rdname GDIPlot
#'
GDIPlot <- function(objCOTAN, genes, cond = "",
                    statType = "S", GDI.df = NULL) {
  sum.raw.norm <- NULL

  logThis("GDI plot", logLevel = 2)

  if (is_empty(GDI.df)) {
    GDI <- calculateGDI(objCOTAN, statType = statType)
  } else {
    GDI <- GDI.df
  }

  GDI$colors <- "none"
  for (n in names(genes)) {
    GDI[rownames(GDI) %in% genes[[n]],]$colors <- n
  }

  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  mycolours <- set_names(col_vector[seq_along(names(genes))], names(genes))

  textdf <- GDI[!GDI$colors == "none",]

  GDI_plot <- ggplot(subset(GDI, colors == "none"),
                     aes(x = sum.raw.norm, y = GDI)) +
              geom_point(alpha = 0.3, color = "#8491B4B2", size = 2.5) +
              geom_point(data = subset(GDI, colors != "none"),
                         aes(x = sum.raw.norm, y = GDI, colour = colors),
                         size = 2.5, alpha = 0.8) +
              geom_hline(yintercept = quantile(GDI$GDI)[4], linetype = "dashed", color = "darkblue") +
              geom_hline(yintercept = quantile(GDI$GDI)[3], linetype = "dashed", color = "darkblue") +
              geom_hline(yintercept = 1.5, linetype = "dotted", color = "red", linewidth = 0.5) +
              scale_color_manual("Status", values = mycolours) +
              scale_fill_manual( "Status", values = mycolours) +
              xlab("log normalized counts") +
              ylab("GDI") +
              geom_label_repel(data = textdf,
                               aes(x = sum.raw.norm, y = GDI,
                                   label = rownames(textdf), fill = colors),
                               label.size = NA, max.overlaps = 40, alpha = 0.8,
                               direction = "both", na.rm = TRUE, seed = 1234) +
              ggtitle(paste("GDI plot ", cond)) +
              plotTheme("GDI", textSize = 10)

  return(GDI_plot)

}
