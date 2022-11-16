#' cleanPlots
#'
#' Retrieves the plots associated to the output of the clena method.
#'
#' @param objCOTAN COTAN object
#' @param pcaCells pca numeric data
#' @param D B cells' group genes
#'
#' @return lists of ggplot2 plots:
#' "pcaCells" is for pca cells,
#' "genes" is for B cells' group genes,
#' "UDE" is for cell UDE,
#' "nu" is for cell nu.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' data("ERCC.cotan")
#' list[objCOTAN, data] <- clean(ERCC.cotan)
#' plots <- cleanPlots(objCOTAN, data[["pcaCells"]], data[["D"]])
#' plot(plots[["UDEPlot"]])
#'
#' @rdname cleanPlots
cleanPlots <- function(objCOTAN, pcaCells, D) {

  #check if the pca plot is clean enough and from the printed genes,
  #if the smallest group of cells are characterized by particular genes

  PC2 <- PC1 <- NULL

  groups <- pcaCells[["groups"]]
  pcaCellsPlot <- ggplot(subset(pcaCells, groups == "A"),
                         aes(x = PC1, y = PC2, colour = groups)) +
                  geom_point(alpha = 0.5, size = 3) +
                  geom_point(data = subset(pcaCells, groups != "A" ),
                             aes(x = PC1, y = PC2, colour = groups),
                             alpha = 0.8, size = 3) +
                  scale_color_manual("groups", values = c("A" = "#8491B4B2", "B"="#E64B35FF")) +
                  plotTheme("pca")

  minN <- max(D[["n"]])- 15
  genesPlot    <- ggplot(D, aes(x = n, y = means)) +
                  ggtitle(label = "B cell group genes mean expression",
                          subtitle = " - B group NOT removed -") +
                  geom_text_repel(data = subset(D, n > minN),
                        aes(n, means, label = rownames(D[D[["n"]] > minN,])),
                        nudge_y = 0.05, nudge_x = 0.05, direction = "x",
                        angle = 90, vjust = 0, segment.size = 0.2) +
                  plotTheme("genes")


  nuEst = round(getNu(objCOTAN), digits = 7)
  UDEPlot      <- ggplot(pcaCells, aes(x = PC1,y = PC2, colour = log(nuEst))) +
                  geom_point(size = 1, alpha = 0.8) +
                  scale_color_gradient2(low = "#E64B35B2", mid = "#4DBBD5B2", high = "#3C5488B2",
                                        midpoint = log(mean(nuEst)), name = "ln(nu)") +
                  ggtitle("Cells PCA coloured by cells efficiency") +
                  plotTheme("UDE")

  nuDf = data.frame("nu" = sort(getNu(objCOTAN)), "n" = c(1:getNumCells(objCOTAN)))
  nuPlot       <- ggplot(nuDf, aes(x = n, y = nu)) +
                  geom_point(colour = "#8491B4B2", size = 1) +
                  plotTheme("common")

  return( list("pcaCells" = pcaCellsPlot, "genes" = genesPlot,
               "UDE" = UDEPlot, "nu" = nuPlot) )
}

