
#' @details `metaGenesPlot()` shows the overlap of the points in the
#'   *meta-genes* in the given \eqn{(x, y)} space (usually the PCA).
#'   The union of all names in the *meta-genes* must be a subset of the names in
#'   the coordinates
#'
#' @param x an `array` giving the first coordinate of the points
#' @param y an `array` giving the second coordinate of the points
#' @param metaGenes a `list` of *meta-genes*, such as the result of
#'   [defineMetaGenes()].
#'
#' @returns `metaGenesPlot()` returns a `list` with:
#'  * "plot" the wanted `ggplot`
#'  * "overlapping" an `array` that for each point gives how many
#'    *meta-genes* it belongs to; centers are returned with the highest value
#'
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_gradient
#'
#' @importFrom rlang set_names
#'
#' @rdname MetaGenes
#'

metaGenesPlot <- function(x, y, metaGenes) {
  pointNames <- names(x)
  colorsValue <- set_names(rep(0, length(pointNames)), pointNames)

  for (center in names(metaGenes)) {
    metaGene <- metaGenes[[center]]
    colorsValue[metaGene] <- colorsValue[metaGene] + 1
  }
  colorsValue[names(metaGenes)] <- max(colorsValue) + 1

  plotDF <- data.frame(x = x, y = y, color = colorsValue)
  myPlot <- ggplot(plotDF) +
    geom_point(aes(x = x, y = y, color = color, size = 0.5)) +
    scale_colour_gradient(low = "grey100", high = "red")

  return(list("plot" = myPlot, "overlapping" = colorsValue))
}
