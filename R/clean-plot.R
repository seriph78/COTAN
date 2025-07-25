
## ------- Clean Plots --------

#' @details `cleanPlots()` creates the plots associated to the output of the
#'   [clean()] method.
#'
#' @param objCOTAN a `COTAN` object
#' @param includePCA a `Boolean` flag to determine whether to calculate the
#'   *PCA* associated with the normalized matrix. When `TRUE` the first four
#'   elements of the returned list will be `NULL`
#'
#' @returns `cleanPlots()` returns a `list` of `ggplot2` plots:
#'   * `"pcaCells"` is for `PCA` cells
#'   * `"pcaCellsData"` is the data of the `PCA` cells (can be plotted)
#'   * `"genes"` is for `B` group cells' genes
#'   * `"UDE"` is for cells' `UDE` against their `PCA`
#'   * `"nu"` is for cell `nu`
#'   * `"zoomedNu"` is the same but zoomed on the left and with an estimate
#'   for the low `nu` threshold that defines problematic cells
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#'
#' @importFrom stats hclust
#' @importFrom stats cutree
#'
#' @importFrom utils head
#'
#' @importFrom Matrix t
#' @importFrom Matrix rowMeans
#'
#' @importFrom BiocSingular runPCA
#' @importFrom BiocSingular IrlbaParam
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
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_label
#' @importFrom ggrepel geom_text_repel
#'
#' @importFrom graphics par
#' @importFrom graphics title
#'
#' @importFrom dplyr filter
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @examples
#' # This creates many infomative plots useful to determine whether
#' # there is still something to drop...
#' # Here we use the tuple-like assignment feature of the `zeallot` package
#' c(pcaCellsPlot, ., genesPlot, UDEPlot, ., zNuPlot) %<-% cleanPlots(objCOTAN)
#' plot(pcaCellsPlot)
#' plot(UDEPlot)
#' plot(zNuPlot)
#'
#' @rdname RawDataCleaning
#'
cleanPlots <- function(objCOTAN, includePCA = TRUE) {
  if (isTRUE(includePCA)) {
    logThis("PCA: START", logLevel = 2L)

    nuNormData <- getNuNormData(objCOTAN)

    logThis("Elaborating PCA - START", logLevel = 3L)

    # re-scale so that all the genes have mean 0.0 and stdev 1.0
    cellsRDM <- runPCA(x = t(nuNormData),
                       rank = 5L,
                       BSPARAM = IrlbaParam(),
                       get.rotation = FALSE)[["x"]]

    logThis("Elaborating PCA - DONE", logLevel = 3L)

    assert_that(identical(rownames(cellsRDM), getCells(objCOTAN)),
                msg = "Issues with pca output")

    distCells <- scale(cellsRDM)
    distCells <- calcDist(distCells, method = "euclidean") # mahlanobis

    cellsRDM <- as.data.frame(cellsRDM)

    logThis("PCA: DONE", logLevel = 2L)

    logThis("Hierarchical clustering: START", logLevel = 2L)

    # hclust cannot operate on more than 2^16 elements
    if (getNumCells(objCOTAN) <= 65500L) {
      hcCells <- hclust(distCells, method = "complete")
      groups <- cutree(hcCells, k = 2L)
    } else {
      groups <- set_names(rep(1L, times = getNumCells(objCOTAN)),
                          getCells(objCOTAN))
      # ensure B group is not empty picking the first cell
      groups[[1L]] <- 2L

      warning("cleanPlots() - More than 65500 cells in the COTAN object: ",
              "'B' group cannot be established and ",
              "is defaulted to include only the first cell", call. = FALSE)
      logThis("Too many cells: cannot establish 'B' group", logLevel = 3L)
    }
    rm(distCells)

    pos1 <- which(groups == 1L)
    pos2 <- which(groups == 2L)

    # ensure "A" group is the most numerous
    if (length(pos1) < length(pos2)) {
      # swap pos1 <-> pos2
      pos3 <- pos1
      pos1 <- pos2
      pos2 <- pos3
    }

    groups[pos1] <- "A"
    groups[pos2] <- "B"

    cellsRDM <- cbind(cellsRDM, groups)

    logThis("Hierarchical clustering: DONE", logLevel = 2L)
    gc()

    toClust <- as.matrix(round(nuNormData, digits = 4L))

    assert_that(identical(colnames(toClust), names(groups)),
                msg = "Error in the cell names")

    # ---- next: to check which genes are specific for the B group of cells
    B <- set_names(as.data.frame(toClust[, pos2]), colnames(toClust)[pos2])

    rm(toClust)

    D <- data.frame("means" = rowMeans(B), "n" = NA)
    rownames(D) <- rownames(B)
    D <- rownames_to_column(D)
    rm(B)

    D <- D[order(D[["means"]], decreasing = TRUE), ]
    rownames(D) <- NULL
    D <- column_to_rownames(D)

    D <- D[D[["means"]] > 0.0, ]
    D[["n"]] <- seq_along(D[["means"]])

    gc()

    #check if the PCA plot is clean enough and from the printed genes,
    #if the smallest group of cells are characterized by particular genes

    cellsRDM_A <- filter(cellsRDM, .data$groups == "A")
    cellsRDM_B <- filter(cellsRDM, .data$groups != "A")

    cellsRDMPlot <- ggplot(cellsRDM_A,
                           aes(x = .data$PC1, y = .data$PC2,
                               colour = .data$groups)) +
                    geom_point(alpha = 0.5, size = 3L) +
                    geom_point(data = cellsRDM_B,
                               aes(x = .data$PC1, y = .data$PC2,
                                   colour = .data$groups),
                               alpha = 0.8, size = 3L) +
                    scale_color_manual(name = "groups",
                                       values = c("A" = "#8491B4B2",
                                                  "B" = "#E64B35FF")) +
                    plotTheme("pca")

    minN <- min(D[["n"]]) + 15.0
    lowD <- D[["n"]] < minN
    genesPlot <- ggplot(D, aes(x = .data$n, y = .data$means)) +
                 geom_point() +
                 ggtitle(label = "B cell group genes mean expression"
                         #subtitle = " - B group NOT removed -"
                         ) +
                 geom_label(data = subset(D, lowD),
                       aes(x = .data$n, y = .data$means,
                           label = rownames(D)[lowD]),
                       #nudge_y = 0.05,
                       nudge_x = 400.0
                       #direction = "x",
                       #angle = 90.0,
                       #vjust = 0L,
                       #segment.size = 0.2
                       ) +
                 plotTheme("genes") +
    theme(plot.title = element_text(hjust = 1.0),
          plot.subtitle = element_text(hjust = 0.95, vjust = -25.0))

    nuEst <- round(getNu(objCOTAN), digits = 7L)
    UDEPlot <- ggplot(cellsRDM, aes(x = .data$PC1, y = .data$PC2,
                                    colour = log(nuEst))) +
               geom_point(size = 1L, alpha = 0.8) +
               scale_color_gradient2(low = "#E64B35B2", mid = "#4DBBD5B2",
                                     high = "#3C5488B2",
                                     midpoint = log(mean(nuEst)),
                                     name = "ln(nu)") +
               ggtitle("Cells PCA coloured by cells efficiency") +
               plotTheme("UDE")
  } else {
    cellsRDMPlot <- NULL
    cellsRDM <- NULL
    genesPlot <- NULL
    UDEPlot <- NULL
  }

  nuDf <- data.frame("nu" = sort(getNu(objCOTAN)),
                     "n" = seq_len(getNumCells(objCOTAN)))
  nuPlot <- ggplot(nuDf, aes(x = .data$n, y = .data$nu)) +
            geom_point(colour = "#8491B4B2", size = 1L) +
            plotTheme("common")

  nuDf <- nuDf[seq_len(min(400L, nrow(nuDf))), ]
  zNuPlot <- ggplot(nuDf, aes(x = .data$n, y = .data$nu)) +
             geom_point(colour = "#8491B4B2", size = 1L) +
             plotTheme("common") +
             xlim(0L, nrow(nuDf)) +
             ylim(0.0, round(max(nuDf[["nu"]]) + 0.05, digits = 2L))

  # estimate the elbow point if any...
  secondDer <- diff(nuDf[seq_len(min(100L, nrow(nuDf))), "nu"],
                    differences = 2L)
  if (min(secondDer) < -0.01) {
    lowUDEThr <- nuDf[(max(which(secondDer < -0.01)) + 1L), "nu"] - 0.005
    zNuPlot <- zNuPlot +
               geom_hline(yintercept = lowUDEThr, linetype = "dashed",
                          color = "darkred") +
               annotate(geom = "text",
                        x = (nrow(nuDf) / 2.0),
                        y = (lowUDEThr - 0.05),
                        label = paste0("cells with nu < ", lowUDEThr,
                                       " should probably be removed"),
                        color = "darkred", size = 4.5)
  }

  return(list("pcaCells" = cellsRDMPlot, "pcaCellsData" = cellsRDM,
              "genes" = genesPlot, "UDE" = UDEPlot,
              "nu" = nuPlot, "zoomedNu" = zNuPlot))
}

## ------- Scree Plots --------
##
#' @details `screePlot()` creates a plots showing the explained variance of the
#'   components of a PCA
#'
#' @param pcaStdDev a `vector` with the standard deviations of the various
#'  components
#'
#' @returns `screePlot()` returns a `ggplot2` plot for the explained variances
#'
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 labs
#'
#' @rdname RawDataCleaning
#'

screePlot <- function(pcaStdDev) {
  # Extract variance explained
  varExplained <- pcaStdDev^2L / sum(pcaStdDev^2L)

  # Convert to data frame for ggplot
  screeDF <- data.frame(PC = seq_along(varExplained), Variance = varExplained)

  # Create ggplot scree plot
  ggplot(screeDF, aes(x = .data$PC, y = .data$Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_point(color = "red", size = 3L) +
    geom_line(group = 1L, color = "red") +
    labs(title = "Scree Plot", x = "Principal Component",
         y = "Explained Variance") +
    plotTheme("pca")
}

