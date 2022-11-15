#' clean
#'
#' Main function that can be used to check and clean the dataset.
#' It also produces and stores the estimators for nu and lambda and fills
#' the raw.norm (raw / nu)
#'
#' @param objCOTAN COTAN object
#' @return lists of objects...
#' "object" is the COTAN object,
#' "data" a list with:
#'        "cluster1" is the first  cell cluster,
#'        "cluster2" is the second cell cluster,
#'        "D" is the B cells' group genes,
#'        "pcaCells" pca numeric data,
#' "plots" a list with ggplot2 plots:
#'         "pcaCells" is for pca cells,
#'         "genes" is for B cells' group genes,
#'         "UDE" is for cell UDE,
#'         "nu" is for cell nu.
#'
#' @export
#'
#' @import gsubfn
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
#' @importFrom Matrix colMeans
#'
#' @examples
#' data("ERCC.cotan")
#' list[objCOTAN, data, plots] <- clean(ERCC.cotan)
#' plot(plots[["UDEPlot"]])
#'
#' @rdname clean
setMethod(
  "clean",
  "COTAN",
  function(objCOTAN) {

    # We want to discard genes having less than 3 non-zero counts per 1000 cells
    threshold <- round(getNumCells(objCOTAN) * 3 / 1000, digits = 0)
    genesToDrop <- getGenes(objCOTAN)[rowSums(getZeroOneProj(objCOTAN)) <= threshold]

    objCOTAN <- dropGenesCells(objCOTAN, genes = genesToDrop)
    print(paste0("Genes selection done:",
                 " dropped [", length(genesToDrop), "] genes"))

    gc()

    list[distCells, pcaCells, objCOTAN] = runEstimatesLinear(objCOTAN)

    gc()

    print("Hierarchical clustering: START")

    hcCells <- hclust(distCells, method = "complete")
    rm(distCells)
    gc()

    groups <- cutree(hcCells, k = 2)

    pos1 <- which(groups == 1)
    pos2 <- which(groups == 2)

    # ensure "A" group is the most numerous
    if( length(pos1) < length(pos2) ) {
      # swap pos1 <-> pos2
      pos3 <- pos1
      pos1 <- pos2
      pos2 <- pos3
    }

    cluster1 <- names(pos1)
    cluster2 <- names(pos2)

    groups[pos1] <- "A"
    groups[pos2] <- "B"

    print("Hierarchical clustering: DONE")

    toClust <- as.matrix(round(getNormalizedData(objCOTAN), digits = 4))

    if (!identical(colnames(toClust), names(groups))) {
      stop("Error in the cell names")
    }

    # ---- next: to check which genes are specific for the B group of cells
    {
      B <- as.data.frame(toClust[, pos2])
      colnames(B) <- colnames(toClust)[pos2]
      #print(utils::head(B, 15))

      rm(toClust)
      gc()

      D <- data.frame("means" = rowMeans(B), "n" = NA)
      rownames(D) <- rownames(B)
      D <- rownames_to_column(D)
      D <- D[order(D[["means"]], decreasing = TRUE), ]
      rownames(D) <- c()
      D <- column_to_rownames(D)

       D <- D[D[["means"]] > 0,]
       D[["n"]] <- seq_along(D[["means"]])
    }

    #check if the pca plot is clean enought and from the printed genes,
    #if the smalest group of cells are caratterised by particular genes

    pcaCells <- cbind(pcaCells, groups)

    PC2 <- PC1 <- NULL

    pcaCellsPlot <- ggplot(subset(pcaCells, groups == "A"),
                           aes(x = PC1, y = PC2, colour = groups)) +
                    geom_point(alpha = 0.5, size = 3) +
                    geom_point(data = subset(pcaCells, groups != "A" ),
                               aes(x = PC1, y = PC2, colour = groups),
                               alpha = 0.8, size = 3) +
                    scale_color_manual("groups", values = c("A" = "#8491B4B2", "B"="#E64B35FF")) +
                    plotTheme("pca")

    # genes plot
    minN <- max(D[["n"]])- 15
    genesPlot    <- ggplot(D, aes(x = n, y = means)) +
                    ggtitle(label = "B cell group genes mean expression",
                            subtitle = " - B group NOT removed -") +
                    geom_text_repel(data = subset(D, n > minN),
                          aes(n, means, label = rownames(D[D[["n"]] > minN,])),
                          nudge_y = 0.05, nudge_x = 0.05, direction = "x",
                          angle = 90, vjust = 0, segment.size = 0.2) +
                    plotTheme("genes")

    # UDE/nu plot
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

    return( list("objCOTAN" = objCOTAN,
                 "data"     = list("cluster1" = cluster1, "cluster2" = cluster2,
                                   "D" = D, "pcaCells" = pcaCells),
                 "plots"    = list("pcaCells" = pcaCellsPlot, "genes" = genesPlot,
                                   "UDE" = UDEPlot, "nu" = nuPlot)) )
  }
)
