#' clean
#'
#' Main function that can be used to check and clean the dataset.
#' It also produces and stores the estimators for nu and lambda and fills
#' the raw.norm (raw / nu)
#'
#' @param objCOTAN COTAN object
#' @return a list of objects containing:
#' "object" is the COTAN object,
#' "cl1" is the first cell cluster,
#' "cl2" is the second cell cluster,
#' "D" is the B cells' group genes,
#' "pca_cells" pca numeric data.
#' "pca.cell.2" is a ggplot2 cell pca plot,
#' "genes.plot" is a ggplot2 B cells' group genes plot,
#' "UDE.plot" is a ggplot2 cell UDE plot,
#' @export
#'
#' @import dplyr
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
#' ttm <- clean(ERCC.cotan)
#'
#' @rdname clean
setMethod(
  "clean",
  "COTAN",
  function(objCOTAN) {

    print("Starting")

    # We want to discard genes having less than 3 non-zero counts per 1000 cells
    threshold <- round(getNumCells(objCOTAN) * 3 / 1000, digits = 0)
    genesToDrop <- getGenes(objCOTAN)[rowSums(getZeroOneProj(objCOTAN)) <= threshold]

    objCOTAN <- dropGenesCells(objCOTAN, genes = genesToDrop)
    print(paste0("Genes selection done:",
                 " dropped [", length(genesToDrop), "] genes"))

    gc()

    list[dist_cells, pca_cells, objCOTAN] = runEstimatesLinear(objCOTAN)
    print("Linear estimations: DONE")

    gc()

    # ** !! **
    print("starting hclust")

    hc_cells <- hclust(dist_cells, method = "complete")
    gc()

    groups <- cutree(hc_cells, k=2)

    if(length(groups[groups == 1]) < length(groups[groups == 2])  ){
      groups[groups == 1] <- "B"
      groups[groups == 2] <- "A"
    }
    else{
      groups[groups == 1] <- "A"
      groups[groups == 2] <- "B"
    }

    cl2 <- names(which(cutree(hc_cells, k = 2) == 2))
    cl1 <- names(which(cutree(hc_cells, k = 2) == 1))

    if (length(cl2) > length(cl1) ) {
      cl2  <-  names(which(cutree(hc_cells, k = 2) == 1))
      cl1  <-  names(which(cutree(hc_cells, k = 2) == 2))
    }

    to_clust <- getNormalizedData(objCOTAN)
    t_to_clust <- Matrix::t(to_clust)
    t_to_clust <- round(t_to_clust,digits = 4)
    t_to_clust <- as.data.frame(as.matrix(t_to_clust))

    if (identical(rownames(t_to_clust),names(groups))) {
      t_to_clust <- cbind(t_to_clust,groups)
    }
    else {
      stop("Error in the cell names")
    }

    # ---- next: to check which genes are specific for the B group of cells

    to_clust <- as.matrix(round(to_clust,digits = 4))

    B <- as.data.frame(to_clust[,colnames(to_clust) %in% cl2])
    rm(to_clust)
    gc()

    colnames(B)<-cl2
    B <- rownames_to_column(B)
    if (dim(B)[2]>2) {
      B <- B[order(rowMeans(B[,2:length(colnames(B))]),
                   decreasing = TRUE), ]
    }
    else{
      B <- B[order(B[,2],decreasing = TRUE), ]
    }

    #print(utils::head(B, 15))

    C <-  arrange(B,rowMeans(B[2:length(colnames(B))]))
    rownames(C) <-  C$rowname
    D <-  data.frame("means" = rowMeans(C[2:length(colnames(C))]),
                     "n" = NA )
    D <- D[D$means>0,]
    D$n <- seq_along(D$means)

    #check if the pca plot is clean enought and from the printed genes,
    #if the smalest group of cells are caratterised by particular genes

    pca_cells  <- cbind(pca_cells, "groups" = t_to_clust$groups)

    PC2 <- PC1 <- NULL

    pca.cell.1 <- ggplot(subset(pca_cells, groups == "A"),
                         aes(x = PC1, y = PC2, colour = groups)) +
      geom_point(alpha = 0.5, size = 3)

    pca.cell.2 <- pca.cell.1 +
      geom_point(data = subset(pca_cells, groups != "A" ),
                 aes(x = PC1, y = PC2, colour = groups),
                 alpha = 0.8, size = 3) +
      scale_color_manual("groups", values = c("A" = "#8491B4B2", "B"="#E64B35FF")) +
      plotTheme("pca")

    # genes plot
    pl <- ggplot(D, aes(x = n, y = means)) + geom_point() +
      geom_text_repel(data = subset(D, n > (max(D$n) - 15)),
                      aes(n, means, label = rownames(D[D$n > (max(D$n)- 15),])),
                      nudge_y = 0.05, nudge_x = 0.05, direction = "x",
                      angle = 90, vjust = 0, segment.size = 0.2) +
      ggtitle(label = "B cell group genes mean expression") +
      plotTheme("genes")

    # UDE/nu plot
    nu_est = round(getNu(objCOTAN), digits = 7)

    plot.nu <- ggplot(pca_cells, aes(x = PC1,y = PC2, colour = log(nu_est))) +
      geom_point(size = 1, alpha = 0.8) +
      scale_color_gradient2(low = "#E64B35B2", mid = "#4DBBD5B2", high = "#3C5488B2",
                            midpoint = log(mean(nu_est)), name = "ln(nu)") +
      ggtitle("Cells PCA coloured by cells efficiency") +
      plotTheme("UDE")

    output <- list("object" = as(objCOTAN, "scCOTAN"),
                   "cl1" = cl1, "cl2" = cl2, "D" = D, "pca_cells" = pca_cells,
                   "pca.cell.2" = pca.cell.2, "genes.plot" = pl, "UDE.plot" = plot.nu)
    return(output)
  }
)
