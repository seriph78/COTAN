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
#'
#' @export
#'
#' @import gsubfn
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
#' list[objCOTAN, data] <- clean(ERCC.cotan)
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

    pcaCells <- cbind(pcaCells, groups)

    return( list("objCOTAN" = objCOTAN,
                 "data"     = list("cluster1" = cluster1, "cluster2" = cluster2,
                                   "D" = D, "pcaCells" = pcaCells)) )
  }
)
