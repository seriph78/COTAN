#' clean
#'
#' Main function that can be used to check and clean the dataset.
#' It also produces and stores the estimators for nu and lambda
#'
#' @param objCOTAN COTAN object
#' @param calcExtraData Boolean flag. When FALSE all PCA/cluster related
#' evaluation will be skipped and the returned "data" list will be empty
#'
#' @return lists of objects...
#' "object" is the COTAN object,
#' "data" a empty list when cleanExtraData is FALSE otherwise a list with:
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
  function(objCOTAN, calcExtraData = TRUE) {

    # We want to discard genes having less than 3 non-zero counts per 1000 cells
    threshold <- round(getNumCells(objCOTAN) * 3 / 1000, digits = 0)
    genesToDrop <- getGenes(objCOTAN)[rowSums(getZeroOneProj(objCOTAN)) <= threshold]

    objCOTAN <- dropGenesCells(objCOTAN, genes = genesToDrop)

    # We want to discard cells having less than 2 non-zero counts per 1000 genes
    threshold <- round(getNumGenes(objCOTAN) * 2 / 1000, digits = 0)
    cellsToDrop <- getCells(objCOTAN)[colSums(getZeroOneProj(objCOTAN)) <= threshold]

    objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)

    print(paste0("Genes/cells selection done:",
                 " dropped [", length(genesToDrop), "] genes",
                 " and [", length(cellsToDrop), "] cells"))

    print(paste0("Working on [", getNumGenes(objCOTAN), "]",
                 " genes and [", getNumCells(objCOTAN), "] cells"))

    objCOTAN <- estimateLambdaLinear(objCOTAN)
    objCOTAN <- estimateNuLinear(objCOTAN)
    objCOTAN <- findHousekeepingGenes(objCOTAN)

    # genesMeans <- getNu(objCOTAN)
    # genesRng <- round(getNumGenes(objCOTAN)     / 2, digits = 0)
    #           : round(getNumGenes(objCOTAN) * 3 / 4, digits = 0)
    # genesMax <- names(sort(genesMeans, decreasing = TRUE)[genesRng])

    gc()

    if (isTRUE(calcExtraData)) {
      print("PCA: START")

      rawNorm <- getNormalizedData(objCOTAN)

      pcaCells <- irlba::prcomp_irlba(t(rawNorm), n = 5)[["x"]]
      rownames(pcaCells) <- getCells(objCOTAN)

      distCells <- stats::dist(scale(pcaCells), method = "euclidean") # mhalanobis

      pcaCells <- as.data.frame(pcaCells)

      print("PCA: DONE")

      print("Hierarchical clustering: START")

      hcCells <- hclust(distCells, method = "complete")
      rm(distCells)

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
      gc()

      toClust <- as.matrix(round(rawNorm, digits = 4))

      if (!identical(colnames(toClust), names(groups))) {
        stop("Error in the cell names")
      }

      # ---- next: to check which genes are specific for the B group of cells
      {
        B <- as.data.frame(toClust[, pos2])
        colnames(B) <- colnames(toClust)[pos2]
        #print(utils::head(B, 15))

        rm(toClust)

        D <- data.frame("means" = rowMeans(B), "n" = NA)
        rownames(D) <- rownames(B)
        D <- rownames_to_column(D)
        rm(B)

        D <- D[order(D[["means"]], decreasing = TRUE), ]
        rownames(D) <- c()
        D <- column_to_rownames(D)

        D <- D[D[["means"]] > 0,]
        D[["n"]] <- seq_along(D[["means"]])
      }
      gc()

      pcaCells <- cbind(pcaCells, groups)

      data <- list("cluster1" = cluster1, "cluster2" = cluster2,
                   "D" = D, "pcaCells" = pcaCells)
    }
    else {
      data <- list()
    }

    return( list("objCOTAN" = objCOTAN, "data" = data) )
  }
)
