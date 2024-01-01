#' pcaGenesCoex
#'
#' @details `pcaGenesCoex()` calculates the **PCA** of `genesCoex`, selecting
#'   the genes with highest **GDI** as features
#'
#' @param objCOTAN a `COTAN` object
#' @param pcaDim the dimension of the **PCA** vectors
#' @param varThreshold genes with **GDI** less than `varThreshold` are excluded
#'   from driving features (`genesCoex` columns)
#' @param idThreshold genes with **GDI** less than `idThreshold` are excluded
#'   from the genes' list (`genesCoex` rows)
#' @param actCoexSpace Boolean, when `TRUE` the transformation to coex space is
#'   applied
#'
#' @returns `pcaGenesCoex()` returns a `list` with two elements:
#'   * `"pca"` the PCA
#'   * `"varExplained"` the variance explained for each component
#'   * `"rotation"` the PCA rotation matrix
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom PCAtools pca
#'
#' @rdname MetaGenes
#'
pcaGenesCoex <- function(objCOTAN, pcaDim,
                         varThreshold = 1.5,
                         idThreshold = 1.25,
                         actCoexSpace = FALSE) {
  GDI <- set_names(calculateGDI(objCOTAN)[["GDI"]], getGenes(objCOTAN))
  GDI <- sort(GDI, decreasing = TRUE)

  varGenes <- names(GDI[GDI >= varThreshold])
  idGenes  <- names(GDI[GDI >= idThreshold])

  genesCoex <- getGenesCoex(objCOTAN, zeroDiagonal = TRUE)[idGenes, varGenes]

  logThis(paste0("Genes coex dimensions: (",
                 paste(dim(genesCoex), collapse = ", "), ")"), logLevel = 2)

  logThis(paste0("GDI mean of feature genes:", mean(GDI[varGenes])),
          logLevel = 3)
  logThis(paste0("GDI mean of selected genes:", mean(GDI[idGenes])),
          logLevel = 3)

  logThis("Use coexpresion space: ", logLevel = 2, appendLF = FALSE)
  if (isTRUE(actCoexSpace)) {
    genesCoex <- tanh(genesCoex * sqrt(getNumCells(objCOTAN)))
    logThis("yes", logLevel = 2)
  } else {
    logThis("no", logLevel = 2)
  }

  logThis("PCA starts", logLevel = 2)
  pcaStruct <- pca(mat = genesCoex, rank = pcaDim, scale  = FALSE,
                   transposed = FALSE, BSPARAM = IrlbaParam())

  rownames(pcaStruct[["rotated"]]) <- rownames(genesCoex)

  cat("Dimensions of transformed matrix:", dim(pcaStruct[["rotated"]]), "\n")
  cat("Dimensions of rotation matrix:", dim(pcaStruct[["loadings"]]), "\n")

  return(list("pca" = pcaStruct[["rotated"]],
              "varExplained" = pcaStruct[["variance"]],
              "rotation" = pcaStruct[["loadings"]]))
}
