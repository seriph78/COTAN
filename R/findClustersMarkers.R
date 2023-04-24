
#' @details `findClustersMarkers()` takes in a `COTAN` object and a
#'   *clusterization* and produces a `data.frame` with the `n` most positively
#'   enriched and the `n` most negatively enriched genes for each *cluster*. The
#'   function also provides whether and the found genes are in the given
#'   `markers` list or not. It also returns the *p-value* and the *adjusted*
#'   *p-value* using the [stats::p.adjust()]
#'
#' @param objCOTAN a `COTAN` object
#' @param n the number of extreme `COEX` values to return
#' @param clusters a *clusterization*
#' @param markers a `list` of marker genes
#' @param coexDF a `data.frame` with *In/Out* `COEX`. E.G. the result of a call
#'   to [DEAOnClusters()]
#' @param pValueDF a `data.frame` with *In/Out* *p-value* based on the `COEX`.
#'   E.G. the result of a call to `DEAOnClusters()`
#' @param deltaExp a `data.frame` with the *delta-expression* in a *cluster*.
#'   E.G. the result of a call to [clustersDeltaExpression()]
#' @param method *p-value* adjustment method. Defaults to `"bonferroni"`
#'
#' @returns `findClustersMarkers()` returns a `data.frame` containing `n`
#'   top/bottom `COEX` scores for each *cluster*
#'
#' @importFrom stats p.adjust
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' clMarkers <- findClustersMarkers(objCOTAN, clusters = clusters)
#'
#' @rdname HandlingClusterizations
#'
findClustersMarkers <- function(
    objCOTAN, n = 10L, clusters = NULL, markers = NULL,
    coexDF = NULL, pValueDF = NULL, deltaExp = NULL,
    method = "bonferroni") {
  logThis("findClustersMarkers - START", logLevel = 2L)

  marks <- unlist(markers)

  assert_that(is_empty(marks) || any(marks %in% getGenes(objCOTAN)),
              msg = "None of the given markers is present in the data")

  if (is.null(coexDF) || is.null(pValueDF)) {
    # picks up the last clusterization if none was given
    DEA <- DEAOnClusters(objCOTAN, clusters = clusters)
    coexDF <- DEA[["coex"]]
    pValueDF <- DEA[["p-value"]]
  }

  if (is.null(deltaExp)) {
    # picks up the last clusterization if none was given
    deltaExp <- clustersDeltaExpression(objCOTAN, clusters = clusters)
  }

  assert_that(all(rownames(deltaExp) == rownames(coexDF) &
                  rownames(pValueDF) == rownames(coexDF) &
                  getGenes(objCOTAN) == rownames(coexDF)),
              msg = paste("Inconsistent data-frames passed in",
                          "for 'coex', 'p-value' or 'delta expression'"))

  retDF <- as.data.frame(matrix(data = NA, nrow = 0L, ncol = 7L))
  colnames(retDF) <- c("CL", "Gene", "Score", "pVal",
                       "adjPVal", "DEA", "IsMarker")

  for (cl in unique(colnames(coexDF))) {
    for (type in c("min", "max")) {
      tmpDF <- as.data.frame(matrix(data = NA, nrow = n, ncol = 7L))
      colnames(tmpDF) <- colnames(retDF)

      # Get the first n minimum/maximum scores for each cluster
      sortedPos <- order(coexDF[, cl], decreasing = (type == "max"))[1L:n]
      extrDF <- coexDF[sortedPos, ]

      tmpDF[["CL"]]       <- cl
      tmpDF[["Gene"]]     <- rownames(extrDF)
      tmpDF[["Score"]]    <- extrDF[, cl]
      tmpDF[["pVal"]]     <- pValueDF[sortedPos, cl]
      tmpDF[["adjPVal"]]  <- p.adjust(tmpDF[["pVal"]], method = method,
                                      n = getNumGenes(objCOTAN))
      tmpDF[["DEA"]]      <- deltaExp[sortedPos, cl]
      tmpDF[["IsMarker"]] <- as.integer(rownames(extrDF) %in% marks)

      retDF <- rbind(retDF, tmpDF)

      rm(extrDF, tmpDF)
    }
  }

  logThis("findClustersMarkers - DONE", logLevel = 2L)

  return(retDF)
}
