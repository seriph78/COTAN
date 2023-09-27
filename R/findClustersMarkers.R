
#' @details `findClustersMarkers()` takes in a `COTAN` object and a
#'   *clusterization* and produces a `data.frame` with the `n` most positively
#'   enriched and the `n` most negatively enriched genes for each *cluster*. The
#'   function also provides whether and the found genes are in the given
#'   `markers` list or not. It also returns the *p-value* and the *adjusted*
#'   *p-value* using the [stats::p.adjust()]
#'
#' @param objCOTAN a `COTAN` object
#' @param n the number of extreme `COEX` values to return
#' @param markers a `list` of marker genes
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName`
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#' @param pValueDF a `data.frame` with *In/Out* *p-value* based on the `COEX`.
#'   E.G. the result of a call to [pValueFromDEA()]
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
    objCOTAN, n = 10L, markers = NULL,
    clName = "", clusters = NULL, coexDF = NULL, pValueDF = NULL,
    deltaExp = NULL, method = "bonferroni") {
  logThis("findClustersMarkers - START", logLevel = 2L)

  marks <- unlist(markers)

  assert_that(is_empty(marks) || any(marks %in% getGenes(objCOTAN)),
              msg = "None of the given markers is present in the data")

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  if (is_empty(coexDF)) {
    if (clName %in% getClusterizations(objCOTAN)) {
      coexDF <- getClusterizationData(objCOTAN, clName = clName)[["coex"]]
    }
    if (is_empty(coexDF)) {
      coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)
    }
  }

  if (is_empty(pValueDF)) {
    pValueDF <- pValueFromDEA(coexDF, getNumCells(objCOTAN))
  }

  if (is_empty(deltaExp)) {
    deltaExp <- clustersDeltaExpression(objCOTAN, clName = clName,
                                        clusters = clusters)
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
