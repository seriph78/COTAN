
#' @details `findClustersMarkers()` takes in a `COTAN` object and a
#'   *clusterization* and produces a `data.frame` with the `n` most positively
#'   enriched and the `n` most negatively enriched genes for each *cluster*. The
#'   function also provides whether and the found genes are in the given
#'   `markers` list or not. It also returns the *adjusted p-value* for
#'   multi-tests using the [stats::p.adjust()]
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
#' @param method *p-value* multi-test adjustment method. Defaults to
#'   `"bonferroni"`; use `"none"` for no adjustment
#'
#' @returns `findClustersMarkers()` returns a `data.frame` containing `n` genes
#'   for each *cluster* scoring top/bottom `COEX` scores. The `data.frame` also
#'   contains:
#'   * `"CL"` the cluster
#'   * `"Gene"` the gene
#'   * `"Score"` the `COEX` score of the gene
#'   * `"adjPVal"` the *p-values* associated to the `COEX`
#'     adjusted for *multi-testing*
#'   * `"DEA"` the differential expression of the gene
#'   * `"IsMarker"` whether the gene is among the given markers
#'   * `"logFoldCh"` the *log-fold-change* of the gene expression inside versus
#'     outside the cluster from [logFoldChangeOnClusters()]
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
    clName = "", clusters = NULL,
    coexDF = NULL, method = "bonferroni") {
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

  adjPValueDF <- pValueFromDEA(coexDF, numCells = getNumCells(objCOTAN),
                               method = method)

  deltaExp <- clustersDeltaExpression(objCOTAN, clName = clName,
                                      clusters = clusters)

  lfcDF <- logFoldChangeOnClusters(objCOTAN, clusters = clusters)

  assert_that(all(rownames(deltaExp)    == rownames(coexDF) &
                  rownames(adjPValueDF) == rownames(coexDF) &
                  getGenes(objCOTAN)    == rownames(coexDF)),
              msg = paste("Inconsistent data-frames passed in",
                          "for 'coex', 'p-value' or 'delta expression'"))

  retDF <- as.data.frame(matrix(data = NA, nrow = 0L, ncol = 7L))
  colnames(retDF) <- c("CL", "Gene", "Score", "adjPVal",
                       "DEA", "IsMarker", "logFoldCh")

  for (cl in unique(colnames(coexDF))) {
    for (type in c("min", "max")) {
      tmpDF <- as.data.frame(matrix(data = NA, nrow = n, ncol = ncol(retDF)))
      colnames(tmpDF) <- colnames(retDF)

      # Get the first n minimum/maximum scores for each cluster
      sortedPos <- order(coexDF[, cl], decreasing = (type == "max"))[1L:n]

      tmpDF[["CL"]]        <- cl
      tmpDF[["Gene"]]      <- rownames(coexDF)[sortedPos]
      tmpDF[["Score"]]     <- coexDF[sortedPos, cl]
      tmpDF[["adjPVal"]]   <- adjPValueDF[sortedPos, cl]
      tmpDF[["DEA"]]       <- deltaExp[sortedPos, cl]
      tmpDF[["IsMarker"]]  <- as.integer(tmpDF[["Gene"]] %in% marks)
      tmpDF[["logFoldCh"]] <- lfcDF[sortedPos, cl]

      retDF <- rbind(retDF, tmpDF)

      rm(tmpDF)
    }
  }

  logThis("findClustersMarkers - DONE", logLevel = 2L)

  return(retDF)
}
