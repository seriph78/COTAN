#' @details `reorderClusterization()` takes in a *clusterizations* and reorder
#'   its labels so that in the new order near labels indicate near clusters
#'   according to a *DEA* (or *Zero-One*) based distance
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName`
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#' @param reverse a flag to the output order
#' @param keepMinusOne a flag to decide whether to keep the cluster `"-1"`
#'   (representing the non-clustered cells) untouched
#' @param useDEA Boolean indicating whether to use the *DEA* to define the
#'   distance; alternatively it will use the average *Zero-One* counts, that is
#'   faster but less precise.
#' @param distance type of distance to use. Default is `"cosine"` for *DEA* and
#'   `"euclidean"` for *Zero-One*. Can be chosen among those supported by
#'   [parallelDist::parDist()]
#' @param hclustMethod It defaults is `"ward.D2"` but can be any of the methods
#'   defined by the [stats::hclust()] function.
#'
#' @returns `reorderClusterization()` returns a `list` with 2 elements:
#'   * `"clusters"` the newly reordered cluster labels array
#'   * `"coex"` the associated `COEX` `data.frame`
#'   * `"permMap"` the reordering mapping
#'
#' @export
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_null
#'
#' @importFrom stats hclust
#' @importFrom stats as.dist
#'
#' @rdname HandlingClusterizations
#'

reorderClusterization <- function(objCOTAN,
                                  clName = "", clusters = NULL, coexDF = NULL,
                                  reverse = FALSE, keepMinusOne = TRUE,
                                  useDEA = TRUE, distance = NULL,
                                  hclustMethod = "ward.D2") {
  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  clDist <- distancesBetweenClusters(objCOTAN, clName = clName,
                                     clusters = clusters, coexDF = coexDF,
                                     useDEA = useDEA, distance = distance)

  dummyList <- list("clusters" = factor(clusters), "coex" = coexDF,
                    "permMap" = set_names(labels(clDist), labels(clDist)))

  minuOnePos <- 0L
  if (keepMinusOne && any(clusters == "-1")) {
    minuOnePos <- which(labels(clDist) == "-1")
    # drop cluster '-1' from the distances
    clDist <- as.dist(as.matrix(clDist)[-minuOnePos, -minuOnePos, drop = FALSE])
  }

  if (length(labels(clDist)) <= 1L) {
    # too few clusters, nothing to reorder: return input as is
    return(dummyList)
  } else {
    rm(dummyList)
  }

  hc <- hclust(clDist, method = hclustMethod)

  # we exploit the rank(x) == order(order(x))
  perm <- order(hc[["order"]])

  if (isTRUE(reverse)) {
    perm <- (length(perm) + 1L) - perm
  }

  clNames <- hc[["labels"]]
  clMap <- set_names(clNames[perm], clNames)

  # handle cluster "-1" separately
  if (minuOnePos != 0L) {
    clMap[["-1"]] <- "-1"
  }

  logThis("Applied reordering to clusterization is:", logLevel = 1L)
  logThis(paste(paste0(names(clMap)), " -> ", paste0(clMap), collapse = ", "),
          logLevel = 1L)

  outputClusters <- factorToVector(factor(clusters))
  outputClusters <- set_names(clMap[outputClusters], names(outputClusters))

  if (is_empty(coexDF) && clName %in% getClusterizations(objCOTAN)) {
    coexDF <- getClusterizationData(objCOTAN, clName = clName)[["coex"]]
  }

  outputCoexDF <- coexDF
  if (!is_empty(coexDF)) {

    minusOnePosInCoex <- 0L
    if (minuOnePos != 0L) {
      minusOnePosInCoex <- which(colnames(coexDF) == "-1")
      outputCoexDF <- coexDF[, -minusOnePosInCoex]
    }

    colnames(outputCoexDF) <- clMap[colnames(outputCoexDF)]

    # Reorder the columns to match wanted order
    outputCoexDF <- outputCoexDF[, hc[["order"]]]

    if (minuOnePos != 0L) {
      outputCoexDF <- setColumnInDF(df = outputCoexDF, colName = "-1",
                                    colToSet = coexDF[, minusOnePosInCoex])
    }
  }

  return(list("clusters" = factor(outputClusters),
              "coex" = outputCoexDF, "permMap" = clMap))
}
