#' @details `reorderClusterization()` takes in a *clusterizations* and reorder
#'   its labels so that in the new order near labels indicate near clusters
#'   according to a `DEA` based distance
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters The *clusterization* to merge. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param coexDF a `data.frame` where each column indicates the `COEX` for each
#'   of the *clusters* of the *clusterization*
#' @param reverse a flag to the output order
#' @param keepMinuOne a flag to decide whether to keep the cluster `"-1"`
#'   (representing the non-clustered cells) untouched
#' @param distance type of distance to use (default is `"cosine"`, `"euclidean"`
#'   and the others from [parallelDist::parDist()] are also available)
#' @param hclustMethod It defaults is `"ward.D2"` but can be any of the methods
#'   defined by the [stats::hclust()] function.
#'
#' @returns `reorderClusterization()` returns a `list` with 2 elements:
#'   * "clusters" the newly reordered cluster labels array
#'   * "coex" the associated `COEX` `data.frame`
#'
#' @export
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_null
#'
#' @importFrom parallelDist parDist
#'
#' @importFrom stats hclust
#'
#' @rdname UniformClusters
#'

reorderClusterization <- function(objCOTAN,
                                  clusters = NULL,
                                  coexDF = NULL,
                                  reverse = FALSE,
                                  keepMinusOne = TRUE,
                                  distance = "cosine",
                                  hclustMethod = "ward.D2") {
  if (is_empty(clusters)) {
    # pick the last clusterization
    c(clusters, coexDF) %<-% getClusterizationData(objCOTAN)
  }

  if (is_empty(coexDF)) {
    coexDF <- DEAOnClusters(objCOTAN, clusters)
  }

  # exclude cluster "-1"
  minusOneClCoex <- NULL
  if (keepMinusOne && any(clusters == "-1")) {
    col <- which(colnames(coexDF) == "-1")
    minusOneClCoex <- coexDF[["-1"]]
    coexDF <- coexDF[, -col]
  }

  # DEA based distance
  coexDist <- parDist(t(as.matrix(coexDF)), method = distance)

  hc <- hclust(coexDist, method = hclustMethod)

  # we exploit the rank(x) == order(order(x))
  perm <- order(hc[["order"]])

  if (isTRUE(reverse)) {
    perm <- (length(perm) + 1L) - perm
  }

  clNames <- hc[["labels"]]
  clMap <- set_names(clNames[perm], clNames)

  if (keepMinusOne && !is_empty(minusOneClCoex)) {
    clMap[["-1"]] <- "-1"
  }

  logThis("Applied reordering to clusterization is:", logLevel = 1L)
  logThis(paste(paste0(names(clMap)), " -> ", paste0(clMap), collapse = ", "),
          logLevel = 1L)

  outputClusters <- factor(clusters)
  outputCoexDF <- coexDF

  levels(outputClusters) <- clMap[levels(outputClusters)]
  colnames(outputCoexDF) <- clMap[colnames(coexDF)]
  # Reorder the columns to match wanted hc[["order"]]
  outputCoexDF <- outputCoexDF[, hc[["order"]]]

  # restore cluster "-1"
  if (keepMinusOne && !is_empty(minusOneClCoex)) {
    outputCoexDF[["-1"]] <- minusOneClCoex
  }

  return(list("clusters" = outputClusters, "coex" = outputCoexDF))
}
