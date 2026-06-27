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
#' @param useDEA `r lifecycle::badge("deprecated")` Legacy cluster-tree scalar.
#'   Boolean indicating whether to use *DEA* profiles to define cluster
#'   distances; when `FALSE`, average *Zero-One* counts are used instead, which
#'   is faster but less precise. Use
#'   `clusterTreeOptions = ClusterTreeOptions(useDEA = ...)` instead.
#' @param distance `r lifecycle::badge("deprecated")` Legacy cluster-tree
#'   scalar. Distance method passed to [parallelDist::parDist()]. The effective
#'   default is `"cosine"` for *DEA* distances and `"euclidean"` for *Zero-One*
#'   distances. Use
#'   `clusterTreeOptions = ClusterTreeOptions(distance = ...)` instead.
#' @param hclustMethod `r lifecycle::badge("deprecated")` Legacy cluster-tree
#'   scalar. Clustering method passed to [stats::hclust()], with default
#'   `"ward.D2"`. Use
#'   `clusterTreeOptions = ClusterTreeOptions(hclustMethod = ...)` instead.
#' @param clusterTreeOptions a `ClusterTreeOptions` object controlling how
#'   distances between clusters are computed and how the hierarchical tree is
#'   built. This is the preferred interface for new code.
#'
#' @returns `reorderClusterization()` returns a `list` with 3 elements:
#'   * `"clusters"` the newly reordered cluster labels array
#'   * `"coex"` the associated `COEX` `data.frame`
#'   * `"permMap"` the reordering mapping
#'
#' @export
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang set_names
#'
#' @importFrom stats hclust
#' @importFrom stats as.dist
#'
#' @rdname HandlingClusterizations
#'
reorderClusterization <- function(objCOTAN,
                                  clName = "",
                                  clusters = NULL,
                                  coexDF = NULL,
                                  reverse = FALSE,
                                  keepMinusOne = TRUE,
                                  useDEA = TRUE,
                                  distance = NULL,
                                  hclustMethod = "ward.D2",
                                  clusterTreeOptions = NULL) {
  callArgs <- names(as.list(match.call(expand.dots = FALSE))[-1L])

  if (is.null(clusterTreeOptions)) {
    .warnDeprecatedPackArgs(
      functionName = .currentFunctionName(),
      callArgs = callArgs,
      packClass = "ClusterTreeOptions",
      replacementArg = "clusterTreeOptions",
      details = paste(
        "Use `clusterTreeOptions = ClusterTreeOptions(...)` to configure",
        "cluster-distance and hierarchical-tree parameters."
      )
    )
  }

  clusterTreeOptions <- resolveClusterTreeOptions(
    useDEA = useDEA,
    distance = distance,
    hclustMethod = hclustMethod,
    clusterTreeOptions = clusterTreeOptions
  )

  # picks up the last clusterization if none was given
  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName,
                           labels = clusters, isCond = FALSE)

  clDist <- distancesBetweenClusters(
    objCOTAN,
    clName = clName,
    clusters = clusters,
    coexDF = coexDF,
    clusterDistanceOptions = clusterTreeOptions
  )

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

  hc <- hclust(clDist, method = clusterTreeOptions@hclustMethod)

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
  logThis(paste(paste0(names(clMap)), " -> ", paste0(clMap),
                collapse = ", "), logLevel = 1L)

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
