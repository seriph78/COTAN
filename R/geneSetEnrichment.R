
#' @details `geneSetEnrichment()` returns a cumulative score of enrichment in a
#'   *cluster* over a gene set. In formulae it calculates
#'   \eqn{\frac{1}{n}\sum_i(1-e^{-\theta X_i})}, where the \eqn{X_i} are the
#'   positive values from [DEAOnClusters()] and \eqn{\theta = -\frac{1}{0.1}
#'   \ln(0.25)}
#'
#' @param clustersCoex the `COEX` `data.frame`
#' @param groupMarkers an optional named `list` with an element for each group
#'   comprised of one or more marker genes
#'
#' @returns `geneSetEnrichment()` returns a `data.frame` with the cumulative
#'   score
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom stringr str_remove
#' @importFrom stringr fixed
#'
#' @rdname HandlingClusterizations
#'
geneSetEnrichment <- function(clustersCoex, groupMarkers = list()) {
  assert_that(isa(groupMarkers, "list"))

  # It might be possible that the column names have the old extra
  # prefix 'cl.'. It will be remove in such cases.
  clTagsMap <- set_names(colnames(clustersCoex),
                         str_remove(colnames(clustersCoex), fixed("cl.")))
  clustersTags <- names(clTagsMap)

  gDf <- as.data.frame(matrix(nrow = length(groupMarkers),
                              ncol = ncol(clustersCoex) + 2L))

  colnames(gDf) <- c(clustersTags, "N. detected", "N. total")
  rownames(gDf) <- names(groupMarkers)

  # TODO: add comment on this constant in the @details above!
  theta <- -1.0 / 0.1 * log(0.25) # ~ 13.86294

  # not_assigned_clusters <- NA
  for (groupName in names(groupMarkers)) {
    genes <- unlist(groupMarkers[[groupName]])

    ex <- clustersCoex[rownames(clustersCoex) %in% genes, , drop = FALSE]
    numDetected <- nrow(ex)

    logThis(paste0("In group ", groupName, " there are ", numDetected,
                   " detected over ", length(genes), " genes"), logLevel = 2L)

    # drop reductions
    ex[ex < 0.0 & !is.na(ex)] <- 0.0
    ex <- 1.0 - exp(-theta * ex)

    for (cl in clustersTags) {
      gDf[groupName, cl] <-
        sum(ex[, clTagsMap[cl], drop = TRUE], na.rm = TRUE) / length(genes)
    }

    gDf[groupName, "N. detected"] <- numDetected
    gDf[groupName, "N. total"] <- length(genes)
  }

  return(gDf)
}
