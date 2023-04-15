
#' @details `geneSetEnrichment()` returns a cumulative score of enrichment in a
#'   *cluster* over a gene set. In formulae it calculates
#'   \eqn{\frac{1}{n}\sum_i(1-e^{-\theta X_i})}, where the \eqn{X_i} are the
#'   positive values from [DEAOnClusters()] and \eqn{\theta = -\frac{1}{0.1}
#'   \ln(0.25)}
#'
#' @param clustersCoex the `COEX` `data.frame`
#' @param groupMarkers a named `list` of arrays of genes
#'
#' @returns `geneSetEnrichment()` returns a `data.frame` with the cumulative
#'   score
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom stringr str_remove
#' @importFrom stringr fixed
#'
#' @rdname HandlingClusterizations
#'
geneSetEnrichment <- function(clustersCoex, groupMarkers) {
  # It might be possible that the column names have the old extra
  # prefix 'cl.'. It will be remove in such cases.
  clTagsMap <- set_names(colnames(clustersCoex),
                         str_remove(colnames(clustersCoex), fixed("cl.")))
  clustersTags <- names(clTagsMap)

  df <- as.data.frame(matrix(nrow = length(groupMarkers),
                             ncol = ncol(clustersCoex) + 2L))

  colnames(df) <- c(clustersTags, "N. detected", "N. total")
  rownames(df) <- names(groupMarkers)

  # TODO: add comment on this constant in the @details above!
  teta <- -1.0 / 0.1 * log(0.25)

  # not_assigned_clusters <- NA
  for (groupName in names(groupMarkers)) {
    genes <- unlist(groupMarkers[[groupName]])

    ex <- clustersCoex[rownames(clustersCoex) %in% genes, , drop = FALSE]
    numDetected <- nrow(ex)

    logThis(paste0("In group ", groupName, " there are ", numDetected,
                   " detected over ", length(genes), " genes"), logLevel = 3L)

    # drop reductions
    ex[ex < 0.0 & !is.na(ex)] <- 0.0
    ex <- 1.0 - exp(-teta * ex)

    for (cl in clustersTags) {
      df[groupName, cl] <-
        sum(ex[, clTagsMap[cl]], na.rm = TRUE) / length(genes)
    }

    df[groupName, "N. detected"] <- numDetected
    df[groupName, "N. total"] <- length(genes)
  }

  return(df)
}
