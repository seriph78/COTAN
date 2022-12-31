
#' geneSetEnrichment
#'
#' @description This function returns a cumulative score of enrichment in a
#'   cluster over a gene set
#'
#' @details Calculates \eqn{\frac{1}{n}\sum_i(1-e^{-\theta X_i})}, where the
#'   \eqn{X_i} are the positive values from [DEAOnClusters()] and \eqn{\theta =
#'   -\frac{1}{0.1} \ln(0.25)}
#'
#' @param objCOTAN a `COTAN` object
#' @param clustersCoex the `data.frame` for the increased or decreased
#'   expression of every gene in the cluster compared to the whole background
#' @param groupMarkers a named `list` of arrays of genes
#'
#' @returns a `data.frame` with the cumulative score
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom stringr str_remove
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = FALSE,
#'                                          outDir = tempdir())
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#' coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)[["coex"]]
#' objCOTAN <- addClusterization(objCOTAN, clName = "clusters",
#'                               clusters = clusters, coexDF = coexDF)
#' groupMarkers <- list(G1 = c("Pcbp2", "Snrpe", "Nfyb"), G2 = c("Prpf40a", "Ergic2"),
#'                      G3 = c("Ncl", "Cd47", "Macrod2", "Fth1", "Supt16"))
#' enrichment <- geneSetEnrichment(clustersCoex = coexDF,
#'                                 groupMarkers = groupMarkers)
#'
#' @rdname geneSetEnrichment
#'
geneSetEnrichment <- function(clustersCoex, groupMarkers) {
  # It might be possible that the column names have the old extra
  # prefix 'cl.'. It will be remove in such cases.
  clTagsMap <- set_names(colnames(clustersCoex),
                         str_remove(colnames(clustersCoex), "cl\\."))
  clustersTags <- names(clTagsMap)

  df <- as.data.frame(matrix(nrow = length(groupMarkers),
                             ncol = ncol(clustersCoex) + 2))

  colnames(df) <- c(clustersTags, "N. detected", "N. total")
  rownames(df) <- names(groupMarkers)

  # TODO: add comment on this constant in the @details above!
  teta <- -1/0.1 * log(0.25)

  # not_assigned_clusters <- NA
  for (groupName in names(groupMarkers)) {
    genes <- unlist(groupMarkers[[groupName]])

    ex <- clustersCoex[rownames(clustersCoex) %in% genes, , drop = FALSE]
    numDetected <- nrow(ex)

    logThis(paste0("In group ", groupName, " there are ", numDetected,
                   " detected over ", length(genes), " genes"), logLevel = 3)

    # drop reductions
    ex[ex < 0 & !is.na(ex)] <- 0
    ex <- 1 - exp(-teta * ex)

    for (cl in clustersTags) {
      df[groupName, cl] <- sum(ex[, clTagsMap[cl]], na.rm = TRUE) / length(genes)
    }

    df[groupName, "N. detected"] <- numDetected
    df[groupName, "N. total"] <- length(genes)
  }

  return(df)
}
