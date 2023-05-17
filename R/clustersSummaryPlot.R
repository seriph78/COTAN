
#' @details `clustersSummaryData()` calculates various statistics about each
#'   cluster (with an optional further `condition` to separate the cells).
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#'
#' @returns `clustersSummaryData()` returns a `data.frame`  with the following
#'   statistics: The calculated statistics are:
#'   * `clName` the *cluster* **labels**
#'   * `condName` the relevant condition (to further sub-divide the *clusters*)
#'   * "CellNumber" the number of cells in the group
#'   * "MeanUDE" the average "UDE" in the group of cells
#'   * "MedianUDE" the median "UDE" in the group of cells
#'   * "ExpGenes25" the number of genes expressed in at the least 25% of the
#'   cells in the group
#'   * "ExpGenes" the number of genes expressed at the least once in any of the
#'   cells in the group
#'   * "CellPercentage" fraction of the cells with respect to the total cells
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' summaryData <- clustersSummaryData(objCOTAN)
#'
#' @rdname HandlingClusterizations
#'
clustersSummaryData <- function(objCOTAN, clName = "", clusters = NULL,
                                condName = "", conditions = NULL) {

  c(clName, clusters) %<-%
    normalizeNameAndLabels(objCOTAN, name = clName, labels = clusters)

  c(condName, conditions) %<-%
    normalizeNameAndLabels(objCOTAN, name = condName,
                           labels = conditions, isCond = TRUE)

  if (isEmptyName(condName)) {
    condName <- "Cond"
  }

  df <- as.data.frame(cbind(factorToVector(clusters),
                            factorToVector(conditions)))
  df <- as.data.frame(table(df))
  assert_that(ncol(df) == 3L, msg = "Internal error creating frequency table")

  colnames(df) <- c(clName, condName, "CellNumber")

  cellsPerc <- round(df[["CellNumber"]] / getNumCells(objCOTAN) * 100.0,
                     digits = 1L)

  df <- setColumnInDF(df, colToSet = cellsPerc, colName = "CellPercentage")

  df <- setColumnInDF(df, NA, colName = "MeanUDE")
  df <- setColumnInDF(df, NA, colName = "MedianUDE")
  df <- setColumnInDF(df, NA, colName = "ExpGenes25")
  df <- setColumnInDF(df, NA, colName = "ExpGenes")

  for (cond in levels(factor(df[[condName]]))) {
    condPosInDF <- df[[condName]] == cond
    condPosInMeta <- conditions == cond

    for (cl in levels(factor(df[[clName]]))) {
      posInDF   <- (df[[clName]] == cl) & condPosInDF
      posInMeta <- (clusters == cl) & condPosInMeta

      nu <- getNu(objCOTAN)[posInMeta]
      df[posInDF, "MeanUDE"]   <- round(mean(nu),   digits = 2L)
      df[posInDF, "MedianUDE"] <- round(median(nu), digits = 2L)

      numTimesGeneIsExpr <- rowSums(getZeroOneProj(objCOTAN)[, posInMeta,
                                                             drop = FALSE])

      df[posInDF, "ExpGenes25"] <- sum(numTimesGeneIsExpr > length(nu) * 0.25)
      df[posInDF, "ExpGenes"]   <- sum(numTimesGeneIsExpr > 0L)
    }
  }

  return(df)
}

#' @details `clustersSummaryPlot()` calculates various statistics about each
#'   cluster via [clustersSummaryData()] and puts them together into a plot.
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#' @param condName The name of a condition in the `COTAN` object to further
#'   separate the cells in more sub-groups. When no condition is given it is
#'   assumed to be the same for all cells (no further sub-divisions)
#' @param conditions The *conditions* to use. If given it will take precedence
#'   on the one indicated by `condName` that will only indicate the relevant
#'   column name in the returned `data.frame`
#' @param plotTitle The title to use for the returned plot
#'
#' @returns `clustersSummaryPlot()` returns a `list` with a `data.frame` and a
#'   `ggplot` objects
#'   * "data" contains the data,
#'   * "plot" is the returned plot
#'
#' @importFrom tidyr gather
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 position_dodge
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_fill_brewer
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @examples
#' dataAndPlot <- clustersSummaryPlot(objCOTAN)
#' plot(dataAndPlot[["plot"]])
#'
#' @rdname HandlingClusterizations
#'
clustersSummaryPlot <- function(objCOTAN, clName = "", clusters = NULL,
                                condName = "", conditions = NULL,
                                plotTitle = "") {

  df <- clustersSummaryData(objCOTAN, clName = clName, clusters = clusters,
                            condName = condName, conditions = conditions)

  plotDF <- df %>%
    gather(keys, values, CellNumber:ExpGenes)

  # normalize col names
  cNames <- colnames(df)
  plotDF[["Cluster"]] <- plotDF[[cNames[[1L]]]]
  plotDF[["Condition"]] <- plotDF[[cNames[[2L]]]]

  plot <- ggplot(plotDF, aes(Cluster, values, fill = Condition)) +
    geom_bar(position = "dodge", stat = "identity", color = "black") +
    geom_text(aes(x = Cluster, y = values, label = values,
                  vjust = 0.5, hjust = -0.1),
              position = position_dodge(width = 1.0)) +
    facet_wrap(~ keys, ncol = 6L, scales = "free") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.4))) +
    coord_flip() +
    theme_classic() +
    ggtitle(plotTitle) +
    scale_fill_brewer(palette = "Accent")

  return(list("data" = df, "plot" = plot))
}



#' @details `clustersTreePlot()` returns the `dendogram` plot where the given
#'   *clusters* are placed on the base of their relative distance. Also if
#'   needed calculates and stores the `DEA` of the relevant *clusterization*.
#'
#' @param objCOTAN a `COTAN` object
#' @param kCuts the number of estimated *cluster* (this defines the height for
#'   the tree cut)
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be returned, as it is probably the most
#'   significant!
#' @param distance type of distance to use (default is `cosine`, `euclidean` is
#'   also available)
#' @param hclustMethod default is "ward.D2" but can be any method defined by
#'   [stats::hclust()] function
#'
#' @returns `clustersTreePlot()` returns a list with 2 objects:
#'  * "dend" a `ggplot2` object representing the `dendrogram` plot
#'  * "objCOTAN" the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @importFrom dendextend set
#' @importFrom dendextend color_labels
#' @importFrom dendextend branches_color
#'
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats as.dendrogram
#'
#' @export
#'
#' @examples
#' treePlotAndObj <- clustersTreePlot(objCOTAN, 2)
#' objCOTAN <- treePlotAndObj[["objCOTAN"]]
#' plot(treePlotAndObj[["dend"]])
#'
#' @rdname HandlingClusterizations
#'
clustersTreePlot <- function(objCOTAN, kCuts,
                             clName = "",
                             distance = "cosine",
                             hclustMethod = "ward.D2") {
  # pick last if no name was given
  clName <- getClusterizationName(objCOTAN, clName = clName)
  c(clusters, coexDF) %<-% getClusterizationData(objCOTAN, clName = clName)
  assert_that(inherits(clusters, "factor"),
              msg = "Internal error - clusters must be factors")

  if (kCuts > nlevels(clusters)) {
    logThis("The number of cuts must be not more than the number of clusters",
            logLevel = 1L)
    kCuts <- nlevels(clusters)
  }

  colVector <- getColorsVector(kCuts)

  if (is_empty(coexDF)) {
    logThis("Coex dataframe is missing: will be calculated and stored",
            logLevel = 1L)
    coexDF <- DEAOnClusters(objCOTAN, clusters = clusters)[["coex"]]
    objCOTAN <- addClusterizationCoex(objCOTAN, clName = clName,
                                      coexDF = coexDF)
  }
  rm(clusters)

  #merge small cluster based on distances
  if (distance == "cosine") {
    # This is the best: cosine dissimilarity
    coexDist <- cosineDissimilarity(as.matrix(coexDF))
  } else if (distance == "euclidean") {
    coexDist <- dist(t(as.matrix(coexDF)))
  } else {
    stop("only 'cosine' and 'euclidean' distances are supported")
  }

  hcNorm <- hclust(coexDist, method = hclustMethod)

  dend <- as.dendrogram(hcNorm)
  dend <- branches_color(dend, k = kCuts, col = colVector, groupLabels = TRUE)
  dend <- color_labels(dend, k = kCuts)
  dend <- set(dend = dend, "branches_lwd", 4L)

  return(list("dend" = dend, "objCOTAN" = objCOTAN))
}
