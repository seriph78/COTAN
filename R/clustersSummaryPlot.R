# Conditions is a named list containing the condition name and in each arrays,
# the cell regarding that condition
# default null: just one condition

#' @details `clustersSummaryPlot()` calculates various statistics about each
#'   cluster (with an optional further `condition` to separate the cells) and
#'   puts them together into a plot. The calculated statistics are:
#'   * "Cluster" the *cluster* **label**
#'   * "Condition" the further element to sub-divide the clusters
#'   * "CellNumber" the number of cells in the group
#'   * "MeanUDE" the average "UDE" in the group of cells
#'   * "MedianUDE" the median "UDE" in the group of cells
#'   * "ExpGenes25" the number of genes expressed in at the least 25% of the
#'     cells in the group
#'   * "ExpGenes" the number of genes expressed at the least once in any of the
#'     cells in the group
#'   * "CellPercentage" fraction of the cells with respect to the total cells
#'
#' @param objCOTAN a `COTAN` object
#' @param condition the name of a column in the `metaCells` `data.frame`
#'   containing the *condition*. This allows to further separate the cells in
#'   more sub-groups. When not given condition is assumed to be the same for all
#'   cells.
#'
#' @param clName The name of the clusterization. If not given the last available
#'   clusterization will be used, as it is probably the most significant!
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
#'
#' @rdname HandlingClusterizations
#'
clustersSummaryPlot <- function(objCOTAN, condition = NULL,
                                clName = "", plotTitle = "") {
  allClustNames <- getClusterizations(objCOTAN, keepPrefix = TRUE)

  emptyName <- !(length(clName) && any(nzchar(clName)))
  if (emptyName) {
    # pick last clusterization
    internalName <- allClustNames[length(allClustNames)]
  } else {
    internalName <- clName
    if (!startsWith(internalName, "CL_")) {
      internalName <- paste0("CL_", clName)
    }
    assert_that(internalName %in% allClustNames,
                msg = "Given cluster name is not among stored clusterizations")
  }

  metaCells <- getMetadataCells(objCOTAN)

  emptyCond <- !(length(condition) && any(nzchar(condition)))
  if (emptyCond || !(condition %in% colnames(metaCells))) {
    if (!emptyCond) {
      warning("Provided condition", condition,
              "is not available: will be ignored")
    }
    emptyCond <- TRUE
    condition <- NULL
  }

  df <- as.data.frame(table(metaCells[, c(internalName, condition),
                                      drop = FALSE]))

  if (emptyCond) {
    colnames(df) <- c("Cluster", "CellNumber")
    condition <- "Cond"
    df <- setColumnInDF(df, "NoCond", colName = condition)
    df <- df[, c("Cluster", condition, "CellNumber")]
  } else {
    colnames(df) <- c("Cluster", condition, "CellNumber")
  }

  assert_that(all(colnames(df) == c("Cluster", condition, "CellNumber")),
              msg = "Internal issue")

  df <- setColumnInDF(df, NA, colName = "MeanUDE")
  df <- setColumnInDF(df, NA, colName = "MedianUDE")
  df <- setColumnInDF(df, NA, colName = "ExpGenes25")
  df <- setColumnInDF(df, NA, colName = "ExpGenes")

  for (cond in unique(df[, condition])) {
    if (!emptyCond) {
      condPosInDF <- df[, condition] == cond
      condPosInMeta <- metaCells[[condition]] == cond
    } else {
      condPosInDF   <- rep_len(TRUE, length.out = nrow(df))
      condPosInMeta <- rep_len(TRUE, length.out = getNumCells(objCOTAN))
    }

    for (cl in levels(df[["Cluster"]])) {
      posInDF   <- (df[["Cluster"]] == cl) & condPosInDF
      posInMeta <- (metaCells[[internalName]] == cl) & condPosInMeta

      nu <- metaCells[["nu"]][posInMeta]
      df[posInDF, "MeanUDE"]   <- round(mean(nu),   digits = 2L)
      df[posInDF, "MedianUDE"] <- round(median(nu), digits = 2L)

      numTimesGeneIsExpr <- rowSums(getZeroOneProj(objCOTAN)[, posInMeta,
                                                             drop = FALSE])

      df[posInDF, "ExpGenes25"] <- sum(numTimesGeneIsExpr > length(nu) * 0.25)
      df[posInDF, "ExpGenes"]   <- sum(numTimesGeneIsExpr > 0L)
    }
  }

  cellsPerc <- round(df[["CellNumber"]] / getNumCells(objCOTAN) * 100.0,
                     digits = 2L)

  df <- setColumnInDF(df, colToSet = cellsPerc, colName = "CellPercentage")

  plotDF <- df %>%
    gather(keys, values, CellNumber, MeanUDE:CellPercentage)

  plotDF[["cond"]] <- plotDF[[condition]]

  plot <- ggplot(plotDF, aes(Cluster, values, fill = cond)) +
    geom_bar(position = "dodge", stat = "identity", color = "black") +
    geom_text(aes(x = Cluster, y = values, label = values,
                  vjust = 0.5, hjust = -0.1),
              position = position_dodge(width = 1.0)) +
    facet_wrap(~ keys, ncol = 6L, scales = "free") +
    scale_y_continuous(expand = expansion(mult = c(.05, .4))) +
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
#' @param kCuts the number of estimated *cluster* (this defines the high for the
#'   tree cut)
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
#' treePlot <- clustersTreePlot(objCOTAN, 2)
#'
#' @rdname HandlingClusterizations
#'
clustersTreePlot <- function(objCOTAN, kCuts,
                             clName = NULL,
                             distance = "cosine",
                             hclustMethod = "ward.D2") {

  emptyName <- !(length(clName) && any(nzchar(clName)))
  if (emptyName) {
    clNames <- getClusterizations(objCOTAN)
    clName <- clNames[length(clNames)]
    rm(clNames)
  }

  c(clusters, coexDF) %<-% getClusterizationData(objCOTAN, clName = clName)

  if (kCuts > length(unique(clusters))) {
    logThis("The number of cuts must be not more than the number of clusters",
            logLevel = 1L)
    kCuts <- length(unique(clusters))
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
