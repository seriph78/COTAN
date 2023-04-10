# Conditions is a named list containing the condition name and in each arrays,
# the cell regarding that condition
# default null: just one condition

#' clustersSummaryPlot
#'
#' @param objCOTAN a `COTAN` object
#' @param condition name of the column in the `metaCells` dataframe containing
#'   the condition. Defaults to `NULL`
#' @param clName The name of the clusterization. If not given the last available
#'   clusterization will be used, as it is probably the most significant!
#' @param plotTitle The title to use for the returned plot
#'
#' @returns a `list` of a `data.frame` and a `ggplot` objects
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
#' data("test.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = test.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = FALSE,
#'                                          outDir = tempdir())
#'
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "clusters",
#'                               clusters = clusters)
#'
#' dataAndPlot <- clustersSummaryPlot(objCOTAN)
#'
#' @rdname clustersSummaryPlot
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






#' ClustersTreePlot
#'
#' @param objCOTAN
#' @param k numeric scalar (OR a vector) with the number of clusters the tree
#'   should be cut into.
#'
#' @returns the dendrogram
#'
#' @export
#' @import RColorBrewer
#'
#' @examples
ClustersTreePlot <- function(objCOTAN, k){
  cluster_data <- getClusterizationData(objCOTAN)[["coex"]]

  ######## This is the best: cosine dissimilarity
  Matrix <- as.matrix(t(cluster_data))
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  tree <- hclust(D_sim,method = "ward.D2")

  ############

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  ########
  dend <- as.dendrogram(tree)
  #colnames(df) <- str_remove_all(colnames(df), pattern = "cl.")
  cut = cutree(tree, k = k)
  dend =branches_color(dend,k=k,col=col_vector[1:k],groupLabels = T)
  dend =color_labels(dend,k=k)
  dend = set(dend = dend, "branches_lwd", 4)


return(dend)

}

