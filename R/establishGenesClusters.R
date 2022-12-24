
#' genesCoexSpace
#'
#' @description To make the `GDI` more specific, it may be desirable to restrict
#'   the set of genes against which `GDI` is computed to a selected subset, with
#'   the recommendation to include a consistent fraction of cell-identity genes,
#'   and possibly focusing on markers specific for the biological question of
#'   interest (for instance neural cortex layering markers). In this case we
#'   denote it as local differentiation index (LDI) relative to the selected
#'   subset.
#'
#' @param objCOTAN a `COTAN` object.
#' @param primaryMarkers A vector of primary marker names.
#' @param numGenesPerMarker The number of genes correlated with the primary
#'   markers that we want to consider. By default this is set to 25.
#'
#' @returns A `list` with:
#'   * the coex `data.frame`
#'   * a named `list` that for each secondary marker, gives the `list` of
#'     primary markers that selected for it
#'   * a `data.frame` with the rank of each gene according to its pValue
#'
#' @export
#'
#' @importFrom stats quantile
#' @importFrom stats pchisq
#'
#' @importFrom tibble rownames_to_column
#'
#' @importFrom Matrix rowMeans
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 12)
#' markers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
#' GCS <- genesCoexpressionSpace(objCOTAN,
#'                               primaryMarkers = markers)
#'                               numGenesPerMarker = 15)
#'
#' @rdname genesCoexSpace
#'
genesCoexSpace <-
  function(objCOTAN, primaryMarkers, numGenesPerMarker = 25) {
  # da sistemare: output tanh of reduced coex matrix
  logThis("Calculating gene coexpression space - START", logLevel = 2)

  {
    genesBelong <- primaryMarkers %in% getGenes(objCOTAN)
    if (!all(genesBelong)) {
      warning(paste0("Genes [",
                     paste0(primaryMarkers[!genesBelong], collapse = ", "),
                     "] are not present in the 'COTAN' object!"))
      primaryMarkers <- primaryMarkers[genesBelong]
    }
    rm(genesBelong)
  }

  pValue <- calculatePValue(objCOTAN, geneSubsetCol = primaryMarkers)

  secondaryMarkers <- primaryMarkers

  for (marker in primaryMarkers) {
    lowestPValueGenesPos <- order(pValue[, marker])[1:numGenesPerMarker]
    secondaryMarkers <- c(secondaryMarkers, getGenes(objCOTAN)[lowestPValueGenesPos])
  }
  # now the secondary markers are unique and sorted accordingly
  # to the original genes' list
  secondaryMarkers <- getGenes(objCOTAN)[getGenes(objCOTAN) %in% unique(secondaryMarkers)]

  logThis(paste("Number of selected secondary markers:", length(secondaryMarkers)), logLevel = 1)

  pValue <- pValue[secondaryMarkers, ]

  rankGenes <- data.frame()
  for (marker in primaryMarkers) {
    rankGenes <- setColumnInDF(rankGenes,
                               colToSet = rank(pValue[, marker], ties.method = "first"),
                               colName = marker, rowNames = secondaryMarkers)
  }

  top10pcCols <- as.integer(max(1, round(length(secondaryMarkers) / 10, digits = 0)))

  pValueSorted <- calculateS(objCOTAN, geneSubsetCol = secondaryMarkers)
  # put on the left the smalles pValues for each row
  pValueSorted <- t(apply(t(pValueSorted), 2, sort, decreasing = TRUE))
  pValueSorted <- pValueSorted[, 1:top10pcCols, drop = FALSE]
  pValueSorted <- pchisq(as.matrix(pValueSorted), df = 1, lower.tail = FALSE)

  LDI <- set_names(as.data.frame(rowMeans(pValueSorted)), "Local.GDI")

  rm(pValueSorted)

  lowLDIGenes <- LDI[["Local.GDI"]] <= quantile(LDI[["Local.GDI"]], probs = 0.1)
  goodGenes <- getGenes(objCOTAN) %in% rownames(LDI)[lowLDIGenes]

  rm(LDI)
  gc()

  GCS <- getGenesCoex(objCOTAN, genes = secondaryMarkers)[goodGenes, ]
  GCS <- tanh(GCS * sqrt(getNumCells(objCOTAN)))

  logThis(paste0("Number of columns (V set - secondary markers): ", dim(GCS)[2]), logLevel = 3)
  logThis(paste0("Number of rows (U set): ", dim(GCS)[1]), logLevel = 3)

  logThis("Calculating gene coexpression space - DONE", logLevel = 2)

  return(list("SecondaryMarkers" = secondaryMarkers, "GCS" = GCS,
              "rankGenes" = rankGenes))
}

#' establishGenesClusters
#'
#' @description This function perform the genes' clustering based on a pool of
#'   gene markers, using the genes' coexpression space
#'
#' @param objCOTAN a `COTAN` object
#' @param groupMarkers a named `list` with an element for each group of one or
#'   more marker genes for each group.
#' @param numMarkers the number of correlated genes to keep as other markers
#'   (default 25)
#' @param kCuts the number of estimated cluster (this defines the high for the
#'   tree cut)
#' @param distance type of distance to use (default "cosine"... "euclidean" is
#'   also available)
#' @param hclustMethod default is "ward.D2" but can be any method defined by
#'   [stats::hclust()] function
#'
#' @returns a `list` of
#'   1. the gene coexpression space dataframe
#'   2. the eigenvalues plot (using eigenvalue from factoextra package)
#'   3. the pca component dataframe
#'   4. the tree plot for the coexpression space genes
#'
#' @export
#'
#' @importFrom factoextra fviz_eig
#'
#' @importFrom dendextend set
#' @importFrom dendextend color_labels
#' @importFrom dendextend color_branches
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats prcomp
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#'
#' @importFrom stringr str_split
#'
#' @examples
#'
#' @rdname establishGenesClusters
#'
establishGenesClusters <-
  function(objCOTAN, groupMarkers, numGenesPerMarker = 25,
           kCuts = 6, distance = "cosine", hclustMethod = "ward.D2"){
  logThis("Establishing gene clusters - START", logLevel = 2)

  stopifnot("Group markers must have names for each group" <- !is.null(names(groupMarkers)))

  maxCol <- brewer.pal.info[c("Set2", "Set1", "Set3"), "maxcolors"]

  stopifnot("kCuts greater than number of possible supported colors" = kCuts <= sum(maxCol))

  # Dropping the genes not present
  filterGenes <- function(markers, genes) {
    markers[markers %in% genes]
  }
  groupMarkers <- lapply(groupMarkers, filterGenes, getGenes(objCOTAN))
  groupMarkers <- groupMarkers[!is.na(groupMarkers)]

  primaryMarkers = unlist(groupMarkers)

  list[secondaryMarkers, GCS, rankGenes] <-
    genesCoexSpace(objCOTAN,
                   numGenesPerMarker = numGenesPerMarker,
                   primaryMarkers = primaryMarkers)


  GCSPca <- prcomp(GCS, center = TRUE, scale. = FALSE)

  SMRelevance <- matrix(nrow = length(secondaryMarkers), ncol = length(groupMarkers),
                        dimnames = list(secondaryMarkers, names(groupMarkers)))

  for (name in names(groupMarkers)) {
    SMRelevance[ , name] <- rowSums(rankGenes[, groupMarkers[[name]], drop = FALSE])
    SMRelevance[groupMarkers[[name]], name] <- 1
  }

  tempCoex = getGenesCoex(objCOTAN, genes = primaryMarkers)[secondaryMarkers, ]

  for (name in names(groupMarkers)) {
    antiCorrelated <- tempCoex[, groupMarkers[[name]]] < 0
    toChange <- rownames(tempCoex)[rowSums(antiCorrelated) > 0]
    SMRelevance[toChange, name] <- 100000
  }

  pos.link  <- set_names(vector("list", length(groupMarkers)), names(groupMarkers))

  for (g in secondaryMarkers) {
    w <- which(SMRelevance[g, ] == min(SMRelevance[g, ]))
    if (length(w) != 1) {
      next
    }
    pos.link[[w]] = c(pos.link[[w]], g)
  }

  plotEigen <- fviz_eig(GCSPca, addlabels = TRUE, ncp = 10)

  if (distance == "cosine") {
    coexDist <- cosineDissimilarity(t(GCS))
  } else if(distance == "euclidean") {
    coexDist <- dist(GCS)
  } else {
    stop("only 'cosine' and 'euclidean' distances are supported")
  }
  hcNorm <- hclust(coexDist, method = hclustMethod)

  dend <- as.dendrogram(hcNorm)

  pca1 <- as.data.frame(GCSPca[["x"]][, 1:10])
  pca1 <- pca1[order.dendrogram(dend), ]

  highlight <- rep("not_marked", nrow(pca1))
  for (n in names(pos.link)) {
    highlight[rownames(pca1) %in% pos.link[[n]]] <- paste0("Genes related to ", n)
  }
  pca1[["highlight"]] <- highlight

  hClust <- cutree(hcNorm, k = kCuts)[rownames(pca1)]
  pca1[["hclust"]] <- hClust

  col_vector = brewer.pal(min(kCuts, maxCol[1]), name = "Set2")
  if (kCuts > maxCol[1]) {
    col_vector <- c(col_vector, brewer.pal(min(kCuts - maxCol[1], maxCol[2]),
                                           name = "Set1"))
  }
  if (kCuts > sum(maxCol[1:2])) {
    col_vector <- c(col_vector, brewer.pal(min(kCuts - sum(maxCol[1:2]), maxCol[3]),
                                           name = "Set3"))
  }

  pca1[["sec_markers"]] <- 0
  pca1[["sec_markers"]][rownames(pca1) %in% secondaryMarkers] <- 1

  colors = rep("#B09C85FF", nrow(pca1))
  c = 1
  for (to.color in unique(highlight)) {
    if (to.color == "not_marked") {
      next
    }
    colors[highlight == to.color] = col_vector[c]
    c <- c + 1
  }
  pca1[["colors"]] <- colors

  pca1 <- pca1[labels(dend),]

  col_branches <- rep("#B09C85FF", nrow(pca1))
  groupLabels  <- hClust
  for (cl in unique(hClust)) {
    clPos = (hClust == cl)
    for(groupName in names(groupMarkers)) {
      for (marker in groupMarkers[[groupName]]) {
        if (!marker %in% rownames(pca1[clPos,])) {
          next
        }
        col_branches[clPos] <- colors[rownames(pca1) %in% marker]
        groupLabels[clPos]  <- groupName
      }
    }
  }

  pca1[["col_branches"]] <- col_branches
  pca1[["groupLabels"]]  <- groupLabels

  uniquePos <- !duplicated(hClust)
  dend <- color_branches(dend, k = kCuts,
                         col         = col_branches[uniquePos],
                         groupLabels = groupLabels [uniquePos])

  dend <- color_labels(dend, labels = rownames(pca1), col = colors)

  logThis("Establishing gene clusters - DONE", logLevel = 2)

  return(list("g.space" = GCS, "plot.eig" = plotEigen,
              "pca_clusters" = pca1, "tree_plot" = dend))
}
