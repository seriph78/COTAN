
#' @title Calculations of genes statistics
#'
#' @description A collection of functions returning various statistics
#'   associated to the genes. In particular the *discrepancy* between the
#'   expected probabilities of zero and their actual occurrences, both at single
#'   gene level or looking at genes' pairs
#'
#' @description To make the `GDI` more specific, it may be desirable to restrict
#'   the set of genes against which `GDI` is computed to a selected subset, with
#'   the recommendation to include a consistent fraction of cell-identity genes,
#'   and possibly focusing on markers specific for the biological question of
#'   interest (for instance neural cortex layering markers). In this case we
#'   denote it as *Local Differentiation Index* (`LDI`) relative to the selected
#'   subset.
#'
#' @name GenesStatistics
NULL

#' @details `genesCoexSpace()` calculates genes groups based on the primary
#'   markers and uses them to prepare the genes' `COEX` space `data.frame`.
#'
#' @param objCOTAN a `COTAN` object.
#' @param primaryMarkers A vector of primary marker names.
#' @param numGenesPerMarker The number of genes correlated with the primary
#'   markers that we want to consider. By default this is set to 25.
#'
#' @returns `genesCoexSpace()` returns a `list` with:
#'  * `"SecondaryMarkers"` a named `list` that for each secondary marker,
#'    gives the `list` of primary markers that selected for it
#'  * `"GCS"` the relevant subset of `COEX` `matrix`
#'  * `"rankGenes"` a `data.frame` with the rank of each gene according to its
#'    *p-value*
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom stats quantile
#' @importFrom stats pchisq
#'
#' @importFrom tibble rownames_to_column
#'
#' @importFrom Rfast rowsums
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 6L, saveObj = FALSE)
#'
#' markers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
#' GCS <- genesCoexSpace(objCOTAN, primaryMarkers = markers,
#'                       numGenesPerMarker = 15)
#'
#' @rdname GenesStatistics
#'
genesCoexSpace <-
  function(objCOTAN, primaryMarkers, numGenesPerMarker = 25L) {
  logThis("Calculating gene co-expression space - START", logLevel = 2L)

  if (TRUE) {
    genesBelong <- primaryMarkers %in% getGenes(objCOTAN)
    if (!all(genesBelong)) {
      warning("Genes [", toString(primaryMarkers[!genesBelong]),
              "] are not present in the 'COTAN' object!")
      primaryMarkers <- primaryMarkers[genesBelong]
    }
    rm(genesBelong)
  }

  pValue <- calculatePValue(objCOTAN, geneSubsetCol = primaryMarkers)

  secondaryMarkers <- primaryMarkers

  for (marker in primaryMarkers) {
    lowestPValueGenesPos <-
      order(pValue[, marker, drop = TRUE])[1L:numGenesPerMarker]
    secondaryMarkers <- c(secondaryMarkers,
                          getGenes(objCOTAN)[lowestPValueGenesPos])
  }
  # now the secondary markers are unique and sorted accordingly
  # to the original genes' list
  secondaryMarkers <-
    getGenes(objCOTAN)[getGenes(objCOTAN) %in% unique(secondaryMarkers)]

  logThis(paste("Number of selected secondary markers:",
                length(secondaryMarkers)), logLevel = 1L)

  pValue <- pValue[secondaryMarkers, , drop = FALSE]

  rankGenes <- data.frame()
  for (marker in primaryMarkers) {
    rankGenes <- setColumnInDF(rankGenes,
                               colToSet = rank(pValue[, marker, drop = TRUE],
                                               ties.method = "first"),
                               colName = marker,
                               rowNames = secondaryMarkers)
  }
  rm(pValue)

  # put on the left the smalles log(-log(pValues)) for each row
  S <- calculateS(objCOTAN, geneSubsetRow = secondaryMarkers)
  LDI <- calculateGDIGivenS(S = S, rowsFraction = 0.1)

  rm(S)
  gc()

  lowLDIGenes <- LDI >= quantile(LDI, probs = 0.9)
  goodGenes <- getGenes(objCOTAN) %in% names(LDI)[lowLDIGenes]

  GCS <- as.matrix(getGenesCoex(objCOTAN,
                                genes = secondaryMarkers))[goodGenes, ]
  GCS <- tanh(GCS * sqrt(getNumCells(objCOTAN)))

  logThis(paste0("Number of columns (V set - secondary markers): ",
                 dim(GCS)[[2L]]), logLevel = 3L)
  logThis(paste0("Number of rows (U set): ", dim(GCS)[[1L]]), logLevel = 3L)

  logThis("Calculating gene co-expression space - DONE", logLevel = 2L)

  return(list("SecondaryMarkers" = secondaryMarkers, "GCS" = GCS,
              "rankGenes" = rankGenes))
}


#' @details `establishGenesClusters()` perform the genes' clustering based on a
#'   pool of gene markers, using the genes' `COEX` space

#' @param objCOTAN a `COTAN` object
#' @param groupMarkers a named `list` with an element for each group comprised
#'   of one or more marker genes
#' @param numGenesPerMarker the number of correlated genes to keep as other
#'   markers (default 25)
#' @param kCuts the number of estimated *cluster* (this defines the height for
#'   the tree cut)
#' @param distance type of distance to use. Default is `"cosine"`. Can be chosen
#'   among those supported by [parallelDist::parDist()]
#' @param hclustMethod default is "ward.D2" but can be any method defined by
#'   [stats::hclust()] function
#'
#' @returns `establishGenesClusters()` a `list` of:
#'   * `"g.space"` the genes' `COEX` space `data.frame`
#'   * `"plot.eig"` the eigenvalues plot
#'   * `"pca_clusters"` the *pca* components `data.frame`
#'   * `"tree_plot"` the tree plot for the genes' `COEX` space
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom dendextend set
#' @importFrom dendextend color_labels
#' @importFrom dendextend color_branches
#'
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#'
#' @importFrom parallelDist parDist
#'
#' @importFrom PCAtools pca
#' @importFrom PCAtools screeplot
#' @importFrom BiocSingular IrlbaParam
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @examples
#' groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
#'                      G2 = c("g-000300"),
#'                      G3 = c("g-000510", "g-000530", "g-000550",
#'                             "g-000570", "g-000590"))
#'
#' resList <-  establishGenesClusters(objCOTAN, groupMarkers = groupMarkers,
#'                                    numGenesPerMarker = 11)
#'
#' @rdname GenesStatistics
#'
establishGenesClusters <-
  function(objCOTAN, groupMarkers, numGenesPerMarker = 25L,
           kCuts = 6L, distance = "cosine", hclustMethod = "ward.D2") {
  logThis("Establishing gene clusters - START", logLevel = 2L)

  assert_that(!is.null(names(groupMarkers)),
              msg = "Group markers must have names for each group")

  colVector <- getColorsVector(kCuts)

  # Dropping the genes not present
  filterGenes <- function(markers, genes) {
    markers[markers %in% genes]
  }
  groupMarkers <- lapply(groupMarkers, filterGenes, getGenes(objCOTAN))
  groupMarkers <- groupMarkers[!is.na(groupMarkers)]

  primaryMarkers <- unlist(groupMarkers)

  c(secondaryMarkers, GCS, rankGenes) %<-%
    genesCoexSpace(objCOTAN,
                   numGenesPerMarker = numGenesPerMarker,
                   primaryMarkers = primaryMarkers)

  GCSPca <- pca(mat = GCS, rank = 10L,
                transposed = TRUE, BSPARAM = IrlbaParam())
  assert_that(identical(rownames(GCSPca[["rotated"]]), rownames(GCS)),
              (ncol(GCSPca[["rotated"]]) == 10L),
              msg = "Issues with pca output")

  SMRelevance <- matrix(nrow = length(secondaryMarkers),
                        ncol = length(groupMarkers),
                        dimnames = list(secondaryMarkers, names(groupMarkers)))

  for (name in names(groupMarkers)) {
    SMRelevance[, name] <- rowSums(rankGenes[, groupMarkers[[name]],
                                             drop = FALSE])
    SMRelevance[groupMarkers[[name]], name] <- 1L
  }

  tempCoex <-
    getGenesCoex(objCOTAN, genes = primaryMarkers)[secondaryMarkers, ,
                                                   drop = FALSE]

  for (name in names(groupMarkers)) {
    antiCorrelated <- tempCoex[, groupMarkers[[name]], drop = FALSE] < 0.0
    toChange <- rownames(tempCoex)[rowSums(antiCorrelated) > 0L]
    SMRelevance[toChange, name] <- 100000.0
  }

  posLink <- set_names(vector("list", length(groupMarkers)),
                        names(groupMarkers))

  for (g in secondaryMarkers) {
    w <- which(SMRelevance[g, ] == min(SMRelevance[g, ]))
    if (length(w) != 1L) {
      next
    }
    posLink[[w]] <- c(posLink[[w]], g)
  }

  plotEigen <- screeplot(GCSPca)

  coexDist <- parDist(as.matrix(GCS), method = distance)

  hcNorm <- hclust(coexDist, method = hclustMethod)

  dend <- as.dendrogram(hcNorm)

  pca1 <- as.data.frame(GCSPca[["rotated"]])
  pca1 <- pca1[order.dendrogram(dend), ]

  highlight <- rep("not_marked", nrow(pca1))
  for (n in names(posLink)) {
    highlight[rownames(pca1) %in% posLink[[n]]] <-
      paste0("Genes related to ", n)
  }
  pca1[["highlight"]] <- highlight

  hClust <- cutree(hcNorm, k = kCuts)[rownames(pca1)]
  pca1[["hclust"]] <- hClust

  pca1[["sec_markers"]] <- 0.0
  pca1[["sec_markers"]][rownames(pca1) %in% secondaryMarkers] <- 1.0

  colors <- rep("#B09C85FF", nrow(pca1))
  c <- 1L
  for (to.color in unique(highlight)) {
    if (to.color == "not_marked") {
      next
    }
    colors[highlight == to.color] <- colVector[c]
    c <- c + 1L
  }
  pca1[["colors"]] <- colors

  pca1 <- pca1[labels(dend), ]

  colBranches <- rep("#B09C85FF", nrow(pca1))
  groupLabels  <- hClust
  for (cl in unique(hClust)) {
    clPos <- (hClust == cl)
    for (groupName in names(groupMarkers)) {
      for (marker in groupMarkers[[groupName]]) {
        if (!marker %in% rownames(pca1[clPos, ])) {
          next
        }
        colBranches[clPos] <- colors[rownames(pca1) %in% marker]
        groupLabels[clPos] <- groupName
      }
    }
  }

  pca1[["col_branches"]] <- colBranches
  pca1[["groupLabels"]]  <- groupLabels

  uniquePos <- !duplicated(hClust)
  dend <- color_branches(dend, k = kCuts,
                         col         = colBranches[uniquePos],
                         groupLabels = groupLabels[uniquePos])

  dend <- color_labels(dend, labels = rownames(pca1), col = colors)

  relPos <- rownames(pca1) %in% colnames(GCS)
  dend <- dend %>%
    dendextend::set("labels",
                    ifelse(labels(dend) %in% rownames(pca1)[relPos],
                           labels(dend), ""))

  logThis("Establishing gene clusters - DONE", logLevel = 2L)

  return(list("g.space" = GCS, "plot.eig" = plotEigen,
              "pca_clusters" = pca1, "tree_plot" = dend))
}
