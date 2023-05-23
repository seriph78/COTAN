
# --------------------- Uniform Clusters ----------------------

#' @details `mergeUniformCellsClusters()` takes in a **uniform**
#'   *clusterization* and iteratively checks whether merging two *near clusters*
#'   would form a **uniform** *cluster* still. This function uses the *cosine
#'   distance* and the [stats::hclust()] function to establish *near clusters
#'   pairs*. It will use the [checkClusterUniformity()] function to check
#'   whether the merged *clusters* are **uniform**. The function will stop once
#'   no *near pairs* of clusters are mergeable.
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters The *clusterization* to merge. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param GDIThreshold the threshold level that discriminates uniform clusters.
#'   It defaults to \eqn{1.4}
#' @param cores number cores used
#' @param distance type of distance to use (default is `"cosine"`, `"euclidean"`
#'   and the others from [parallelDist::parDist()] are also available)
#' @param hclustMethod It defaults is `"ward.D2"` but can be any of the methods
#'   defined by the [stats::hclust()] function.
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output. The effective
#'   output will be paced in a sub-folder.
#'
#' @returns a `list` with "clusters", "coexDF" and "pValueDF"
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom Matrix t
#'
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats as.dendrogram
#'
#' @importFrom stringr str_split
#' @importFrom stringr fixed
#'
#' @importFrom dendextend get_nodes_attr
#'
#' @examples
#' data("test.dataset")
#'
#' objCOTAN <- automaticCOTANObjectCreation(raw = test.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = FALSE)
#'
#' groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
#'                      G2 = c("g-000300", "g-000330"),
#'                      G3 = c("g-000510", "g-000530", "g-000550",
#'                             "g-000570", "g-000590"))
#' gdiPlot <- GDIPlot(objCOTAN, genes = groupMarkers, cond = "test")
#' plot(gdiPlot)
#'
#' ## Here we override the default GDI threshold as a way to speed-up
#' ## calculations as higher threshold implies less stringent uniformity
#' ## It real applications it might be appropriate to change the threshold
#' ## in cases of relatively low genes/cells number, or in cases when an
#' ## rough clusterization is needed in the early satges of the analysis
#' ##
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    GDIThreshold = 1.5, saveObj = FALSE)
#'
#' checkClusterUniformity(objCOTAN, GDIThreshold = 1.5,
#'                        cluster = clusters[1],
#'                        cells = getCells(objCOTAN)[clusters %in% clusters[1]],
#'                        cores = 12,
#'                        saveObj = FALSE)
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "uniformClusters",
#'                               clusters = clusters)
#'
#' mergedList <- mergeUniformCellsClusters(objCOTAN, GDIThreshold = 1.5,
#'                                         clusters = clusters,
#'                                         cores = 12,
#'                                         distance = "cosine",
#'                                         hclustMethod = "ward.D2",
#'                                         saveObj = FALSE)
#'
#' objCOTAN <- addClusterization(objCOTAN, clName = "mergedUniformClusters",
#'                               clusters = mergedList[["clusters"]],
#'                               coexDF = mergedList[["coexDF"]])
#'
#' @rdname UniformClusters
#'

mergeUniformCellsClusters <- function(objCOTAN,
                                      clusters = NULL,
                                      GDIThreshold = 1.4,
                                      cores = 1L,
                                      distance = "cosine",
                                      hclustMethod = "ward.D2",
                                      saveObj = TRUE,
                                      outDir = ".") {
  logThis("Merging cells' uniform clustering: START", logLevel = 2L)

  outputClusters <- clusters
  if (is_empty(outputClusters)) {
    # pick the last clusterization
    outputClusters <- getClusterizationData(objCOTAN)[["clusters"]]
  }

  outputClusters <- factorToVector(outputClusters)

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  outDirCond <- file.path(outDir, cond)
  if (!file.exists(outDirCond)) {
    dir.create(outDirCond)
  }

  mergeOutDir <- file.path(outDirCond, "leafs_merge")
  if (isTRUE(saveObj) && !file.exists(mergeOutDir)) {
    dir.create(mergeOutDir)
  }

  notMergeable <- vector(mode = "character")
  iter <- 0L
  repeat {
    iter <- iter + 1L
    logThis(paste0("Start merging smallest clusters: iteration ", iter),
            logLevel = 3L)

    oldNumClusters <- length(unique(outputClusters))

    c(coexDF, pValueDF) %<-% DEAOnClusters(objCOTAN, clusters = outputClusters)
    gc()

    ## To drop the cells with only zeros
    ## TODO: To be fixed! Where it come from?
    coexDF <- coexDF[, colSums(coexDF != 0.0) > 0L]

    # merge small cluster based on distances
    coexDist <- parDist(t(as.matrix(coexDF)), method = distance)

    hcNorm <- hclust(coexDist, method = hclustMethod)

    dend <- as.dendrogram(hcNorm)

    # This checks if any little two pair of leaf clusters of the dendogram
    # could be merged

    id <- NULL
    {
      members <- get_nodes_attr(dend, "members")
      for (i in seq_along(members)) {
        if (members[i] == 2L) {
          id <- c(id, i + 1L, i + 2L)
        }
      }
    }

    logThis(paste0("Created leafs ID for merging: ",
                   paste(get_nodes_attr(dend, "label", id = id),
                         collapse = " ")), logLevel = 2L)

    p <- 1L
    while (p < length(id)) {
      logThis("*", logLevel = 1L, appendLF = FALSE)

      cl1 <- get_nodes_attr(dend, "label", id = id[p + 0L])
      cl2 <- get_nodes_attr(dend, "label", id = id[p + 1L])

      mergedClName <- paste0(min(cl1, cl2), "_", max(cl1, cl2), "-merge")

      logThis(mergedClName, logLevel = 3L)

      p  <- p + 2L

      if (mergedClName %in% notMergeable) {
        logThis(paste0("Clusters ", cl1, " and ", cl2,
                       " already analyzed and not mergeable: skip."),
                logLevel = 3L)
        next
      }

      mergedCluster <- names(outputClusters)[outputClusters %in% c(cl1, cl2)]

      clusterIsUniform <- checkClusterUniformity(objCOTAN,
                                                 cluster = mergedClName,
                                                 cells = mergedCluster,
                                                 GDIThreshold = GDIThreshold,
                                                 cores = cores,
                                                 saveObj = saveObj,
                                                 outDir = mergeOutDir)

      gc()

      if (!clusterIsUniform) {
        logThis(paste("Merging clusters", cl1, "and", cl2,
                      "results in a too high GDI"), logLevel = 1L)

        notMergeable <- c(notMergeable, mergedClName)
      } else {
        logThis(paste("Clusters", cl1, "and", cl2, "can be merged"),
                logLevel = 1L)

        outputClusters[mergedCluster] <- mergedClName
      }
    }

    if (length(unique(outputClusters)) == oldNumClusters) {
      # no merges happened: stop!
      break
    }
  }

  logThis(paste("The final merged clusterization contains [",
                length(unique(outputClusters)), "] different clusters:",
                toString(unique(sort(outputClusters)))), logLevel = 3L)

  # replace the clusters' tags
  # non merged clusters keep their tag, while merged one use the lesser one!
  # there might be holes in the enumeration!
  {
    clTags <- levels(factor(outputClusters))
    # here we use the fact that the merged clusters have
    # the 'smaller' cluster as the first in the pair
    clTagsMap <- str_split(clTags, pattern = fixed("_"), simplify = TRUE)[, 1L]
    clTagsMap <- set_names(clTagsMap, clTags)

    outputClusters <- clTagsMap[outputClusters]
    outputClusters <- set_names(outputClusters, getCells(objCOTAN))

    colnames(coexDF)   <- clTagsMap[colnames(coexDF)]
    colnames(pValueDF) <- clTagsMap[colnames(pValueDF)]
  }

  outputClusters <- factor(outputClusters)

  logThis("Merging cells' uniform clustering: DONE", logLevel = 2L)

  return(list("clusters" = outputClusters,
              "coexDF" = coexDF, "pValueDF" = pValueDF))
}
