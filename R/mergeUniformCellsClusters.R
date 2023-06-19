
# --------------------- Uniform Clusters ----------------------

#' @details `mergeUniformCellsClusters()` takes in a **uniform**
#'   *clusterization* and iteratively checks whether merging two *near clusters*
#'   would form a **uniform** *cluster* still. This function uses the *cosine
#'   distance* to establish the *nearest clusters pairs*. It will use the
#'   [checkClusterUniformity()] function to check whether the merged *clusters*
#'   are **uniform**. The function will stop once no *near pairs* of clusters
#'   are mergeable in a single batch
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters The *clusterization* to merge. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param GDIThreshold the threshold level that discriminates uniform clusters.
#'   It defaults to \eqn{1.4}
#' @param batchSize Number pairs to test in a single round. If none of them
#'   succeeds the merge stops
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
#' @returns a `list` with:
#'   * `"clusters"` the merged cluster labels array
#'   * `"coex"` the associated `COEX` `data.frame`
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom Matrix t
#'
#' @importFrom parallelDist parDist
#'
#' @importFrom stats hclust
#' @importFrom stats as.dendrogram
#' @importFrom stats cophenetic
#'
#' @importFrom stringr str_split
#' @importFrom stringr fixed
#'
#' @importFrom dendextend get_nodes_attr
#'
#' @importFrom zeallot `%<-%`
#' @importFrom zeallot `%->%`
#'
#' @examples
#' data("test.dataset")
#'
#' objCOTAN <- automaticCOTANObjectCreation(raw = test.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12L,
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
#'
#' splitList <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                     GDIThreshold = 1.5, saveObj = FALSE)
#'
#' clusters <- splitList[["clusters"]]
#'
#' firstCluster <- getCells(objCOTAN)[clusters %in% clusters[[1L]]]
#' checkClusterUniformity(objCOTAN,
#'                        GDIThreshold = 1.5,
#'                        cluster = clusters[[1L]],
#'                        cells = firstCluster,
#'                        cores = 12L,
#'                        saveObj = FALSE)
#'
#' objCOTAN <- addClusterization(objCOTAN,
#'                               clName = "split",
#'                               clusters = clusters)
#'
#' objCOTAN <- addClusterizationCoex(objCOTAN,
#'                                   clName = "split",
#'                                   coexDF = splitList[["coex"]])
#'
#' expect_identical(reorderClusterization(objCOTAN)[["clusters"]], clusters)
#'
#' mergedList <- mergeUniformCellsClusters(objCOTAN,
#'                                         GDIThreshold = 1.5,
#'                                         batchSize = 5L,
#'                                         clusters = clusters,
#'                                         cores = 12L,
#'                                         distance = "cosine",
#'                                         hclustMethod = "ward.D2",
#'                                         saveObj = FALSE)
#'
#' objCOTAN <- addClusterization(objCOTAN,
#'                               clName = "merged",
#'                               clusters = mergedList[["clusters"]],
#'                               coexDF = mergedList[["coex"]])
#'
#' @rdname UniformClusters
#'

mergeUniformCellsClusters <- function(objCOTAN,
                                      clusters = NULL,
                                      GDIThreshold = 1.4,
                                      batchSize = 10L,
                                      cores = 1L,
                                      distance = "cosine",
                                      hclustMethod = "ward.D2",
                                      saveObj = TRUE,
                                      outDir = ".") {
  logThis("Merging cells' uniform clustering: START", logLevel = 2L)

  outputClusters <- clusters
  if (is_empty(outputClusters)) {
    # pick the last clusterization
    outputClusters <- getClusters(objCOTAN)
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

  mergedName <- function(cl1, cl2) {
    return(paste0(min(cl1, cl2), "_", max(cl1, cl2), "-merge"))
  }

  testPairListMerge <- function(pList) {
    logThis(paste0("Clusters pairs for merging:\n",
                   paste(pList, collapse = " ")), logLevel = 1L)

    for (p in pList) {
      logThis("*", logLevel = 1L, appendLF = FALSE)

      c(cl1, cl2) %<-% p

      if (all(outputClusters != cl1) || all(outputClusters != cl2)) {
        logThis(paste0("Clusters ", cl1, " or ", cl2,
                       " is now missing due to previous merges: skip."),
                logLevel = 3L)
        next
      }

      mergedClName <- mergedName(cl1, cl2)

      logThis(mergedClName, logLevel = 3L)

      if (mergedClName %in% notMergeable) {
        logThis(paste0("Clusters ", cl1, " and ", cl2,
                       " already analyzed and not mergeable: skip."),
                logLevel = 3L)
        next
      }

      mergedCluster <- names(outputClusters)[outputClusters %in% c(cl1, cl2)]

      clusterIsUniform <-
        checkClusterUniformity(objCOTAN,
                               cluster = mergedClName,
                               cells = mergedCluster,
                               GDIThreshold = GDIThreshold,
                               cores = cores,
                               saveObj = saveObj,
                               outDir = mergeOutDir)[["isUniform"]]

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

    return(list("oc" = outputClusters, "nm" = notMergeable))
  }

  iter <- 0L
  repeat {
    iter <- iter + 1L
    logThis(paste0("Start merging nearest clusters: iteration ", iter),
            logLevel = 3L)

    oldNumClusters <- length(unique(outputClusters))

    coexDF <- DEAOnClusters(objCOTAN, clusters = outputClusters)
    gc()

    ## To drop the cells with only zeros
    ## TODO: To be fixed! Where it come from?
    coexDF <- coexDF[, colSums(coexDF != 0.0) > 0L]

    # merge small cluster based on their DEA based distances
    coexDist <- parDist(t(as.matrix(coexDF)), method = distance)

    if (saveObj) {
      pdf(file.path(mergeOutDir, paste0("dend_iter_", iter, "_plot.pdf")))

      hcNorm <- hclust(coexDist, method = hclustMethod)
      plot(as.dendrogram(hcNorm))

      dev.off()
    }

    # We will check whether it is possible to merge a list of cluster pairs.
    # These pairs correspond to N lowest distances as calculated before
    # If none of them can be merges, the loop stops

    allLabels <- colnames(coexDF)

    # create all pairings with different clusters
    pList <- rbind(rep((1L:oldNumClusters), each  = oldNumClusters),
                   rep((1L:oldNumClusters), times = oldNumClusters))
    pList <- pList[, pList[1L, ] < pList[2L, ], drop = FALSE]
    pList <- matrix(allLabels[pList], nrow = 2L)

    # reorder the pairings using the distance and pick only those necessary
    pList <- as.list(as.data.frame(pList))

    # reorder based on distance
    pList <- pList[order(coexDist)]

    # drop the already tested pairs
    pNamesList <- lapply(pList, function(p) mergedName(p[[1L]], p[[2L]]))
    pList <- pList[!pNamesList %in% notMergeable]

    # take the first N remaining
    numPairsToTest <- min(batchSize, length(pList))
    pList <- pList[1L:numPairsToTest]

    c(outputClusters, notMergeable) %<-% testPairListMerge(pList)

    newNumClusters <- length(unique(outputClusters))
    if (newNumClusters == 1L) {
      # nothing left to do: stop!
      break
    }

    if (newNumClusters == oldNumClusters) {
      logThis(msg = paste("None of the", numPairsToTest,
                          "nearest cluster pairs could be merged"),
              logLevel = 3L)

      # No merges happened -> too low probability of new merges...
      break
    } else {
      logThis(paste0("Executed ", (oldNumClusters - newNumClusters),
                     " merges out of ", numPairsToTest), logLevel = 3L)
    }
  }

  logThis(paste0("The final merged clusterization contains [",
                 length(unique(outputClusters)), "] different clusters: ",
                 toString(sort(unique(outputClusters)))), logLevel = 1L)

  # replace the clusters' tags with completely new ones
  {
    clTags <- sort(unique(outputClusters))

    clTagsMap <- paste0(seq_along(clTags))
    clTagsMap <- factorToVector(niceFactorLevels(clTagsMap))
    clTagsMap <- set_names(clTagsMap, clTags)

    outputClusters <- clTagsMap[outputClusters]
    outputClusters <- set_names(outputClusters, getCells(objCOTAN))

    colnames(coexDF)   <- clTagsMap[colnames(coexDF)]
  }

  c(outputClusters, coexDF) %<-%
    reorderClusterization(objCOTAN, clusters = outputClusters, coexDF = coexDF,
                          reverse = FALSE, keepMinusOne = FALSE,
                          distance = distance, hclustMethod = hclustMethod)

  logThis("Merging cells' uniform clustering: DONE", logLevel = 2L)

  return(list("clusters" = outputClusters, "coex" = coexDF))
}
