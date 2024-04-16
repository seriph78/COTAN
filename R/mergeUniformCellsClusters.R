
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
#'   It defaults to \eqn{1.43}
#' @param batchSize Number pairs to test in a single round. If none of them
#'   succeeds the merge stops
#' @param notMergeable An array of names of merged clusters that are already
#'   known for not being uniform. Useful to restart the *merging* process after
#'   an interruption.
#' @param cores number of cores to use. Default is 1.
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
#' @param useDEA Boolean indicating whether to use the *DEA* to define the
#'   distance; alternatively it will use the average *Zero-One* counts, that is
#'   faster but less precise.
#' @param distance type of distance to use. Default is `"cosine"` for *DEA* and
#'   `"euclidean"` for *Zero-One*. Can be chosen among those supported by
#'   [parallelDist::parDist()]
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
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @examples
#' data("test.dataset")
#'
#' objCOTAN <- automaticCOTANObjectCreation(raw = test.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 6L,
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
#' splitList <- cellsUniformClustering(objCOTAN, cores = 6L,
#'                                     optimizeForSpeed = TRUE,
#'                                     deviceStr = "cuda",
#'                                     initialResolution = 0.8,
#'                                     GDIThreshold = 1.46, saveObj = FALSE)
#'
#' clusters <- splitList[["clusters"]]
#'
#' firstCluster <- getCells(objCOTAN)[clusters %in% clusters[[1L]]]
#' checkClusterUniformity(objCOTAN,
#'                        GDIThreshold = 1.46,
#'                        cluster = clusters[[1L]],
#'                        cells = firstCluster,
#'                        cores = 6L,
#'                        optimizeForSpeed = TRUE,
#'                        deviceStr = "cuda",
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
#' identical(reorderClusterization(objCOTAN)[["clusters"]], clusters)
#'
#' mergedList <- mergeUniformCellsClusters(objCOTAN,
#'                                         GDIThreshold = 1.46,
#'                                         batchSize = 5L,
#'                                         clusters = clusters,
#'                                         cores = 6L,
#'                                         optimizeForSpeed = TRUE,
#'                                         deviceStr = "cpu",
#'                                         distance = "cosine",
#'                                         hclustMethod = "ward.D2",
#'                                         saveObj = FALSE)
#'
#' objCOTAN <- addClusterization(objCOTAN,
#'                               clName = "merged",
#'                               clusters = mergedList[["clusters"]],
#'                               coexDF = mergedList[["coex"]])
#'
#' identical(reorderClusterization(objCOTAN), mergedList)
#'
#' @rdname UniformClusters
#'

mergeUniformCellsClusters <- function(objCOTAN,
                                      clusters = NULL,
                                      GDIThreshold = 1.43,
                                      batchSize = 10L,
                                      notMergeable = NULL,
                                      cores = 1L,
                                      optimizeForSpeed = TRUE,
                                      deviceStr = "cuda",
                                      useDEA = TRUE,
                                      distance = NULL,
                                      hclustMethod = "ward.D2",
                                      saveObj = TRUE,
                                      outDir = ".") {
  logThis("Merging cells' uniform clustering: START", logLevel = 2L)

  assert_that(estimatorsAreReady(objCOTAN),
              msg = paste("Estimators lambda, nu, dispersion are not ready:",
                          "Use proceeedToCoex() to prepare them"))

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

  if (is_empty(notMergeable)) {
    notMergeable <- vector(mode = "character")
  }

  mergedName <- function(cl1, cl2) {
    return(paste0(min(cl1, cl2), "_", max(cl1, cl2), "-merge"))
  }

  iter <- 0L
  allCheckResults <- data.frame()
  errorCheckResults <-
    list("isUniform" = FALSE, "fractionAbove" = NA, "firstPercentile" = NA)


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

      checkResults <- tryCatch(
        checkClusterUniformity(objCOTAN,
                               clusterName = mergedClName,
                               cells = mergedCluster,
                               GDIThreshold = GDIThreshold,
                               cores = cores,
                               optimizeForSpeed = optimizeForSpeed,
                               deviceStr = deviceStr,
                               saveObj = saveObj,
                               outDir = mergeOutDir),
        error = function(err) {
          logThis(paste("While checking cluster uniformity", err),
                  logLevel = 0L)
          logThis("Marking pair as not mergable", logLevel = 1L)
          return(errorCheckResults)
        })

      gc()

      allCheckResults <- rbind(allCheckResults, checkResults)
      rownames(allCheckResults)[[nrow(allCheckResults)]] <- mergedClName

      if (!checkResults[["isUniform"]]) {
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

  repeat {
    iter <- iter + 1L
    logThis(paste0("Start merging nearest clusters: iteration ", iter),
            logLevel = 3L)

    oldNumClusters <- length(unique(outputClusters))

    clDist <- distancesBetweenClusters(objCOTAN, clusters = outputClusters,
                                       useDEA = useDEA, cores = cores,
                                       distance = distance)
    gc()

    if (isTRUE(saveObj)) tryCatch({
        pdf(file.path(mergeOutDir, paste0("dend_iter_", iter, "_plot.pdf")))

        hcNorm <- hclust(clDist, method = hclustMethod)
        plot(as.dendrogram(hcNorm))

        dev.off()
      },
      error = function(err) {
        logThis(paste("While saving dendogram plot", err), logLevel = 0L)
      }
    )

    # We will check whether it is possible to merge a list of cluster pairs.
    # These pairs correspond to N lowest distances as calculated before
    # If none of them can be merges, the loop stops

    allLabels <- labels(clDist)
    assert_that(length(allLabels) == oldNumClusters,
                msg = "Internal error - distance has no labels")

    # create all pairings with different clusters
    pList <- rbind(rep((1L:oldNumClusters), each  = oldNumClusters),
                   rep((1L:oldNumClusters), times = oldNumClusters))
    pList <- pList[, pList[1L, ] < pList[2L, ], drop = FALSE]
    pList <- matrix(allLabels[pList], nrow = 2L)

    # reorder the pairings using the distance and pick only those necessary
    pList <- as.list(as.data.frame(pList))

    # reorder based on distance
    pList <- pList[order(clDist)]

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

    if (isTRUE(saveObj)) tryCatch({
        outFile <- file.path(mergeOutDir,
                             paste0("merge_clusterization_", iter, ".csv"))
        write.csv(outputClusters, file = outFile)

        outFile <- file.path(mergeOutDir,
                             paste0("non_mergeable_clusters_", iter, ".csv"))
        write.csv(notMergeable, file = outFile)

        outFile <- file.path(mergeOutDir,
                             paste0("all_check_results_", iter, ".csv"))
        write.csv(allCheckResults, file = outFile)
      },
      error = function(err) {
        logThis(paste("While saving current clusterization", err),
                logLevel = 0L)
      }
    )
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

    allCheckResults <- allCheckResults[clTags, , drop = FALSE]
    allCheckResults <- allCheckResults[clTagsMap[clTags], , drop = FALSE]
  }

  outputCoexDF <-
    tryCatch(DEAOnClusters(objCOTAN, clusters = outputClusters, cores = cores),
             error = function(err) {
               logThis(paste("Calling DEAOnClusters", err), logLevel = 0L)
               return(NULL)
             })

  c(outputClusters, outputCoexDF) %<-% tryCatch(
    reorderClusterization(objCOTAN, clusters = outputClusters,
                          coexDF = outputCoexDF, reverse = FALSE,
                          keepMinusOne = FALSE, useDEA = useDEA, cores = cores,
                          distance = distance, hclustMethod = hclustMethod),
    error = function(err) {
      logThis(paste("Calling reorderClusterization", err), logLevel = 0L)
      return(list(outputClusters, outputCoexDF))
    })

  if (isTRUE(saveObj)) tryCatch({
    outFile <- file.path(outDirCond, "merge_check_results.csv")
    write.csv(allCheckResults, file = outFile)
  })

  logThis("Merging cells' uniform clustering: DONE", logLevel = 2L)

  return(list("clusters" = outputClusters, "coex" = outputCoexDF))
}
