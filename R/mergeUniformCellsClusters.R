
# --------------------- Uniform Clusters ----------------------

#' @details `mergeUniformCellsClusters()` takes in a **uniform**
#'   *clusterization* and iteratively checks whether merging two *near clusters*
#'   would form a **uniform** *cluster* still. Multiple thresholds will be used
#'   from \eqn{1.37} up to the given one in order to prioritize merge of the
#'   best fitting pairs.
#'
#'   This function uses the *cosine distance* to establish the *nearest clusters
#'   pairs*. It will use the [checkClusterUniformity()] function to check
#'   whether the merged *clusters* are **uniform**. The function will stop once
#'   no *tested pairs* of clusters are mergeable after testing all pairs in a
#'   single batch
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters The *clusterization* to merge. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param checkers a `list` of objects that defines the method and the
#'   *increasing* thresholds to discriminate whether to merge two *clusters* if
#'   deemed *uniform transcript*. See [UniformTranscriptCheckers] for more
#'   details
#' @param GDIThreshold legacy. The threshold level that is used in a
#'   [SimpleGDIUniformityCheck-class]. It defaults to \eqn{1.43}
#' @param batchSize Number pairs to test in a single round. If none of them
#'   succeeds the merge stops. Defaults to \eqn{2 (\#cl)^{2/3}}
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
#' @param allCheckResults An optional `data.frame` with the results of previous
#'   checks about the merging of clusters. Useful to restart the *merging*
#'   process after an interruption.
#' @param initialIteration the number associated tot he first iteration; it
#'   defaults to 1. Useful in case of restart of the procedure to avoid
#'   intermediate data override
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
#' @importFrom rlang is_named2
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.cur
#'
#' @importFrom Matrix t
#'
#' @importFrom stats hclust
#' @importFrom stats as.dendrogram
#' @importFrom stats cophenetic
#'
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#'
#' @importFrom dendextend get_nodes_attr
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom methods new
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
#'
#' gdiPlot <- GDIPlot(objCOTAN, genes = groupMarkers, cond = "test")
#' plot(gdiPlot)
#'
#' ## Here we override the default checker as a way to reduce the number of
#' ## clusters as higher thresholds imply less stringent uniformity checks
#' ##
#' ## In real applications it might be appropriate to do so in the cases when
#' ## the wanted resolution is lower such as in the early stages of the analysis
#' ##
#'
#' checker <- new("AdvancedGDIUniformityCheck")
#' identical(checker@firstCheck@GDIThreshold, 1.297)
#'
#' checker2 <- shiftCheckerThresholds(checker, 0.1)
#' identical(checker2@firstCheck@GDIThreshold, 1.397)
#'
#' splitList <- cellsUniformClustering(objCOTAN, cores = 6L,
#'                                     optimizeForSpeed = TRUE,
#'                                     deviceStr = "cuda",
#'                                     initialResolution = 0.8,
#'                                     checker = checker2,
#'                                     saveObj = FALSE)
#'
#' clusters <- splitList[["clusters"]]
#'
#' firstCluster <- getCells(objCOTAN)[clusters %in% clusters[[1L]]]
#'
#' checkerRes <-
#'   checkClusterUniformity(objCOTAN, checker = checker2,
#'                          cluster = clusters[[1L]], cells = firstCluster,
#'                          cores = 6L, optimizeForSpeed = TRUE,
#'                          deviceStr = "cuda", saveObj = FALSE)
#'
#' objCOTAN <- addClusterization(objCOTAN,
#'                               clName = "split",
#'                               clusters = clusters,
#'                               coexDF = splitList[["coex"]],
#'                               override = FALSE)
#'
#' identical(reorderClusterization(objCOTAN)[["clusters"]], clusters)
#'
#' ## It is possible to pass a list of checkers tot the merge function that will
#' ## be applied each to the *resulting* merged *clusterization* obtained using
#' ## the previous checker. This ensures that the most similar clusters are
#' ## merged first improving the overall performance
#'
#' mergedList <- mergeUniformCellsClusters(objCOTAN,
#'                                         checkers = c(checker, checker2),
#'                                         batchSize = 2L,
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
#'                               coexDF = mergedList[["coex"]],
#'                               override = TRUE)
#'
#' identical(reorderClusterization(objCOTAN), mergedList[["clusters"]])
#'
#' @rdname UniformClusters
#'

mergeUniformCellsClusters <- function(objCOTAN,
                                      clusters = NULL,
                                      checkers = NULL,
                                      GDIThreshold = NaN,
                                      batchSize = 0L,
                                      cores = 1L,
                                      optimizeForSpeed = TRUE,
                                      deviceStr = "cuda",
                                      useDEA = TRUE,
                                      distance = NULL,
                                      hclustMethod = "ward.D2",
                                      allCheckResults = data.frame(),
                                      initialIteration = 1L,
                                      saveObj = TRUE,
                                      outDir = ".") {
  # returns merged name given underlying names
  toMergedName <- function(clName1, clName2) {
    return(paste0(min(clName1, clName2), "_", max(clName1, clName2), "-merge"))
  }

  # returns underlying names given merged name
  fromMergedName <- function(mergedClName, currentClNames) {
    partialMatch <- vapply(currentClNames, function(clName, mergedName) {
      return(str_detect(mergedName, clName))
    }, FUN.VALUE = logical(1L), mergedClName)

    c(clName1, clName2) %<-% c("", "")
    if (sum(partialMatch) == 2L) {
      c(clName1, clName2) %<-% currentClNames[partialMatch]
    }

    return(c("clName1" = clName1, "clName2" = clName2))
  }

  # updates checker results with thresholds from current checker
  updateChecker <- function(mergedClName, checker, allCheckResults) {
    if (all(names(allCheckResults) != mergedClName)) {
      return(checker)
    } else {
      return(checkObjIsUniform(currentC = checker,
                               previousC = allCheckResults[[mergedClName]],
                               objCOTAN = NULL))
    }
  }

  # select all clusters pairings to be checked for this batch
  selectPairsList <- function(pList, batchSize, checker, allCheckResults) {
    # drop the already tested pairs
    pNamesList <- lapply(pList, function(p) toMergedName(p[[1L]], p[[2L]]))

    untestedPairs <-
      vapply(pNamesList, function(pName, check, allRes) {
        # we exploit the convention that a non-avalable result is marked by a
        # zero clusterSize
        updateChecker(pName, check, allRes)@clusterSize == 0L
      }, FUN.VALUE = logical(1L), checker, allCheckResults)

    pList <- pList[untestedPairs]

    # take the first N remaining
    numPairsToTest <- min(batchSize, length(pList))
    return(pList[seq_len(numPairsToTest)])
  }

  # runs the UT check on all pairings of the batch
  testPairListMerge <- function(pList, outputClusters,
                                checker, allCheckResults) {
    logThis(paste("Updating check results for the", length(allCheckResults),
                  "already tested pairs to align to new checker"),
            logLevel = 1L)

    if (TRUE) {
      allNames <- names(allCheckResults)
      allCheckResults <-
        lapply(seq_along(allCheckResults),
               function(r, check, allRes) {
                 updateChecker(names(allCheckResults)[[r]], check, allRes)
               }, checker, allCheckResults)
      names(allCheckResults) <- allNames
      rm(allNames)
    }

    logThis(paste0(length(pList), " new clusters pairs to be tested for",
                   " merging:\n", paste(pList, collapse = " ")), logLevel = 1L)

    for (p in pList) {
      logThis("*", logLevel = 1L, appendLF = FALSE)

      c(cl1, cl2) %<-% p

      if (!(any(outputClusters %in% cl1) && any(outputClusters %in% cl2))) {
        logThis(paste0("Clusters ", cl1, " or ", cl2,
                       " is now missing due to previous merges: skip"),
                logLevel = 3L)
        next
      }

      mergedClName <- toMergedName(cl1, cl2)

      logThis(mergedClName, logLevel = 3L)

      checkResults <- updateChecker(mergedClName, checker, allCheckResults)
      if (checkResults@clusterSize != 0L) {
        logThis(paste("Clusters", cl1, "and", cl2, "already analyzed: skip"),
                logLevel = 3L)
        next
      }
      # else we have insufficient information about the pair [re]calculate

      mergedCluster <- names(outputClusters)[outputClusters %in% c(cl1, cl2)]

      checkResults <- tryCatch(
        checkClusterUniformity(objCOTAN,
                               clusterName = mergedClName,
                               cells = mergedCluster,
                               checker = checker,
                               cores = cores,
                               optimizeForSpeed = optimizeForSpeed,
                               deviceStr = deviceStr,
                               saveObj = saveObj,
                               outDir = mergeOutDir),
        error = function(err) {
          logThis(paste("While checking cluster uniformity", err),
                  logLevel = 0L)
          logThis("marking pair as not mergable", logLevel = 1L)
          return(checker)
        })

      gc()

      if (mergedClName %in% names(allCheckResults)) {
        # can happen if previous check information is insufficient to update the
        # result
        allCheckResults[[mergedClName]] <- checkResults
      } else {
        allCheckResults <- append(allCheckResults, checkResults)
        names(allCheckResults)[length(allCheckResults)] <- mergedClName
      }

      logThis(paste("Clusters", cl1, "and", cl2,
                    ifelse(checkResults@isUniform, "can", "cannot"),
                    "be merged"), logLevel = 1L)
    }

    return(allCheckResults)
  }

  # sorts all results by significance then merge all possible pairings in order
  mergeAllClusters <- function(outputClusters, allCheckResults) {
    # filters out missing clusters
    allClNames <- unique(outputClusters)
    assert_that(!any(vapply(allClNames, isEmptyName, logical(1L))),
                is_named2(allCheckResults))

    posToKeep <-
      vapply(seq_along(allCheckResults),
             function(r) {
               c(clName1, clName2) %<-%
                 fromMergedName(names(allCheckResults)[r],
                                allClNames)
               return(clName1 %in% allClNames && clName2 %in% allClNames)
             }, logical(1L))
    checkRes <- allCheckResults[posToKeep]

    # filter out non uniform pairs
    # it is assumed that the checkRes have been already
    # updated to align to current checker
    posToKeep <- vapply(seq_along(checkRes),
                        function(r) return(checkRes[[r]]@isUniform),
                        logical(1L))
    checkRes <- checkRes[posToKeep]

    if (is_empty(checkRes)) {
      logThis(paste("No clusters will be merged"), logLevel = 2L)

      return(outputClusters)
    }

    shifts <- vapply(checkRes, calculateThresholdShiftToUniformity, double(1L))
    checkRes <- checkRes[order(shifts)]

    #operate the merges
    for (mergedName in names(checkRes)) {
      allClNames <- unique(outputClusters)
      c(clName1, clName2) %<-%
        fromMergedName(mergedName, allClNames)
      if (!(clName1 %in% allClNames && clName2 %in% allClNames)) {
        logThis(paste0("One or both of the clusters ", clName1, ", ", clName2,
                       " is no more in the clusterization"), logLevel = 3L)
        next
      }
      logThis(paste("Clusters", clName1, "and", clName2, "will be merged"),
              logLevel = 2L)
      outputClusters <-
        mergeClusters(outputClusters,
                      names = c(clName1, clName2),
                      mergedName = toMergedName(clName1, clName2))
      outputClusters <- factorToVector(outputClusters)
    }
    return(outputClusters)
  }

  # start analysis

  logThis("Merging cells' uniform clustering: START", logLevel = 2L)

  assert_that(estimatorsAreReady(objCOTAN),
              msg = paste("Estimators `lambda`, `nu`, `dispersion` are not",
                          "ready: Use proceedToCoex() to prepare them"))

  if (is_empty(checkers)) {
    GDIThreshold <- ifelse(is.finite(GDIThreshold), GDIThreshold, 1.43)

    # keep old behaviour
    thresholdGap <- max(GDIThreshold - 1.37, 0.0)
    numThresholds <- ceiling(thresholdGap / 0.03)
    allThresholds <- 1.37 +
      c(0L, seq_len(numThresholds)) * thresholdGap / numThresholds

    checkers <-
      lapply(allThresholds,
             function(threshold) {
               new("SimpleGDIUniformityCheck",
                   check = new("GDICheck",
                               GDIThreshold = threshold,
                               maxRatioBeyond = 0.01))
             })

    rm(thresholdGap, numThresholds, allThresholds)
  } else {
    assert_that(!is.finite(GDIThreshold),
                msg = paste("Either `checker` object[s] or",
                            "a legacy `GDIThreshold` must be given"))
    if (!is(checkers, "list")) {
      checkers <- list(checkers)
    }
    assert_that(all(vapply(checkers, is, logical(1L), "BaseUniformityCheck")))
  }
  rm(GDIThreshold)
  logThis(paste("The merge algorithm will use", length(checkers), "passes"),
          logLevel = 1L)

  assert_that(isa(allCheckResults, "data.frame"),
              (is_empty(allCheckResults) ||
                 identical(colnames(allCheckResults),
                           colnames(checkersToDF(checkers[[1L]])))),
              msg = "Previous results passed in are of wrong type or columns")

  allCheckResults <- dfToCheckers(allCheckResults, class(checkers[[1L]]))

  if (is_empty(clusters)) {
    # pick the last clusterization
    clusters <- getClusters(objCOTAN)
  }
  outputClusters <-
    factorToVector(asClusterization(clusters, getCells(objCOTAN)))

  if (batchSize == 0L) {
    # default is twice the (2/3) power of the number of clusters
    numCl <- length(unique(outputClusters))
    batchSize <- as.integer(ceiling(1.2 * numCl^(2.0 / 3.0)))
    rm(numCl)
  }

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  outDirCond <- file.path(outDir, cond)
  if (isTRUE(saveObj) && !dir.exists(outDirCond)) {
    dir.create(outDirCond)
  }

  mergeOutDir <- file.path(outDirCond, "leafs_merge")
  if (isTRUE(saveObj) && !dir.exists(mergeOutDir)) {
    dir.create(mergeOutDir)
  }

  iter <- initialIteration - 1L

  for (checker in checkers) {
    logThis(paste0("Start merging nearest clusters - the main threshold is: ",
                   getCheckerThreshold(checker)),
            logLevel = 2L)
    repeat {
      iter <- iter + 1L
      logThis(paste0("Start merging nearest clusters: iteration ", iter),
              logLevel = 3L)

      firstBatch <- is_empty(allCheckResults)
      oldNumClusters <- length(unique(outputClusters))

      clDist <- distancesBetweenClusters(objCOTAN, clusters = outputClusters,
                                         useDEA = useDEA, distance = distance)
      gc()

      if (isTRUE(saveObj)) tryCatch({
          pdf(file.path(mergeOutDir,
                        paste0("dend_iter_", iter, "_tau_",
                               getCheckerThreshold(checker), "_plot.pdf")))

          hcNorm <- hclust(clDist, method = hclustMethod)
          plot(as.dendrogram(hcNorm))
        }, error = function(err) {
          logThis(paste("While saving dendogram plot", err), logLevel = 0L)
        }, finally = {
          # Check for active device
          if (dev.cur() > 1L) {
            dev.off()
          }
        })

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

      pList <- selectPairsList(pList, batchSize, checker, allCheckResults)

      allCheckResults <-
        testPairListMerge(pList, outputClusters, checker, allCheckResults)

      if (!firstBatch) {
        outputClusters <- mergeAllClusters(outputClusters, allCheckResults)
      }

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
                               paste0("all_check_results_", iter, ".csv"))
          write.csv(checkersToDF(allCheckResults), file = outFile, na = "NaN")
        },
        error = function(err) {
          logThis(paste("While saving current clusterization", err),
                  logLevel = 0L)
        })
      if (firstBatch) {
        logThis("Finished the first batch - no merges were executed",
                logLevel = 3L)
      } else if (newNumClusters == oldNumClusters) {
        logThis(paste("None of the remaining tested",
                      "cluster pairs could be merged"), logLevel = 3L)

        # No merges happened -> too low probability of new merges...
        break
      } else {
        logThis(paste("Executed", (oldNumClusters - newNumClusters), "merges"),
                logLevel = 3L)
      }
    }
    logThis(paste0("Executed all merges for threshold ",
                   getCheckerThreshold(checker), " out of ",
                   length(allCheckResults), " checks"), logLevel = 3L)
  }

  logThis(paste0("The final merged clusterization contains [",
                 length(unique(outputClusters)), "] different clusters: ",
                 toString(sort(unique(outputClusters)))), logLevel = 1L)

  # replace the clusters' tags with completely new ones
  if (TRUE) {
    clTags <- sort(unique(outputClusters))

    clTagsMap <- paste0(seq_along(clTags))
    clTagsMap <- factorToVector(niceFactorLevels(clTagsMap))
    clTagsMap <- set_names(clTagsMap, clTags)

    # keep `-1` tag if it has not been merged
    if ("-1" %in% clTags) {
      clTagsMap[["-1"]] <- "-1"
    }

    outputClusters <- set_names(clTagsMap[outputClusters], getCells(objCOTAN))

    checksTokeep <- names(allCheckResults) %in% clTags
    allCheckResults <- allCheckResults[checksTokeep]
    names(allCheckResults) <- clTagsMap[names(allCheckResults)]
  }

  outputCoexDF <-
    tryCatch(DEAOnClusters(objCOTAN, clusters = outputClusters),
             error = function(err) {
               logThis(paste("Calling DEAOnClusters", err), logLevel = 0L)
               return(NULL)
             })

  c(outputClusters, outputCoexDF, permMap) %<-% tryCatch(
    reorderClusterization(objCOTAN, clusters = outputClusters,
                          coexDF = outputCoexDF, reverse = FALSE,
                          keepMinusOne = TRUE, useDEA = useDEA,
                          distance = distance, hclustMethod = hclustMethod),
    error = function(err) {
      logThis(paste("Calling reorderClusterization", err), logLevel = 0L)
      return(list(outputClusters, outputCoexDF))
    })
  names(allCheckResults) <- permMap[names(allCheckResults)]

  if (isTRUE(saveObj)) tryCatch({
    outFile <- file.path(outDirCond, "merge_check_results.csv")
    write.csv(checkersToDF(allCheckResults), file = outFile, na = "NaN")
  })

  logThis("Merging cells' uniform clustering: DONE", logLevel = 2L)

  return(list("clusters" = outputClusters, "coex" = outputCoexDF))
}
