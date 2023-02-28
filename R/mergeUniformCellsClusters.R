
#' mergeUniformCellsClusters
#'
#' @description This function takes in input a `COTAN` object along with a
#'   uniform clusterization and, through the *cosine distance* and the `hclust`,
#'   checks if merging two leaf clusters will form a still uniform cluster
#'   (this is done iteratively). All structures are saved on disk on request.
#'
#' @details This function uses [checkClusterUniformity()] to check whether
#'   merged clusters' uniformity
#'
#' @param objCOTAN a `COTAN` object
#' @param clusters The clusterization to merge. If not given the last available
#'   clusterization will be used, as it is probably the most significant!
#' @param cores number cores used
#' @param distance type of distance to use (default is `cosine`, `euclidean` is
#'   also available)
#' @param hclustMethod default is "ward.D2" but can be any method defined by
#'   [stats::hclust()] function
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output.
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
#' @importFrom dendextend get_nodes_attr
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
#'
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#'
#' mergedList <- mergeUniformCellsClusters(objCOTAN,
#'                                         clusters = clusters,
#'                                         cores = 12,
#'                                         distance = "cosine",
#'                                         hclustMethod = "ward.D2",
#'                                         saveObj = FALSE,
#'                                         outDir = tempdir())
#' mergedClusters <- mergedList[["clusters"]]
#'
#' @rdname mergeUniformCellsClusters
#'

# @param GEO GEO or other data set code
# @param sc.method scRNAseq method
# @param cond sample condition name
# @param markers a `list` of marker genes. Default `NULL`

mergeUniformCellsClusters <- function(objCOTAN,
                                      clusters = NULL,
                                      cores = 1,
                                      distance = "cosine",
                                      hclustMethod = "ward.D2",
                                      saveObj = FALSE,
                                      outDir = ".") {
  logThis("Merging cells' uniform clustering: START", logLevel = 2)

  outputClusters <- clusters
  if (is_empty(outputClusters)) {
    # pick last clusterization
    outputClusters <- getClusterizationData(objCOTAN)[["clusters"]]
  }

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  outDirCond <- file.path(outDir, cond)
  if(!file.exists(outDirCond)){
    dir.create(outDirCond)
  }

  mergeOutDir <- file.path(outDirCond, "leafs_merge")
  if (isTRUE(saveObj) && !file.exists(mergeOutDir)) {
    dir.create(mergeOutDir)
  }

  notMergeable <- c()
  iter <- 0
  repeat {
    iter <- iter + 1
    logThis(paste0("Start merging smallest clusters: iteration ", iter), logLevel = 3)

    oldNumClusters = length(unique(outputClusters))

    c(coexDF, pValueDF) %<-% DEAOnClusters(objCOTAN, clusters = outputClusters)
    gc()

    ## To drop the cells with only zeros
    ## TODO: To be fixed! Where it come from?
    coexDF <- coexDF[, colSums(coexDF != 0) > 0]

    #merge small cluster based on distances
    if (distance == "cosine") {
      coexDist <- cosineDissimilarity(as.matrix(coexDF))
    } else if(distance == "euclidean") {
      coexDist <- dist(t(as.matrix(coexDF)))
    } else {
      stop("only 'cosine' and 'euclidean' distances are supported")
    }

    hcNorm <- hclust(coexDist, method = hclustMethod)
    #plot(hcNorm)

    dend <- as.dendrogram(hcNorm)

    # This checks if any little two pair of leaf clusters of the dendogram
    # could be merged

    id <- c()
    {
      members <- get_nodes_attr(dend, "members")
      for (i in seq_along(members)) {
        if (members[i] == 2) {
          id <- c(id, i + 1, i + 2)
        }
      }
    }

    logThis(paste0("Created leafs ID for merging: ",
                   paste(get_nodes_attr(dend, "label", id = id),
                         collapse = " ")), logLevel = 2)

    p = 1
    while (p < length(id)) {
      logThis("*", logLevel = 1, appendLF = FALSE)

      cl1 <- get_nodes_attr(dend, "label", id = id[p+0])
      cl2 <- get_nodes_attr(dend, "label", id = id[p+1])

      mergedCluster <- paste0(min(cl1, cl2), "_", max(cl1, cl2), "-merge")

      logThis(mergedCluster, logLevel = 3)

      p  <- p+2

      if (mergedCluster %in% notMergeable) {
        logThis(paste0("Clusters ", cl1, " and ", cl2,
                       " already analyzed and not mergeable: skip."),
                logLevel = 3)
        next
      }

      mergedCells <- names(outputClusters)[outputClusters %in% c(cl1, cl2)]

      clusterIsUniform <- checkClusterUniformity(objCOTAN,
                                                 cluster = mergedCluster,
                                                 cells = mergedCells,
                                                 cores = cores,
                                                 saveObj = saveObj,
                                                 outDir = mergeOutDir)

      gc()

      if (!clusterIsUniform) {
        logThis(paste("Merging clusters", cl1, "and", cl2,
                      "results in a too high GDI"), logLevel = 3)

        notMergeable <- c(notMergeable, mergedCluster)
      } else {
        logThis(paste("Clusters", cl1, "and", cl2, "can be merged"), logLevel = 3)

        outputClusters[mergedCells] <- mergedCluster
      }
    }

    if (length(unique(outputClusters)) == oldNumClusters) {
      # no merges happened: stop!
      break
    }
  }

  logThis(paste("The final merged clusterization contains [",
                length(unique(outputClusters)), "] different clusters:",
                paste0(unique(sort(outputClusters)), collapse = ", ")), logLevel = 3)

  # replace the clusters' tags
  # non merged clusters keep their tag, while merged one use the lesser one!
  # there might be holes in the enumeration!
  {
    clTags <- levels(factor(outputClusters))
    # here we use the fact that the merged clusters have the smaller cluster as first
    clTagsMap <- str_split(clTags, pattern = "_", simplify = TRUE)[,1]
    clTagsMap <- set_names(clTagsMap, clTags)

    outputClusters <- clTagsMap[outputClusters]
    outputClusters <- set_names(outputClusters, getCells(objCOTAN))

    colnames(coexDF)   <- clTagsMap[colnames(coexDF)]
    colnames(pValueDF) <- clTagsMap[colnames(pValueDF)]
  }

  logThis("Merging cells' uniform clustering: DONE", logLevel = 2)

  return(list("clusters" = outputClusters,
              "coexDF" = coexDF, "pValueDF" = pValueDF))
}
