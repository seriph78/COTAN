# --------------------- Uniform Clusters ----------------------

#' @title Uniform Clusters
#'
#' @description This group of functions takes in input a `COTAN` object and
#'   handle the task of dividing the dataset into **Uniform Clusters**, that is
#'   *clusters* that have an homogeneous genes' expression. This condition is
#'   checked by calculating the `GDI` of the *cluster* and verifying that no
#'   more than a small fraction of the genes have their `GDI` level above the
#'   given `GDIThreshold`
#'
#' @name UniformClusters
NULL


### ------ Seurat Clustering -------

#' @title Get a clusterization running the `Seurat` package
#'
#' @description The function uses the \pkg{Seurat} to clusterize the given
#'   counts raw data.
#'
#' @details The parameter resolution is set at 0.5 initially, but in case of too
#'   few clusters it can be raised up to 2.5.
#'
#' @param rawData The raw counts
#' @param initialResolution The resolution to use at first in the
#'   *clusterization* algorithm
#' @param minNumClusters The minimum number of *clusters* expected from this
#'   *clusterization*. In cases it is not reached, it will increase the
#'   resolution of the *clusterization*
#' @param useCoexEigen Boolean to determine whether to project the data `matrix`
#'   onto the first eigenvectors of the **COEX** `matrix` or instead restrict
#'   the data `matrix` to the selected genes before applying the `PCA` reduction
#' @param dataMethod selects the method to use to create the `data.frame` to
#'   pass to the [UMAPPlot()]. See [getDataMatrix()] for more details.
#' @param genesSel Decides whether and how to perform the gene-selection
#'   (defaults to `"HVG_Seurat"`). See [genesSelector()] for more details.
#' @param numGenes the number of genes to select using the above method. Will be
#'   ignored when an explicit list of genes has been passed in
#' @param numReducedComp the number of calculated **RDM** components
#'
#' @returns a list with:
#'   * `"SeuratClusters"` a `Seurat` *clusterization*
#'   * `"CellsRDM"` the Reduced Data Matrix
#'   * `"Resolution"` the used resolution
#'   * `"UsedMaxResolution"` whether the maximum resolution has been reached
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat CreateDimReducObject
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat Idents
#'
#' @importFrom assertthat assert_that
#'
#' @noRd
#'

seuratClustering <- function(objCOTAN,
                             initialResolution, minNumClusters,
                             useCoexEigen, dataMethod,
                             genesSel, numGenes,
                             numReducedComp) {
  tryCatch({
    logThis("Creating new clusterization: START", logLevel = 2L)

    assert_that(numReducedComp <= getNumGenes(objCOTAN))

    numReducedCompToCalc <- numReducedComp + 15L
    cellsRDM <- calculateReducedDataMatrix(
      objCOTAN, useCoexEigen = FALSE,
      dataMethod = "LogNorm", numComp = numReducedCompToCalc,
      genesSel = genesSel, numGenes = numGenes)

    assert_that(identical(dim(cellsRDM),
                          c(getNumCells(objCOTAN), numReducedCompToCalc)),
                msg = "Returned PCA matrix has wrong dimensions")

    # Create the Seurat object
    srat <- CreateSeuratObject(counts = getRawData(objCOTAN),
                               project = "clusterization")

    # Add RDM manually
    srat[["pca"]] <-
      CreateDimReducObject(embeddings = cellsRDM,
                           key = "PC_",
                           assay = "RNA")

    srat <- FindNeighbors(srat, dims = seq_len(numReducedComp))

    resolution <- initialResolution
    resolutionStep <- 0.5
    maxResolution <- initialResolution + 10.0 * resolutionStep

    seuratClusters <- NULL
    usedMaxResolution <- FALSE
    repeat {
      seuratClusters <- Idents(FindClusters(
        srat,
        resolution = resolution,
        algorithm = 2L,      # Louvain (refined)
        random.seed = 137    # controls igraph::cluster_louvain()
      ))

      # The next lines are necessary to make cluster smaller while
      # the number of residual cells decrease and to stop clustering
      # if the algorithm has gone for too long
      usedMaxResolution <- (resolution + 0.1 * resolutionStep) > maxResolution
      if (nlevels(seuratClusters) > minNumClusters || usedMaxResolution) {
        break
      }

      logThis(paste("Number of clusters is too small.",
                    "Reclustering at resolution higher than:", resolution),
              logLevel = 3L)

      resolution <- resolution + resolutionStep
    }

    logThis(paste("Used resolution for Seurat clusterization is:", resolution),
            logLevel = 2L)

    logThis("Creating new clusterization: DONE", logLevel = 2L)

    rm(srat)
    gc()

    # returned objects
    return(list("SeuratClusters" = seuratClusters, "CellsRDM" = cellsRDM,
                "Resolution" = resolution,
                "UsedMaxResolution" = usedMaxResolution))
  },
  error = function(e) {
    logThis(msg = paste("Seurat clusterization failed with",
                        getNumCells(objCOTAN),
                        "cells with the following error:"), logLevel = 1L)
    logThis(msg = conditionMessage(e), logLevel = 1L)
    return(list("SeuratClusters" = NULL, "CellsRDM" = NULL,
                "Resolution" = NaN, "UsedMaxResolution" = FALSE))
  })
}

## -------- Cells Uniform Clustering --------

#' @details `cellsUniformClustering()` finds a **Uniform** *clusterizations* by
#'   means of the `GDI`. Once a preliminary *clusterization* is obtained from
#'   the `Seurat-package` methods, each *cluster* is checked for **uniformity**
#'   via the function [checkClusterUniformity()]. Once all *clusters* are
#'   checked, all cells from the **non-uniform** clusters are pooled together
#'   for another iteration of the entire process, until all *clusters* are
#'   deemed **uniform**. In the case only a few cells are left out (\eqn{\leq
#'   50}), those are flagged as `"-1"` and the process is stopped.
#'
#' @param objCOTAN a `COTAN` object
#' @param checker the object that defines the method and the threshold to
#'   discriminate whether a *cluster* is *uniform transcript*. See
#'   [UniformTranscriptCheckers] for more details
#' @param GDIThreshold legacy. The threshold level that is used in a
#'   [SimpleGDIUniformityCheck-class]. It defaults to \eqn{1.40}
#' @param initialResolution a number indicating how refined are the clusters
#'   before checking for **uniformity**. It defaults to \eqn{0.8}, the same as
#'   [Seurat::FindClusters()]
#' @param maxIterations max number of re-clustering iterations. It defaults to
#'   \eqn{25}
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
#'   faster but less precise
#' @param distance type of distance to use. Default is `"cosine"` for *DEA* and
#'   `"euclidean"` for *Zero-One*. Can be chosen among those supported by
#'   [parallelDist::parDist()]
#' @param useCoexEigen Boolean to determine whether to project the data `matrix`
#'   onto the first eigenvectors of the **COEX** `matrix` or instead restrict
#'   the data `matrix` to the selected genes before applying the `PCA` reduction
#' @param dataMethod selects the method to use to create the `data.frame` to
#'   pass to the [UMAPPlot()]. See [getDataMatrix()] for more details.
#' @param genesSel Decides whether and how to perform the gene-selection
#'   (defaults to `"HVG_Seurat"`). See [genesSelector()] for more details.
#' @param numGenes the number of genes to select using the above method. Will be
#'   ignored when an explicit list of genes has been passed in
#' @param numReducedComp the number of calculated **RDM** components
#' @param hclustMethod It defaults is `"ward.D2"` but can be any of the methods
#'   defined by the [stats::hclust()] function.
#' @param initialClusters an existing *clusterization* to use as starting point:
#'   the *clusters* deemed **uniform** will be kept and the remaining cells will
#'   be processed as normal
#' @param initialIteration the number associated tot he first iteration; it
#'   defaults to 1. Useful in case of restart of the procedure to avoid
#'   intermediate data override
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output. The effective
#'   output will be paced in a sub-folder.
#'
#' @returns `cellsUniformClustering()` returns a `list` with 2 elements:
#'   * `"clusters"` the newly found cluster labels array
#'   * `"coex"` the associated `COEX` `data.frame`
#'
#' @export
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_null
#'
#' @importFrom stringr str_equal
#' @importFrom stringr str_detect
#' @importFrom stringr str_pad
#'
#' @importFrom utils as.roman
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom methods new
#' @importFrom methods validObject
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.cur
#'
#' @importFrom ggplot2 annotate
#'
#' @importFrom withr with_options
#'
#' @rdname UniformClusters
#'

cellsUniformClustering <- function(objCOTAN,
                                   checker = NULL,
                                   GDIThreshold = NaN,
                                   initialResolution = 0.8,
                                   maxIterations = 25L,
                                   cores = 1L,
                                   optimizeForSpeed = TRUE,
                                   deviceStr = "cuda",
                                   useDEA = TRUE,
                                   distance = NULL,
                                   useCoexEigen = FALSE,
                                   dataMethod = "",
                                   genesSel = "HVG_Seurat",
                                   numGenes = 2000L,
                                   numReducedComp = 25L,
                                   hclustMethod = "ward.D2",
                                   initialClusters = NULL,
                                   initialIteration = 1L,
                                   saveObj = TRUE,
                                   outDir = ".") {
  logThis("Creating cells' uniform clustering: START", logLevel = 2L)

  assert_that(estimatorsAreReady(objCOTAN),
              msg = paste("Estimators `lambda`, `nu`, `dispersion` are not",
                          "ready: Use proceedToCoex() to prepare them"))

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  outDirCond <- file.path(outDir, cond)
  if (isTRUE(saveObj) && !dir.exists(outDirCond)) {
    dir.create(outDirCond)
  }

  splitOutDir <- file.path(outDirCond, "reclustering")
  if (isTRUE(saveObj) && !dir.exists(splitOutDir)) {
    dir.create(splitOutDir)
  }

  outputClusters <- set_names(rep(NA, length = getNumCells(objCOTAN)),
                              getCells(objCOTAN))

  if (!is.integer(initialIteration)) {
    initialIteration <- 1L
  }

  iter <- initialIteration - 1L
  iterReset <- -1L
  numClustersToRecluster <- 0L
  srat <- NULL
  allCheckResults <- list()

  if (is.null(checker)) {
    GDIThreshold <- ifelse(is.finite(GDIThreshold), GDIThreshold, 1.40)
    checker <- new("SimpleGDIUniformityCheck",
                   check = new("GDICheck",
                               GDIThreshold = GDIThreshold,
                               maxRatioBeyond = 0.01))
  } else {
    assert_that(!is.finite(GDIThreshold),
                msg = paste("Either a `checker` object or",
                            "a legacy `GDIThreshold` must be given"))
  }

  if (isEmptyName(dataMethod)) {
    dataMethod = "LogNormalized"
  }

  repeat {
    iter <- iter + 1L

    logThis(paste0("In iteration ", iter, " "), logLevel = 1L, appendLF = FALSE)
    logThis(paste("the number of cells to re-cluster is",
                  sum(is.na(outputClusters)), "cells belonging to",
                  numClustersToRecluster, "clusters"), logLevel = 2L)

    # create COTAN sub-object
    cellsToDrop <- getCells(objCOTAN)[!is.na(outputClusters)]
    subObj <- dropGenesCells(objCOTAN, cells = cellsToDrop)

    if (str_equal(genesSel, "HGDI")) {
      subObj <-
        proceedToCoex(subObj, calcCoex = TRUE, cores = cores,
                      optimizeForSpeed = optimizeForSpeed,
                      deviceStr = deviceStr,
                      saveObj = FALSE, outDir = outDirCond)
    }

    #Step 1
    minNumClusters <- floor(1.2 * numClustersToRecluster) + 1L
    c(testClusters, cellsRDM, resolution, usedMaxResolution) %<-%
      seuratClustering(subObj, initialResolution = initialResolution,
                       minNumClusters = minNumClusters,
                       useCoexEigen = useCoexEigen,
                       dataMethod = dataMethod,
                       numReducedComp = numReducedComp,
                       genesSel = genesSel, numGenes = numGenes)

    if (is_null(testClusters)) {
      logThis(paste("NO new possible uniform clusters!",
                    "Unclustered cell left:", sum(is.na(outputClusters))),
              logLevel = 1L)
      break
    }

    ## print umap graph
    if (isTRUE(saveObj)) tryCatch({
      outFile <- file.path(splitOutDir, paste0("pdf_umap_", iter, ".pdf"))
      logThis(paste("Creating PDF UMAP in file: ", outFile), logLevel = 2L)
      pdf(outFile)

      plot(UMAPPlot(dataIn = cellsRDM,
                    clusters = getCondition(subObj),
                    title = paste0("Cells number: ", nrow(cellsRDM))))

      plot(UMAPPlot(dataIn = cellsRDM,
                    clusters = testClusters,
                    title = paste0("Cells number: ", nrow(cellsRDM), "\n",
                                   "Cl. resolution: ", resolution)))
    }, error = function(err) {
      logThis(paste("While saving seurat UMAP plot", err), logLevel = 1L)
    }, finally = {
      # Check for active device
      if (dev.cur() > 1L) {
        dev.off()
      }
    })

    # Step 2

    # The next lines check for each new cell cluster if it is homogeneous.
    # In cases they are not, save the corresponding cells in cellsToRecluster

    numClustersToRecluster <- 0L
    cellsToRecluster <- vector(mode = "character")

    if (iter == initialIteration && !is_null(initialClusters)) {
      logThis("Using passed in clusterization", logLevel = 3L)
      testClusters <- asClusterization(initialClusters, getCells(objCOTAN))
    }
    allCells <- names(testClusters)
    testClList <- toClustersList(testClusters)

    globalClName <- ""

    minimumUTClusterSize <- 20L
    maxClusterSize <- max(lengths(testClList))

    if (maxClusterSize >= minimumUTClusterSize) {
      for (clName in names(testClList)) {
        logThis("*", logLevel = 1L, appendLF = FALSE)
        logThis(paste0(" checking uniformity of cluster '", clName,
                       "' of ", length(testClList), " clusters"),
                logLevel = 2L)

        globalClName <-
          paste0(str_pad(iter, width = 2L, pad = "0"), "_",
                 str_pad(clName, width = 4L, pad = "0"))

        cells <- testClList[[clName]]
        if (length(cells) < 20L) {
          logThis(paste("cluster", globalClName, "has too few cells:",
                        "will be reclustered!"), logLevel = 2L)

          numClustersToRecluster <- numClustersToRecluster + 1L
          cellsToRecluster <- c(cellsToRecluster, cells)
        } else {
          checkResults <- tryCatch(
            checkClusterUniformity(objCOTAN = objCOTAN,
                                   clusterName = globalClName,
                                   cells = cells,
                                   checker = checker,
                                   cores = cores,
                                   optimizeForSpeed = optimizeForSpeed,
                                   deviceStr = deviceStr,
                                   saveObj = saveObj,
                                   outDir = splitOutDir),
            error = function(err) {
              logThis(paste("while checking cluster uniformity", err),
                      logLevel = 0L)
              logThis("marking cluster as not uniform", logLevel = 2L)
              return(checker)
            })

          invisible(validObject(checkResults))

          allCheckResults <- append(allCheckResults, checkResults)
          names(allCheckResults)[length(allCheckResults)] <- globalClName

          if (!checkResults@isUniform) {
            logThis(paste("cluster", globalClName, "has too high GDI:",
                          "will be reclustered!"), logLevel = 2L)

            numClustersToRecluster <- numClustersToRecluster + 1L
            cellsToRecluster <- c(cellsToRecluster, cells)
          } else {
            logThis(paste("cluster", globalClName, "is uniform"), logLevel = 2L)
          }
          logThis("", logLevel = 2L)

          rm(checkResults)
          gc()
        }
      }
    } else {
      # all clusters are too small: nothing to do
      logThis(paste("All clusters in iteration ", iter, "have too few cells"),
              logLevel = 2L)
      numClustersToRecluster <- length(testClList)
      cellsToRecluster <- allCells
    }

    logThis("", logLevel = 1L)
    logThis(paste("Found", length(testClList) - numClustersToRecluster,
                  "uniform and ", numClustersToRecluster,
                  "non-uniform clusters"), logLevel = 2L)

    if (numClustersToRecluster == length(testClList)) {
      warning("In iteration '", iter, "' no uniform clusters found!")
      # Another iteration can be attempted as the minimum number of clusters
      # will be higher. This happens unless the resolution already reached
      # its maximum. In the latter case we simply stop here.
      if (isTRUE(usedMaxResolution) || maxClusterSize < minimumUTClusterSize) {
        logThis("Max resolution reached", logLevel = 1L)
        if (iterReset != -1L) {
          logThis("Cannot clusterize anything more", logLevel = 2L)
          break
        } else {
          # try to recluster remaining cells from scratch
          logThis("Trying to recluster remaining cells from scratch",
                  logLevel = 2L)
          numClustersToRecluster <- 0L
          iterReset <- iter
        }
      }
    } else {
      iterReset <- -1L
    }

    # Step 3: save the already uniform clusters keeping track of the iteration
    if (TRUE) {
      flagInUniformCl <- !allCells %in% cellsToRecluster
      outputClusters[allCells[flagInUniformCl]] <-
        paste0(str_pad(iter, width = 2L, pad = "0"), "_",
               str_pad(testClusters[flagInUniformCl], width = 4L, pad = "0"))
    }

    if (isTRUE(saveObj)) tryCatch({
        outFile <- file.path(splitOutDir,
                             paste0("partial_clusterization_", iter, ".csv"))
        write.csv(outputClusters, file = outFile)

        outFile <- file.path(splitOutDir,
                             paste0("all_check_results_", iter, ".csv"))
        write.csv(checkersToDF(allCheckResults), file = outFile, na = "NaN")
    },
      error = function(err) {
        logThis(paste("While saving current clusterization", err),
                logLevel = 0L)
      }
    )

    if (sum(is.na(outputClusters)) != length(cellsToRecluster)) {
      warning("Some problems in cells reclustering")
      break
    }

    if (length(cellsToRecluster) < 40L
        || iter > maxIterations + initialIteration) {
      logThis(paste("NO new possible uniform clusters!",
                    "Unclustered cell left:", length(cellsToRecluster)),
              logLevel = 1L)
      break
    }

    rm(cellsToRecluster)
    gc()
  } # End repeat

  logThis(paste("The final raw clusterization contains [",
                length(unique(outputClusters)), "] different clusters:",
                toString(unique(sort(outputClusters)))), logLevel = 3L)

  # replace the clusters' tags
  if (TRUE) {
    clTags <- unique(sort(outputClusters))

    clTagsMap <- paste0(seq_along(clTags))
    clTagsMap <- factorToVector(niceFactorLevels(clTagsMap))
    clTagsMap <- set_names(clTagsMap, clTags)

    unclusteredCells <- is.na(outputClusters)
    outputClusters[!unclusteredCells] <-
      clTagsMap[outputClusters[!unclusteredCells]]
    outputClusters[unclusteredCells] <- "-1"
    outputClusters <- set_names(outputClusters, getCells(objCOTAN))

    checksTokeep <- names(allCheckResults) %in% clTags
    allCheckResults <- allCheckResults[checksTokeep]
    names(allCheckResults) <- clTagsMap[names(allCheckResults)]
    if (any(unclusteredCells)) {
      checker@clusterSize <- sum(unclusteredCells)
      allCheckResults <- append(allCheckResults, checker)
      names(allCheckResults)[length(allCheckResults)] <- "-1"
    }
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

  outputList <- list("clusters" = factor(outputClusters), "coex" = outputCoexDF)

  if (isTRUE(saveObj)) tryCatch({
      outFile <- file.path(outDirCond, "split_check_results.csv")
      write.csv(checkersToDF(allCheckResults), file = outFile, na = "NaN")
    },
    error = function(err) {
      logThis(paste("While saving results csv", err), logLevel = 1L)
    }
  )

  logThis("Creating cells' uniform clustering: DONE", logLevel = 2L)

  return(outputList)
}
