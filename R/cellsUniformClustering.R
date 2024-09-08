
#' @title Get a clusterization running the `Seurat` package
#'
#' @description The function uses the [Seurat-package] to clusterize the given
#'   counts raw data.
#'
#' @details The parameter resolution is set at 0.5 initially, but in case of too
#'   few clusters it can be raised up to 2.5.
#'
#' @param rawData The raw counts
#' @param cond The sample condition
#' @param iter The current iteration
#' @param initialResolution The resolution to use at first in the
#'   *clusterization* algorithm
#' @param minNumClusters The minimum number of *clusters* expected from this
#'   *clusterization*. In cases it is not reached, it will increase the
#'   resolution of the *clusterization*.
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDirCond an existing directory for the analysis output.
#'
#' @returns a list with a `Seurat` object along a Boolean on whether maximum
#'   resolution has been used
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat DimPlot
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @importFrom ggplot2 annotate
#'
#' @importFrom withr with_options
#'
#' @noRd
#'
seuratClustering <- function(rawData, cond, iter, initialResolution,
                             minNumClusters, saveObj, outDirCond) {
  tryCatch({
    logThis("Creating Seurat object: START", logLevel = 2L)

    srat <- CreateSeuratObject(counts = rawData,
                               project = paste0(cond, "_reclustering_", iter),
                               min.cells = if (iter == 1L) 3L else 1L,
                               min.features = if (iter == 1L) 4L else 2L)
    srat <- NormalizeData(srat)
    srat <- FindVariableFeatures(srat, selection.method = "vst",
                                 nfeatures = 2000L)
    srat <- ScaleData(srat, features = rownames(srat))

    maxRows <- nrow(srat@meta.data) - 1L
    srat <- RunPCA(srat, features = VariableFeatures(object = srat),
                   npcs = min(50L, maxRows))

    srat <- FindNeighbors(srat, dims = 1L:min(25L, maxRows))

    resolution <- initialResolution
    resolutionStep <- 0.5
    maxResolution <- initialResolution + 10.0 * resolutionStep
    usedMaxResolution <- FALSE
    repeat {
      srat <- FindClusters(srat, resolution = resolution, algorithm = 2L)

      # The next lines are necessary to make cluster smaller while
      # the number of residual cells decrease and to stop clustering
      # if the algorithm gives too many singletons.
      usedMaxResolution <- (resolution + 0.1 * resolutionStep) > maxResolution
      numClusters <- nlevels(factor(srat[["seurat_clusters", drop = TRUE]]))
      if (numClusters > minNumClusters || usedMaxResolution) {
        break
      }

      logThis(paste("Number of clusters is too small.",
                    "Reclustering at resolution higher than:", resolution),
              logLevel = 3L)

      resolution <- resolution + resolutionStep
    }

    logThis(paste("Used resolution for Seurat clusterization is:", resolution),
            logLevel = 2L)

    # disable annoying warning about Seurat::RunUMAP()
    srat <- with_options(list(Seurat.warn.umap.uwot = FALSE),
                         RunUMAP(srat, umap.method = "uwot", metric = "cosine",
                                 dims = 1L:min(c(50L, maxRows))))

    if (isTRUE(saveObj)) tryCatch({
        outFile <- file.path(outDirCond, paste0("pdf_umap_", iter, ".pdf"))
        logThis(paste("Creating PDF UMAP in file: ", outFile), logLevel = 2L)
        pdf(outFile)

        if (iter == 1L) {
          plot(DimPlot(srat, reduction = "umap", label = FALSE,
                       group.by = "orig.ident"))
        }

        plot(DimPlot(srat, reduction = "umap", label = TRUE) +
             annotate(geom = "text", x = 0.0, y = 30.0, color = "black",
                      label = paste0("Cells number: ", ncol(rawData), "\n",
                                     "Cl. resolution: ", resolution)))
      }, error = function(err) {
        logThis(paste("While saving seurat UMAP plot", err), logLevel = 1L)
      }, finally = {
        # Check for active device
        if (dev.cur() > 1L) {
          dev.off()
        }
      })
    gc()

    logThis("Creating Seurat object: DONE", logLevel = 2L)

    # returned objects
    return(list("SeuratObj" = srat, "UsedMaxResolution" = usedMaxResolution))
  },
  error = function(e) {
    logThis(msg = paste("Seurat clusterization failed with", ncol(rawData),
                        "cells with the following error:"), logLevel = 1L)
    logThis(msg = conditionMessage(e), logLevel = 1L)
    return(list("SeuratObj" = NULL, "UsedMaxResolution" = FALSE))
  })
}

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
#' @param cores number of cores to use. Default is 1.
#' @param maxIterations max number of re-clustering iterations. It defaults to
#'   \eqn{25}
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
#' @param initialResolution a number indicating how refined are the clusters
#'   before checking for **uniformity**. It defaults to \eqn{0.8}, the same as
#'   [Seurat::FindClusters()]
#' @param initialClusters an existing *clusterization* to use as starting point:
#'   the *clusters* deemed **uniform** will be kept and the rest processed as
#'   normal
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
#' @returns `cellsUniformClustering()` returns a `list` with 2 elements:
#'   * `"clusters"` the newly found cluster labels array
#'   * `"coex"` the associated `COEX` `data.frame`
#'
#' @export
#'
#' @importFrom rlang set_names
#' @importFrom rlang is_null
#'
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
#' @rdname UniformClusters
#'

cellsUniformClustering <- function(objCOTAN,
                                   checker = NULL,
                                   GDIThreshold = NaN,
                                   cores = 1L,
                                   maxIterations = 25L,
                                   optimizeForSpeed = TRUE,
                                   deviceStr = "cuda",
                                   initialClusters = NULL,
                                   initialResolution = 0.8,
                                   useDEA = TRUE,
                                   distance = NULL,
                                   hclustMethod = "ward.D2",
                                   saveObj = TRUE, outDir = ".") {
  logThis("Creating cells' uniform clustering: START", logLevel = 2L)

  assert_that(estimatorsAreReady(objCOTAN),
              msg = paste("Estimators lambda, nu, dispersion are not ready:",
                          "Use proceeedToCoex() to prepare them"))

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  outDirCond <- file.path(outDir, cond)
  if (isTRUE(saveObj) && !dir.exists(outDirCond)) {
    dir.create(outDirCond)
  }

  splitOutDir <- file.path(outDirCond, "reclustering")
  if (isTRUE(saveObj) && !dir.exists(splitOutDir)) {
    dir.create(splitOutDir)
  }

  saveSeuratObj <- saveObj && FALSE

  outputClusters <- set_names(rep(NA, length = getNumCells(objCOTAN)),
                              getCells(objCOTAN))

  iter <- 0L
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
    assert_that(is.finite(GDIThreshold),
                msg = paste("Either a `checker` object or",
                            "a legacy `GDIThreshold` must be given"))
  }

  repeat {
    iter <- iter + 1L

    logThis(paste0("In iteration ", iter, " "), logLevel = 1L, appendLF = FALSE)
    logThis(paste("the number of cells to re-cluster is",
                  sum(is.na(outputClusters)), "cells belonging to",
                  numClustersToRecluster, "clusters"), logLevel = 2L)

    #Step 1
    minNumClusters <- floor(1.2 * numClustersToRecluster) + 1L
    c(objSeurat, usedMaxResolution) %<-%
      seuratClustering(rawData = getRawData(objCOTAN)[, is.na(outputClusters)],
                       cond = cond, iter = iter,
                       initialResolution = initialResolution,
                       minNumClusters = minNumClusters,
                       saveObj = saveObj, outDirCond = splitOutDir)

    if (is_null(objSeurat)) {
      logThis(paste("NO new possible uniform clusters!",
                    "Unclustered cell left:", sum(is.na(outputClusters))),
              logLevel = 1L)
      break
    }

    if (saveSeuratObj && iter == 1L) tryCatch({
      # save the Seurat object to file to be reloaded later
      saveRDS(objSeurat,
              file.path(outDirCond, "Seurat_obj_with_cotan_clusters.RDS"))
      },
      error = function(err) {
        logThis(paste("While saving seurat object", err), logLevel = 1L)
      })

    metaData <- objSeurat@meta.data

    rm(objSeurat)
    gc()

    # Step 2

    # The next lines check for each new cell cluster if it is homogeneous.
    # In cases they are not, save the corresponding cells in cellsToRecluster

    numClustersToRecluster <- 0L
    cellsToRecluster <- vector(mode = "character")

    testClusters <- factor(metaData[["seurat_clusters"]])
    allCells <- rownames(metaData)
    names(testClusters) <- allCells
    if (iter == 1L && !is_null(initialClusters)) {
      assert_that(setequal(names(initialClusters), getCells(objCOTAN)),
                  msg = "Given clusterization has the wrong set of cells")
      logThis("Using passed in clusterization", logLevel = 3L)
      testClusters <- factor(initialClusters)
      allCells <- names(initialClusters)
    }
    testClList <- toClustersList(testClusters)

    globalClName <- ""

    for (clName in names(testClList)) {
      logThis("*", logLevel = 1L, appendLF = FALSE)
      logThis(paste0(" checking uniformity of cluster '", clName,
                     "' of ", length(testClList), " clusters"),
              logLevel = 2L)

      if (clName == "singleton") {
        next
      }

      globalClName <-
        paste0(str_pad(iter, width = 2L, pad = "0"), "_",
               str_pad(clName, width = 4L, pad = "0"))

      cells <- testClList[[clName]]
      if (length(cells) < 20L) {
        logThis(paste("cluster", globalClName, "has too few cells:",
                      "will be reclustered!"), logLevel = 1L)

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
            logThis("marking cluster as not uniform", logLevel = 1L)
            return(checker)
          })

        invisible(validObject(checkResults))

        allCheckResults <- append(allCheckResults, checkResults)
        names(allCheckResults)[length(allCheckResults)] <- globalClName

        if (!checkResults@isUniform) {
          logThis(paste("cluster", globalClName, "has too high GDI:",
                        "will be reclustered!"), logLevel = 1L)

          numClustersToRecluster <- numClustersToRecluster + 1L
          cellsToRecluster <- c(cellsToRecluster, cells)
        } else {
          logThis(paste("cluster", globalClName, "is uniform"), logLevel = 1L)
        }

        rm(checkResults)
        gc()
      }
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
      if (isTRUE(usedMaxResolution)) {
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

    if (length(cellsToRecluster) < 40L || iter >= maxIterations) {
      logThis(paste("NO new possible uniform clusters!",
                    "Unclustered cell left:", length(cellsToRecluster)),
              logLevel = 1L)
      break
    }
    rm(cellsToRecluster)
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
      checker@clusterSize <- length(unclusteredCells)
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

      if (saveSeuratObj) {
        clusterizationName <-
          paste0(as.roman(length(getClusterizations(objCOTAN)) + 1L))

        srat <-
          readRDS(file.path(outDirCond, "Seurat_obj_with_cotan_clusters.RDS"))

        if (!setequal(rownames(srat@meta.data), names(outputClusters))) {
          warning("List of cells got corrupted")
          return(outputList)
        }

        srat@meta.data <-
          setColumnInDF(srat@meta.data,
                        colToSet = outputClusters[rownames(srat@meta.data)],
                        colName = paste0("COTAN_", clusterizationName))

        if (!dim(srat)[[2L]] == getNumCells(objCOTAN)) {
          warning("Number of cells got wrong")
          return(outputList)
        }

        logThis("Cluster, UMAP and Saving the Seurat dataset", logLevel = 2L)

        srat <- FindNeighbors(srat, dims = 1L:25L)
        srat <- FindClusters(srat, resolution = 0.5, algorithm = 2L)
        srat <- RunUMAP(srat, umap.method = "uwot",
                        metric = "cosine", dims = 1L:25L)

        saveRDS(srat, file.path(outDirCond,
                                "Seurat_obj_with_cotan_clusters.RDS"))
      }
    },
    error = function(err) {
      logThis(paste("While saving seurat object", err), logLevel = 1L)
    }
  )

  rm(srat)
  gc()

  logThis("Creating cells' uniform clustering: DONE", logLevel = 2L)

  return(outputList)
}
