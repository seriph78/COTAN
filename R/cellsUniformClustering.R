
#' Get a clusterization running the `Seurat` package
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
#' @param outDir an existing directory for the analysis output.
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
                             minNumClusters, saveObj, outDir) {
  ret <- tryCatch({
    logThis("Creating Seurat object: START", logLevel = 2L)

    srat <- CreateSeuratObject(counts = as.data.frame(rawData),
                               project = paste0(cond, "_reclustering_", iter),
                               min.cells = if (iter == 0L) 3L else 1L,
                               min.features = if (iter == 0L) 4L else 2L)
    srat <- NormalizeData(srat)
    srat <- FindVariableFeatures(srat, selection.method = "vst",
                                 nfeatures = 2000L)
    srat <- ScaleData(srat, features = rownames(srat))

    maxRows <- nrow(srat@meta.data) - 1L
    srat <- RunPCA(srat, features = VariableFeatures(object = srat),
                   npcs = min(50L, maxRows))

    srat <- FindNeighbors(srat, dims = 1L:min(25L, maxRows))

    resolution <- initialResolution
    repeat {
      srat <- FindClusters(srat, resolution = resolution, algorithm = 2L)

      # The next lines are necessary to make cluster smaller while
      # the number of residual cells decrease and to stop clustering
      # if the algorithm gives too many singletons.
      if ((minNumClusters <
             nlevels(factor(srat[["seurat_clusters", drop = TRUE]]))) ||
          (resolution > initialResolution + 5.0)) {
        break
      }

      logThis(paste("Number of clusters is too small.",
                    "Reclustering at higher resolution:", resolution),
              logLevel = 3L)

      resolution <- resolution + 0.5
    }

    logThis(paste("Used resolution for Seurat clusterization is:", resolution),
            logLevel = 2L)

    # disable annoying warning about Seurat::RunUMAP()
    srat <- with_options(list(Seurat.warn.umap.uwot = FALSE),
                         RunUMAP(srat, umap.method = "uwot", metric = "cosine",
                                 dims = 1L:min(c(50L, maxRows))))

    if (isTRUE(saveObj)) {
      logThis(paste0("Creating PDF UMAP in file:",
                     file.path(outDir, "pdf_umap.pdf")), logLevel = 2L)
      pdf(file.path(outDir, "pdf_umap.pdf"))

      if (iter == 0L) {
        plot(DimPlot(srat, reduction = "umap", label = FALSE,
                     group.by = "orig.ident"))
      }

      plot(DimPlot(srat, reduction = "umap", label = TRUE) +
           annotate(geom = "text", x = 0.0, y = 30.0, color = "black",
                    label = paste0("Cells number: ", ncol(rawData), "\n",
                                   "Cl. resolution: ", resolution)))

      dev.off()
      gc()
    }

    logThis("Creating Seurat object: DONE", logLevel = 2L)

    # returned objects
    list("SeuratObj" = srat,
         "UsedMaxResolution" = resolution > initialResolution + 1.5)
  },
  error = function(e) {
    logThis(msg = paste("Seurat clusterization failed with", ncol(rawData),
                        "cells with the following error:"), logLevel = 1L)
    logThis(msg = conditionMessage(e), logLevel = 1L)
    return(list("SeuratObj" = NULL, "UsedMaxResolution" = FALSE))
  })

  return(ret)
}

# --------------------- Uniform Clusters ----------------------

#' Uniform Clusters
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
#' @param GDIThreshold the threshold level that discriminates uniform
#'   *clusters*. It defaults to \eqn{1.4}
#' @param cores number of cores used
#' @param initialResolution a number indicating how refined are the clusters
#'   before checking for **uniformity**. It defaults to \eqn{0.8}, the same as
#'   [Seurat::FindClusters()]
#' @param maxIterations max number of re-clustering iterations. It defaults to
#'   \eqn{25}
#' @param distance type of distance to use (default is `"kullback"`, `"cosine"`
#'   and the others from [parallelDist::parDist()] are also available)
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
#' @importFrom zeallot `%<-%`
#' @importFrom zeallot `%->%`
#'
#' @rdname UniformClusters
#'

cellsUniformClustering <- function(objCOTAN,  GDIThreshold = 1.4,
                                   cores = 1L,
                                   maxIterations = 25L,
                                   initialResolution = 0.8,
                                   distance = "kullback",
                                   hclustMethod = "ward.D2",
                                   saveObj = TRUE, outDir = ".") {
  logThis("Creating cells' uniform clustering: START", logLevel = 2L)

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  outDirCond <- file.path(outDir, cond)
  if (!file.exists(outDirCond)) {
    dir.create(outDirCond)
  }

  outputClusters <- set_names(rep(NA, length = getNumCells(objCOTAN)),
                              getCells(objCOTAN))

  iter <- 0L
  numClustersToRecluster <- 0L
  srat <- NULL

  repeat {
    outDirIter <- file.path(outDirCond, paste0("reclustering_", iter))
    if (!file.exists(outDirIter)) {
      dir.create(file.path(outDirIter))
    }

    logThis(paste0("In iteration ", iter, " "), logLevel = 1L, appendLF = FALSE)
    logThis(paste("the number of cells to re-cluster is",
                  sum(is.na(outputClusters)), "cells belonging to",
                  numClustersToRecluster, "clusters"), logLevel = 2L)

    #Step 1
    c(objSeurat, usedMaxResolution) %<-%
      seuratClustering(rawData = getRawData(objCOTAN)[, is.na(outputClusters)],
                       cond = cond, iter = iter,
                       initialResolution = initialResolution,
                       minNumClusters = numClustersToRecluster + 1L,
                       saveObj = saveObj, outDir = outDirIter)

    if (is_null(objSeurat)) {
      logThis(paste("NO new possible uniform clusters!",
                    "Unclustered cell left:", sum(is.na(outputClusters))),
              logLevel = 1L)
      break
    }

    if (saveObj && iter == 0L) {
      # save the Seurat object on the global raw data
      srat <- objSeurat
    }

    metaData <- objSeurat@meta.data

    rm(objSeurat)
    gc()

    # Step 2

    # The next lines check for each new cell cluster if it is homogeneous.
    # In cases they are not, save the corresponding cells in cellsToRecluster

    numClustersToRecluster <- 0L
    cellsToRecluster <- vector(mode = "character")
    sratClusters <- levels(factor(metaData[["seurat_clusters"]]))
    for (cl in sratClusters) {
      logThis("*", logLevel = 1L, appendLF = FALSE)
      logThis(paste0(" checking uniformity of cluster '", cl,
                     "' of ", length(sratClusters), " clusters"), logLevel = 2L)
      if (cl != "singleton") {
        cells <- rownames(metaData[metaData[["seurat_clusters"]] == cl, ])
        if (length(cells) < 10L) {
          logThis(paste("cluster", cl, "has too few cells:",
                        "will be reclustered!"), logLevel = 1L)
          cellsToRecluster <- c(cellsToRecluster, cells)
        } else {
          clusterIsUniform <-
            checkClusterUniformity(objCOTAN = objCOTAN,
                                   cluster = cl,
                                   cells = cells,
                                   cores = cores,
                                   GDIThreshold = GDIThreshold,
                                   saveObj = saveObj,
                                   outDir = outDirIter)[["isUniform"]]
          if (!clusterIsUniform) {
            logThis(paste("cluster", cl, "has too high GDI:",
                          "will be reclustered!"), logLevel = 1L)

            numClustersToRecluster <- numClustersToRecluster + 1L
            cellsToRecluster <- c(cellsToRecluster, cells)
          } else {
            logThis(paste("cluster", cl, "is uniform"), logLevel = 1L)
          }
        }
      }
    }
    logThis("", logLevel = 1L)
    logThis(paste("Found", length(sratClusters) - numClustersToRecluster,
                  "uniform and ", numClustersToRecluster,
                  "non-uniform clusters"), logLevel = 2L)

    if (numClustersToRecluster == length(sratClusters)) {
      warning("In iteration", iter, "no uniform clusters found!")
      # Another iteration can be attempted as the minimum number of clusters
      # will be higher. This happens unless the resolution already reached
      # its maximum. In the latter case we simply stop here.
      if (isTRUE(usedMaxResolution)) {
        break
      }
    }

    # Step 3: save the already uniform clusters keeping track of the iteration
    {
      flagInUniformCl <- !rownames(metaData) %in% cellsToRecluster
      goodClusters <- metaData[flagInUniformCl, "seurat_clusters", drop = FALSE]
      outputClusters[rownames(goodClusters)] <-
        paste0(str_pad(iter, width = 2L, pad = "0"), "_",
               str_pad(goodClusters[[1L]], width = 4L, pad = "0"))
    }

    if (sum(is.na(outputClusters)) != length(cellsToRecluster)) {
      warning("Some problems in cells reclustering")
      break
    }

    if (length(cellsToRecluster) <= 50L || iter >= maxIterations) {
      logThis(paste("NO new possible uniform clusters!",
                    "Unclustered cell left:", length(cellsToRecluster)),
              logLevel = 1L)
      break
    }
    rm(cellsToRecluster)

    iter <- iter + 1L
  } # End repeat

  logThis(paste("The final raw clusterization contains [",
                length(unique(outputClusters)), "] different clusters:",
                toString(unique(sort(outputClusters)))), logLevel = 3L)

  # replace the clusters' tags
  {
    clTags <- unique(sort(outputClusters))

    clTagsMap <- paste0(seq_along(clTags))
    clTagsMap <- factorToVector(niceFactorLevels(clTagsMap))
    clTagsMap <- set_names(clTagsMap, clTags)

    unclusteredCells <- is.na(outputClusters)
    outputClusters[!unclusteredCells] <-
      clTagsMap[outputClusters[!unclusteredCells]]
    outputClusters[unclusteredCells] <- "-1"
    outputClusters <- set_names(outputClusters, getCells(objCOTAN))
  }

  c(outputClusters, outputCoexDF) %<-%
    reorderClusterization(objCOTAN, clusters = outputClusters,
                          coexDF = NULL, keepMinusOne = TRUE,
                          distance = distance, hclustMethod = hclustMethod)

  if (saveObj) {
    clusterizationName <-
      paste0(as.roman(length(getClusterizations(objCOTAN)) + 1L))

    if (!setequal(rownames(srat@meta.data), names(outputClusters))) {
      warning("List of cells got corrupted")
      return(outputClusters)
    }

    srat@meta.data <-
      setColumnInDF(srat@meta.data,
                    colToSet = outputClusters[rownames(srat@meta.data)],
                    colName = paste0("COTAN_", clusterizationName))

    if (!dim(srat)[[2L]] == getNumCells(objCOTAN)) {
      warning("Number of cells got wrong")
      return(outputClusters)
    }

    logThis("Cluster, UMAP and Saving the Seurat dataset", logLevel = 2L)

    srat <- FindNeighbors(srat, dims = 1L:25L)
    srat <- FindClusters(srat, resolution = 0.5, algorithm = 2L)
    srat <- RunUMAP(srat, umap.method = "uwot",
                    metric = "cosine", dims = 1L:25L)

    saveRDS(srat, file.path(outDirCond, "Seurat_obj_with_cotan_clusters.RDS"))
  }

  rm(srat)
  gc()

  logThis("Creating cells' uniform clustering: DONE", logLevel = 2L)

  return(list("clusters" = factor(outputClusters), "coex" = outputCoexDF))
}
