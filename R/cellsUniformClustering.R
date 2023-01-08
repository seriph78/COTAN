
#' get clusterization running the `Seurat` package
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
#' @param minNumClusters The minimum number of clusters expected from this
#'   clusterization. In cases it is not reached, it will increase the resolution
#'   of the clusterization.
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
#' @noRd
#'
seuratClustering <- function(rawData, cond, iter, minNumClusters,
                             saveObj, outDir) {
  logThis("Creating Seurat object: START", logLevel = 2)

  srat <- CreateSeuratObject(counts = as.data.frame(rawData),
                             project = paste0(cond, "_reclustering_", iter),
                             min.cells = if (iter == 0) 3 else 1,
                             min.features = if (iter == 0) 4 else 2)
  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat <- ScaleData(srat, features = rownames(srat))

  maxRows = nrow(srat@meta.data) - 1
  srat <- RunPCA(srat, features = VariableFeatures(object = srat),
                 npcs = min(50, maxRows))

  srat <- FindNeighbors(srat, dims = 1:min(25, maxRows))

  resolution <- 0.5
  repeat {
    srat <- FindClusters(srat, resolution = resolution, algorithm = 2)

    # The next lines are necessary to make cluster smaller while
    # the number of residual cells decrease and to stop clustering
    # if the algorithm gives too many singletons.
    if (minNumClusters < length(unique(srat$seurat_clusters)) ||
        resolution > 2) {
      break
    }

    logThis(paste("Number of clusters is too small.",
                  "Reclustering at higher resolution:", resolution),
            logLevel = 3)

    resolution <- resolution + 0.5
  }

  logThis(paste("Used resolution for Seurat clusterization is:", resolution),
          logLevel = 2)

  # disable annoying warning about Seurat::RunUMAP()
  options(Seurat.warn.umap.uwot = FALSE)
  srat <- RunUMAP(srat, umap.method = "uwot", metric = "cosine",
                  dims = 1:min(c(50, maxRows)))

  if (isTRUE(saveObj))
  {
    logThis(paste0("Creating PDF UMAP in file:",
                   file.path(outDir, "pdf_umap.pdf")), logLevel = 2)
    pdf(file.path(outDir, "pdf_umap.pdf"))

    if (iter == 0) {
      plot(DimPlot(srat, reduction = "umap", label = TRUE, group.by = "orig.ident"))
    }

    plot(DimPlot(srat, reduction = "umap",label = TRUE) +
         annotate(geom="text", x=0, y=30, color = "black",
                  label = paste0("Cells number: ", nrow(rawData), "\n",
                                 "Cl. resolution: ", resolution) ))

    dev.off()
    gc()
  }

  logThis("Creating Seurat object: DONE", logLevel = 2)

  return(list("SeuratObj" = srat, "UsedMaxResolution" = resolution > 2))
}


#' cellsUniformClustering
#'
#' @description The function finds a clusterizations that is **uniform** by GDI
#'
#' @details Once a preliminary *clusterization* is obtained from `Seurat`
#'   package methods, each cluster is checked for uniformity by looking at its
#'   GDI. If The GDI is too high, the cluster is deemed **non-uniform**. All
#'   cells from **non-uniform** clusters are then pooled together for another
#'   iteration of the entire process, until all clusters are deemed **uniform**.
#'
#' @param objCOTAN a `COTAN` object
#' @param cores number of cores (NB for windows system no more that 1 can be
#'   used)
#' @param maxIterations Max number of re-clustering iterations. Defaults to 25.
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output.
#'
#' @returns the new clusterization
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom stringr str_detect
#' @importFrom stringr str_pad
#'
#' @importFrom utils as.roman
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
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#' objCOTAN <- addClusterization(objCOTAN, clName = "clusters",
#'                               clusters = clusters)
#'
#' @rdname cellsUniformClustering
#'
cellsUniformClustering <-
  function(objCOTAN, cores = 1, maxIterations = 25,
           saveObj = TRUE, outDir = ".") {
  logThis("Creating cells' uniform clustering: START", logLevel = 2)

  cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

  outDirCond <- file.path(outDir, cond)
  if(!file.exists(outDirCond)){
    dir.create(outDirCond)
  }

  outputClusters <- set_names(rep(NA, length = getNumCells(objCOTAN)),
                              getCells(objCOTAN))

  iter <- 0
  numClustersToRecluster <- 0
  srat <- NULL

  repeat {
    outDirIter <- file.path(outDirCond, paste0("reclustering_", iter))
    if (!file.exists(outDirIter)) {
      dir.create(file.path(outDirIter))
    }

    logThis(paste0("In iteration ", iter, " "), logLevel = 1, appendLF = FALSE)
    logThis(paste("the number of cells to re-cluster is",
                  sum(is.na(outputClusters)), "cells belonging to",
                  numClustersToRecluster, "clusters"), logLevel = 2)

    #Step 1
    list[objSeurat, usedMaxResolution] <-
      seuratClustering(rawData = getRawData(objCOTAN)[,is.na(outputClusters)],
                       cond = cond, iter = iter,
                       minNumClusters = numClustersToRecluster + 1,
                       saveObj = saveObj, outDir = outDirIter)

    if (saveObj && iter == 0) {
      # save the Seurat object on the global raw data
      srat <- objSeurat
    }

    metaData <- objSeurat@meta.data

    rm(objSeurat)
    gc()

    # Step 2

    # The next lines check for each new cell cluster if it is homogeneous.
    # In cases they are not, save the corresponding cells in cellsToRecluster

    numClustersToRecluster <- 0
    cellsToRecluster <- vector(mode = "character")
    sratClusters <- unique(metaData[["seurat_clusters"]])
    for (cl in sratClusters) {
      logThis("*", logLevel = 1, appendLF = FALSE)
      logThis(paste0(" checking uniformity of cluster '", cl,
                     "' of ", length(sratClusters), " clusters" ), logLevel = 2)
      if (cl != "singleton") {
        cells <- rownames(metaData[metaData[["seurat_clusters"]] == cl,])
        if (length(cells) < 10) {
          logThis(paste("cluster", cl, "has too few cells: will be reclustered!"),
                  logLevel = 3)
          cellsToRecluster <- c(cellsToRecluster, cells)
        } else {
          clusterIsUniform <- checkClusterUniformity(objCOTAN = objCOTAN,
                                                     cluster = cl,
                                                     cells = cells,
                                                     cores = cores,
                                                     saveObj = saveObj,
                                                     outDir = outDirIter)
          if (!clusterIsUniform) {
            logThis(paste("cluster", cl, "has too high GDI: will be reclustered!"),
                    logLevel = 3)
            numClustersToRecluster <- numClustersToRecluster + 1
            cellsToRecluster <- c(cellsToRecluster, cells)
          } else {
            logThis(paste("cluster", cl, "is uniform"), logLevel = 3)
          }
        }
      }
    }
    logThis("", logLevel = 1)
    logThis(paste("Found", length(sratClusters) - numClustersToRecluster,
                  "uniform and ", numClustersToRecluster,
                  "non-uniform clusters"), logLevel = 2)

    if (numClustersToRecluster == length(sratClusters)) {
      warning(paste("In iteration", iter, "no uniform clusters found!"))
      # Another iteration can be attempted as the minimum number of clusters
      # will should higher. This happens unless the resolution already reached
      # its maximum. In the latter case we simply stop here.
      if (isTRUE(usedMaxResolution)) {
        break
      }
    }

    # Step 3: save the already uniform clusters keeping track of the iteration
    {
      flagInUniformCl <- !rownames(metaData) %in% cellsToRecluster
      goodClusters <- metaData[flagInUniformCl, ]["seurat_clusters"]
      outputClusters[rownames(goodClusters)] <-
        paste0(str_pad(iter, width = 2, pad = "0"), "_",
               str_pad(goodClusters[[1]], width = 4, pad = "0"))
    }

    if (sum(is.na(outputClusters)) != length(cellsToRecluster)) {
      warning("Some problems in cells reclustering")
      break
    }

    if (length(cellsToRecluster) <= 50 || iter >= maxIterations) {
      logThis(paste("NO new possible uniform clusters!",
                    "Unclustered cell left:", length(cellsToRecluster)), logLevel = 1)
      break
    }
    rm(cellsToRecluster)

    iter <- iter + 1
  } # End repeat

  logThis(paste("The final raw clusterization contains [",
                length(unique(outputClusters)), "] different clusters:",
                paste0(unique(outputClusters), collapse = ", ")), logLevel = 3)

  # replace the clusters' tags
  {
    clTags <- unique(sort(outputClusters))
    clTagsMap <- set_names(seq_along(clTags), clTags)

    unclusteredCells <- is.na(outputClusters)
    outputClusters[!unclusteredCells] <- clTagsMap[outputClusters[!unclusteredCells]]
    outputClusters <- set_names(as.character(outputClusters), getCells(objCOTAN))
    outputClusters[unclusteredCells] = "not_clustered"
  }

  if (saveObj) {
    clusterizationName <- paste0(as.roman(length(getClusterizations(objCOTAN)) + 1))

    if (!setequal(rownames(srat@meta.data), names(outputClusters))) {
      warning("List of cells got corrupted")
      return(outputClusters)
    }

    srat@meta.data <- setColumnInDF(srat@meta.data,
                                    colToSet = outputClusters[rownames(srat@meta.data)],
                                    colName = paste0("COTAN_", clusterizationName))

    if (!dim(srat)[2] == getNumCells(objCOTAN)) {
      warning("Number of cells got wrong")
      return(outputClusters)
    }

    logThis("Cluster, UMAP and Saving the Seurat dataset", logLevel = 2)

    srat <- FindNeighbors(srat, dims = 1:25)
    srat <- FindClusters(srat, resolution = 0.5, algorithm = 2)
    srat <- RunUMAP(srat, umap.method = "uwot", metric = "cosine", dims = 1:25)

    saveRDS(srat, file.path(outDirCond, "Seurat_obj_with_cotan_clusters.RDS"))
  }

  rm(srat)
  gc()

  logThis("Creating cells' uniform clustering: DONE", logLevel = 2)

  return(outputClusters)
}
