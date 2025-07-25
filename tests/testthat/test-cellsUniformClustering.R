tm <- tempdir()
stopifnot(file.exists(tm))

library(zeallot)

options(parallelly.fork.enable = TRUE)

test_that("Cell Uniform Clustering", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, calcCoex = TRUE,
                       cores = 6L, saveObj = TRUE, outDir = tm)

  initialResolution <- 0.8
  c(sClusters, cellsRDM, resolution, usedMaxResolution) %<-%
    seuratClustering(obj, initialResolution = initialResolution,
                     minNumClusters = 3L, useCoexEigen = TRUE,
                     dataMethod = "SignLogL", numReducedComp = 25L,
                     genesSel = "", numGenes = -1L) #irrelevant inputs

  expect_identical(nlevels(sClusters), 5L)
  expect_identical(as.vector(table(sClusters)), c(360L, 252L, 213L, 192L, 183L))
  expect_identical(dim(cellsRDM), c(1200L, 25L + 15L))
  expect_identical(resolution, 1.3)
  expect_identical(usedMaxResolution, FALSE)

  GDIThreshold <- 1.46
  suppressWarnings({
    c(clusters, coexDF) %<-%
      cellsUniformClustering(obj, GDIThreshold = GDIThreshold,
                             initialResolution = initialResolution,
                             dataMethod = "LogNormalized", useCoexEigen = FALSE,
                             genesSel = "HVG_Seurat", numGenes = 2000L,
                             numReducedComp = 25L,
                             cores = 6L, optimizeForSpeed = TRUE,
                             deviceStr = "cuda", saveObj = TRUE, outDir = tm)
  })

  expect_true(file.exists(file.path(tm, "test", "reclustering",
                                    "pdf_umap_1.pdf")))
  expect_true(file.exists(file.path(tm, "test", "reclustering",
                                    "partial_clusterization_1.csv")))
  expect_true(file.exists(file.path(tm, "test", "reclustering",
                                    "all_check_results_1.csv")))
  expect_true(file.exists(file.path(tm, "test", "split_check_results.csv")))
  expect_true(file.exists(file.path(tm, "test", "split_clusterization.csv")))

  gc()

  expect_identical(nlevels(clusters), 4L)
  expect_setequal(levels(clusters), colnames(coexDF))

  obj <- addClusterization(obj, clName = "clusters",
                           clusters = clusters, coexDF = coexDF)

  expect_equal(getClusters(obj), clusters, ignore_attr = TRUE)
  expect_identical(reorderClusterization(obj)[["clusters"]], clusters)

  firstCl <- clusters[[1L]]
  advChecker <- new("AdvancedGDIUniformityCheck")

  suppressWarnings({
      checkerRes <- checkClusterUniformity(
        objCOTAN = obj, checker = advChecker,
        clusterName = paste0("Cluster_", firstCl),
        cells = names(clusters)[clusters == firstCl],
        optimizeForSpeed = TRUE, deviceStr = "cuda",
        saveObj = TRUE, outDir = tm)
  })

  expect_true(checkerRes@isUniform)
  expect_lte(checkerRes@firstCheck@fractionBeyond,
             checkerRes@firstCheck@maxRatioBeyond)
  expect_lte(checkerRes@firstCheck@quantileAtRatio,
             checkerRes@firstCheck@GDIThreshold)
  expect_true(((checkerRes@secondCheck@fractionBeyond >=
                  checkerRes@secondCheck@maxRatioBeyond) &&
               (checkerRes@thirdCheck@fractionBeyond <=
                  checkerRes@thirdCheck@maxRatioBeyond)) ||
              (checkerRes@fourthCheck@thresholdRank <=
                 checkerRes@fourthCheck@maxRankBeyond))
  expect_true(((checkerRes@secondCheck@quantileAtRatio >=
                 checkerRes@secondCheck@GDIThreshold) &&
               (checkerRes@thirdCheck@quantileAtRatio <=
                 checkerRes@thirdCheck@GDIThreshold)) ||
              (checkerRes@fourthCheck@quantileAtRank <=
                 checkerRes@fourthCheck@GDIThreshold))
  expect_identical(checkerRes@clusterSize, sum(clusters == firstCl))
  expect_identical(checkerRes@firstCheck@GDIThreshold,
                   advChecker@firstCheck@GDIThreshold)
  expect_identical(checkerRes@secondCheck@maxRatioBeyond,
                   advChecker@secondCheck@maxRatioBeyond)
  expect_identical(checkerRes@thirdCheck@maxRatioBeyond,
                   advChecker@thirdCheck@maxRatioBeyond)
  expect_identical(checkerRes@fourthCheck@maxRankBeyond,
                   advChecker@fourthCheck@maxRankBeyond)

  clusters2 <- factor(clusters, levels = c("-1", levels(clusters)))
  clusters2[1L:50L] <- "-1"
  coexDF2 <- DEAOnClusters(obj, clusters = clusters2)
  c(rClusters2, rCoexDF2, permMap2) %<-%
    reorderClusterization(obj, reverse = TRUE, keepMinusOne = TRUE,
                          clusters = clusters2, coexDF = coexDF2)

  expect_identical(levels(rClusters2)[rClusters2[1L:50L]],
                   levels(clusters2)[clusters2[1L:50L]])
  expect_identical(rClusters2 == "-1", clusters2 == "-1")
  # this is an happenstance
  expect_identical(colnames(rCoexDF2), rev(levels(rClusters2)))
  expect_identical(permMap2, set_names(paste0(c(2L, 1L, 4L, 3L, -1L)),
                                       nm = paste0(c(1L:4L, -1L))))

  clusters3 <- factor(clusters, levels = c(levels(clusters), "-1"))
  clusters3[51L:100L] <- "-1"
  c(clusters3, coexDF3, permMap3) %<-%
    reorderClusterization(obj, useDEA = FALSE,
                          reverse = FALSE, keepMinusOne = TRUE,
                          clusters = clusters3, coexDF = coexDF2)

  expect_identical(levels(clusters3)[clusters3[51L:100L]],
                   levels(clusters2)[clusters2[1L:50L]])
  expect_identical((clusters3 == "-1")[51L:150L], (clusters2 == "-1")[1L:100L])
  # this is an happenstance
  expect_identical(colnames(coexDF3)[-5L], levels(clusters3)[-1L])
  expect_identical(permMap3, set_names(paste0(c(2L, 3L, 1L, 4L, -1L)),
                                       nm = paste0(c(1L:4L, -1L))))

  exactClusters <- set_names(rep(1L:2L, each = 600L), nm = getCells(obj))

  suppressWarnings({
    splitList <- cellsUniformClustering(
      obj, checker = shiftCheckerThresholds(advChecker, 0.1),
      initialResolution = initialResolution,
      initialClusters = exactClusters,
      cores = 6L, optimizeForSpeed = FALSE, deviceStr = "cpu",
      genesSel = "HGDI", saveObj = TRUE, outDir = tm)
  })

  expect_identical(splitList[["clusters"]], factor(exactClusters))

  clMarkersDF <- findClustersMarkers(obj)

  expect_identical(colnames(clMarkersDF),
                   c("CL", "Gene", "DEA", "adjPVal", "IsMarker", "logFoldCh"))
  expect_identical(nrow(clMarkersDF), 10L * 2L * length(unique(clusters)))
  expect_type(clMarkersDF[["Gene"]],     "character")
  expect_type(clMarkersDF[["IsMarker"]], "integer")
  expect_identical(sum(clMarkersDF[["IsMarker"]]), 0L)
  expect_gt(min(clMarkersDF[["DEA"]] * clMarkersDF[["logFoldCh"]]), 0.0)

  topGenesNum <- as.integer(substring(clMarkersDF[["Gene"]], 6L))
  highPos <- (1L:80L) %in% c(11L:20L, 31L:40L, 41L:50L, 61L:70L)
  expect_gt(min(topGenesNum[ highPos]), 480L)
  expect_lt(max(topGenesNum[!highPos]), 241L)

  primaryMarkers <-
    c("g-000010", "g-000020", "g-000030", "g-000300", "g-000330",
      "g-000510", "g-000530", "g-000550", "g-000570", "g-000590")
  clMarkersDF2 <- findClustersMarkers(obj, markers = primaryMarkers)

  expect_identical(colnames(clMarkersDF2), colnames(clMarkersDF))
  expect_identical(clMarkersDF2[, -5L], clMarkersDF[, -5L])
  expect_gt(sum(clMarkersDF2[["IsMarker"]]), 0L)

  clMarkersDF3 <- findClustersMarkers(obj, clusters = clusters)

  expect_identical(clMarkersDF3, clMarkersDF)

  # Test cluster/gene contingency tables
  c(observedCT, expectedCT) %<-%
    clusterGeneContingencyTables(objCOTAN = obj, gene = "g-000200",
                                 cells = toClustersList(exactClusters)[[1L]])

  expect_identical(rownames(observedCT), rownames(expectedCT))
  expect_identical(rownames(observedCT), c("g-000200.yes", "g-000200.no"))
  expect_identical(colSums(observedCT),
                   set_names(c(600.0, 600.0), c("cells.in", "cells.out")))
  expect_identical(colSums(expectedCT),
                   set_names(c(600.0, 600.0), c("cells.in", "cells.out")))
  expect_equal(rowSums(observedCT), rowSums(expectedCT), tolerance = 5e-6)
  expect_identical(observedCT < expectedCT,
                   matrix(c(TRUE, FALSE, FALSE, TRUE), nrow = 2L,
                          dimnames = dimnames(observedCT)))
  expect_identical(observedCT[1L, 1L],
                   sum(getZeroOneProj(obj)["g-000200", 1L:600L]))
  expect_equal(expectedCT[2L, 1L],
               sum(getProbabilityOfZero(obj)["g-000200", 1L:600L]),
               tolerance = 5e-6)

  # Test the low GDI (homogeneity) for each defined clusters

  for (cl in sample(unique(clusters[!is.na(clusters)]), size = 2L)) {
    print(paste("Tested cluster:", cl))

    cellsToDrop <- which(clusters != cl)

    temp.obj <- dropGenesCells(objCOTAN = obj,
                               cells = getCells(obj)[cellsToDrop])

    suppressWarnings({
      temp.obj <- proceedToCoex(temp.obj, cores = 6L, saveObj = FALSE)
    })
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_lte(nrow(GDI_data[GDI_data[["GDI"]] >= GDIThreshold, ]),
               0.01 * nrow(GDI_data))
  }
})
