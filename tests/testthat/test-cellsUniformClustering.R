tm <- tempdir()
stopifnot(file.exists(tm))

library(zeallot)

options(parallelly.fork.enable = TRUE)

test_that("Cell Uniform Clustering", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(objCOTAN = obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(objCOTAN = obj, calcCoex = TRUE,
                       cores = 6L, saveObj = TRUE, outDir = tm)

  initialResolution <- 1.3
  c(sClusters, cellsRDM, resolution, usedMaxResolution) %<-%
    seuratClustering(objCOTAN = obj, initialResolution = initialResolution,
                     minNumClusters = 3L, useCoexEigen = TRUE,
                     dataMethod = "AdjBin", numReducedComp = 25L,
                     genesSel = "", numGenes = -1L) #irrelevant inputs

  expect_identical(nlevels(sClusters), 4L)
  expect_identical(as.vector(table(sClusters)), c(321L, 284L, 214L, 181L))
  expect_identical(dim(cellsRDM), c(1000L, 25L + 15L))
  expect_identical(resolution, 1.3)
  expect_identical(usedMaxResolution, FALSE)

  # Make it a less strict check as it is only for testing
  checker <- new("AdvancedGDIUniformityCheck")
  checker <- shiftCheckerThresholds(checker, 0.1)

  suppressWarnings({
    splitData <-
      cellsUniformClustering(objCOTAN = obj,
                             checker = checker,
                             initialResolution = initialResolution,
                             dataMethod = "LogLikelihood",
                             useCoexEigen = TRUE,
                             genesSel = "HGDI",
                             numGenes = 2000L,
                             numReducedComp = 50L,
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

  clusters <- splitData[["clusters"]]
  coexDF <- splitData[["coex"]]
  expect_identical(nlevels(clusters), 4L)
  expect_setequal(levels(clusters), colnames(coexDF))

  obj <- addClusterization(objCOTAN = obj, clName = "clusters",
                           clusters = clusters, coexDF = coexDF)

  expect_equal(getClusters(objCOTAN = obj), clusters, ignore_attr = TRUE)
  expect_identical(reorderClusterization(objCOTAN = obj)[["clusters"]],
                   clusters)

  firstCl <- clusters[[1L]]

  suppressWarnings({
      checkerRes <- checkClusterUniformity(
        objCOTAN = obj, checker = checker,
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
                   checker@firstCheck@GDIThreshold)
  expect_identical(checkerRes@secondCheck@maxRatioBeyond,
                   checker@secondCheck@maxRatioBeyond)
  expect_identical(checkerRes@thirdCheck@maxRatioBeyond,
                   checker@thirdCheck@maxRatioBeyond)
  expect_identical(checkerRes@fourthCheck@maxRankBeyond,
                   checker@fourthCheck@maxRankBeyond)

  clusters2 <- factor(clusters, levels = c("-1", levels(clusters)))
  clusters2[1L:50L] <- "-1"
  coexDF2 <- DEAOnClusters(objCOTAN = obj, clusters = clusters2)
  reorderRes2 <-
    reorderClusterization(objCOTAN = obj, reverse = TRUE, keepMinusOne = TRUE,
                          clusters = clusters2, coexDF = coexDF2)
  rClusters2 <- reorderRes2[["clusters"]]

  expect_identical(levels(rClusters2)[rClusters2[1L:50L]],
                   levels(clusters2)[clusters2[1L:50L]])
  expect_identical(rClusters2 == "-1", clusters2 == "-1")
  # this is an happenstance
  expect_identical(colnames(reorderRes2[["coex"]]), rev(levels(rClusters2)))
  expect_identical(reorderRes2[["permMap"]],
                   set_names(paste0(c(4L:1L, -1L)),
                             nm = paste0(c(1L:4L, -1L))))

  clusters3 <- factor(clusters, levels = c(levels(clusters), "-1"))
  clusters3[51L:100L] <- "-1"
  reorderRes3 <-
    reorderClusterization(objCOTAN = obj, useDEA = FALSE,
                          reverse = FALSE, keepMinusOne = TRUE,
                          clusters = clusters3, coexDF = coexDF2)
  clusters3 <- reorderRes3[["clusters"]]

  expect_identical(levels(clusters3)[clusters3[51L:100L]],
                   levels(clusters2)[clusters2[1L:50L]])
  expect_identical((clusters3 == "-1")[51L:150L], (clusters2 == "-1")[1L:100L])
  # this is an happenstance
  expect_identical(colnames(reorderRes3[["coex"]])[-5L], levels(clusters3)[-1L])
  expect_identical(reorderRes3[["permMap"]],
                   set_names(paste0(c(1L:4L, -1L)),
                             nm = paste0(c(1L:4L, -1L))))

  clSize <- getNumCells(obj) / 2L
  exactClusters <- set_names(rep(1L:2L, each = clSize), nm = getCells(obj))

  suppressWarnings({
    splitData2 <-
      cellsUniformClustering(objCOTAN = obj,
                             checker = checker,
                             initialResolution = initialResolution,
                             initialClusters = exactClusters,
                             cores = 6L,
                             optimizeForSpeed = TRUE,
                             deviceStr = "cuda",
                             saveObj = TRUE,
                             outDir = tm)
  })

  expect_identical(splitData2[["clusters"]], factor(exactClusters))

  clMarkersDF <- findClustersMarkers(objCOTAN = obj)

  expect_identical(colnames(clMarkersDF),
                   c("CL", "Gene", "DEA", "adjPVal", "IsMarker", "logFoldCh"))
  expect_identical(nrow(clMarkersDF), 10L * 2L * length(unique(clusters)))
  expect_type(clMarkersDF[["Gene"]],     "character")
  expect_type(clMarkersDF[["IsMarker"]], "integer")
  expect_identical(sum(clMarkersDF[["IsMarker"]]), 0L)
  expect_gt(min(clMarkersDF[["DEA"]] * clMarkersDF[["logFoldCh"]]), 0.0)

  topGenesNum <- as.integer(substring(clMarkersDF[["Gene"]], 6L))
  expect_gt(min(topGenesNum), 150L)

  primaryMarkers <-
    c("g-000010", "g-000020", "g-000138", "g-000300", "g-000330",
      "g-000510", "g-000530", "g-000550", "g-000579", "g-000590")
  clMarkersDF2 <- findClustersMarkers(objCOTAN = obj, markers = primaryMarkers)

  expect_identical(colnames(clMarkersDF2), colnames(clMarkersDF))
  expect_identical(clMarkersDF2[, -5L], clMarkersDF[, -5L])
  expect_gt(sum(clMarkersDF2[["IsMarker"]]), 0L)

  clMarkersDF3 <- findClustersMarkers(objCOTAN = obj, clusters = clusters)

  expect_identical(clMarkersDF3, clMarkersDF)

  # Test cluster/gene contingency tables
  contingencyTables <-
    clusterGeneContingencyTables(objCOTAN = obj, gene = "g-000138",
                                 cells = toClustersList(exactClusters)[[1L]])
  observedCT <- contingencyTables[["observed"]]
  expectedCT <- contingencyTables[["expected"]]

  expect_identical(rownames(observedCT), rownames(expectedCT))
  expect_identical(rownames(observedCT), c("g-000138.yes", "g-000138.no"))
  expect_identical(colSums(observedCT),
                   set_names(c(clSize, clSize), c("cells.in", "cells.out")))
  expect_identical(colSums(expectedCT),
                   set_names(c(clSize, clSize), c("cells.in", "cells.out")))
  expect_equal(rowSums(observedCT), rowSums(expectedCT), tolerance = 5e-6)
  expect_identical(observedCT < expectedCT,
                   matrix(c(TRUE, FALSE, FALSE, TRUE), nrow = 2L,
                          dimnames = dimnames(observedCT)))
  expect_identical(observedCT[1L, 1L],
                   sum(getZeroOneProj(objCOTAN = obj)["g-000138", 1L:clSize]))
  expect_equal(expectedCT[2L, 1L],
               sum(getProbabilityOfZero(objCOTAN = obj)["g-000138", 1L:clSize]),
               tolerance = 5e-6)

  # Test the low GDI (homogeneity) for each defined clusters
  simpleChecker <- checker@thirdCheck
  for (cl in levels(clusters)[[1L]]) {
    print(paste("Tested cluster:", cl))

    cellsToDrop <- names(clusters)[clusters != cl]

    tmpObj <- dropGenesCells(objCOTAN = obj, cells = cellsToDrop)

    suppressWarnings({
      tmpObj <- proceedToCoex(objCOTAN = tmpObj, cores = 6L, saveObj = FALSE)
    })

    GDIData <- calculateGDI(objCOTAN = tmpObj)

    expect_lte(sum(GDIData[["GDI"]] >= simpleChecker@GDIThreshold),
               simpleChecker@maxRatioBeyond * nrow(GDIData))
  }

  rm(obj, splitData, splitData2)

  gc()
})
