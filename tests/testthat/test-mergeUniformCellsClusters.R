tm <- tempdir()
stopifnot(file.exists(tm))

library(stringr)
library(zeallot)

options(parallelly.fork.enable = TRUE)

test_that("Merge Uniform Cells Clusters", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(objCOTAN = obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(objCOTAN = obj,
                       calcCoex = FALSE, cores = 6L, saveObj = FALSE)

  clusters <- factor(readRDS(file.path(getwd(), "split.clusters.test.RDS")))
  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))

  obj <- addClusterization(objCOTAN = obj,
                           clName = "clusters", clusters = clusters)

  coexDF <- DEAOnClusters(objCOTAN = obj)

  obj <- addClusterizationCoex(objCOTAN = obj,
                               clName = "clusters", coexDF = coexDF)

  expect_setequal(colnames(coexDF), levels(clusters))
  expect_identical(rownames(coexDF), getGenes(objCOTAN = obj))

  lfcDF <- logFoldChangeOnClusters(objCOTAN = obj,
                                   clusters = clusters,
                                   floorLambdaFraction = 0.01)

  expect_setequal(colnames(lfcDF), clusters)
  expect_identical(rownames(lfcDF), getGenes(objCOTAN = obj))
  expect_gte(min(colSums(lfcDF > 0.0)), 280L)
  expect_lte(max(colSums(lfcDF > 0.0)), 320L)
  expect_lt(max(colMeans(lfcDF)),  0.01)
  expect_gt(min(colMeans(lfcDF)), -0.12)

  adjustmentMethod <- "bonferroni"

  pValDF <- pValueFromDEA(coexDF, numCells = getNumCells(objCOTAN = obj),
                          adjustmentMethod = "none")
  adjPValDF <- pValueFromDEA(coexDF, numCells = getNumCells(objCOTAN = obj),
                             adjustmentMethod = adjustmentMethod)

  expect_setequal(colnames(pValDF), levels(clusters))
  expect_identical(rownames(pValDF), getGenes(objCOTAN = obj))
  expect_setequal(colnames(adjPValDF), levels(clusters))
  expect_identical(rownames(adjPValDF), getGenes(objCOTAN = obj))

  coexDF_exp <- readRDS(file.path(getwd(), "coex.clusters.test.RDS"))
  pValDF_exp <- readRDS(file.path(getwd(), "pvalues.clusters.test.RDS"))

  expect_equal(coexDF[genes.names.test, ], coexDF_exp, tolerance = 1.0e-12)
  expect_equal(pValDF[genes.names.test, ], pValDF_exp, tolerance = 1.0e-12)

  deltaExpression <- clustersDeltaExpression(objCOTAN = obj)

  expect_identical(rownames(deltaExpression), getGenes(objCOTAN = obj))
  expect_setequal(colnames(deltaExpression), levels(clusters))

  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
                       G2 = c("g-000300", "g-000330", "g-000660"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  e.df <- geneSetEnrichment(clustersCoex = coexDF, groupMarkers = groupMarkers)

  expect_identical(nrow(e.df), length(groupMarkers))
  expect_identical(ncol(e.df), ncol(coexDF) + 2L)
  expect_lte(max(e.df[, 1L:(ncol(e.df) - 2L)]), 1L)
  expect_gte(min(e.df[, 1L:(ncol(e.df) - 2L)]), 0L)
  expect_gte(min(e.df[["N. total"]] - e.df[["N. detected"]]), 0L)
  expect_equal(e.df[["N. total"]], lengths(groupMarkers), ignore_attr = TRUE)

  checker <- new("AdvancedGDIUniformityCheck")
  checkers <- list(checker, shiftCheckerThresholds(checker, 0.1))

  suppressWarnings({
    c(mergedClusters, mergedCoexDF) %<-%
      mergeUniformCellsClusters(
        objCOTAN = obj, clusters = clusters,
        checkers = checkers, allCheckResults = data.frame(),
        batchSize = 3L, distance = "cosine", hclustMethod = "ward.D2",
        cores = 6L, optimizeForSpeed = TRUE, deviceStr = "cuda",
        saveObj = TRUE, outDir = tm)
  })
  expect_true(file.exists(file.path(tm, "test", "leafs_merge",
                                    "merge_clusterization_1.csv")))
  expect_true(file.exists(file.path(tm, "test", "leafs_merge",
                                    "merge_clusterization_3.csv")))
  expect_true(file.exists(file.path(tm, "test", "leafs_merge",
                                    "all_check_results_1.csv")))
  expect_true(file.exists(file.path(tm, "test", "leafs_merge",
                                    "all_check_results_3.csv")))
  expect_true(file.exists(file.path(tm, "test", "merge_check_results.csv")))
  expect_true(file.exists(file.path(tm, "test", "merge_clusterization.csv")))

  allCheckDF <- read.csv(file.path(tm, "test", "leafs_merge",
                                   "all_check_results_3.csv"),
                         header = TRUE, row.names = 1L)
  expect_identical(nrow(allCheckDF), 1L)
  expect_identical(ncol(allCheckDF), 35L)

  allCheckRes <- dfToCheckers(allCheckDF)

  expect_length(allCheckRes, nrow(allCheckDF))
  expect_named(allCheckRes)
  expect_true(all(stringr::str_ends(names(allCheckRes), fixed("-merge"))))
  expect_true(all(vapply(allCheckRes, methods::is,
                         logical(1L), "AdvancedGDIUniformityCheck")))
  expect_false(allCheckRes[[1L]]@isUniform)
  expect_false(is.finite(allCheckRes[[1L]]@firstCheck@fractionBeyond))
  expect_true(is.finite(allCheckRes[[1L]]@firstCheck@quantileAtRatio))

  expect_identical(nlevels(mergedClusters), nlevels(clusters))
  expect_setequal(colnames(mergedCoexDF), mergedClusters)
  expect_identical(rownames(mergedCoexDF), getGenes(objCOTAN = obj))

  obj <- addClusterization(objCOTAN = obj,
                           clName = "merge", clusters = mergedClusters,
                           coexDF = mergedCoexDF, override = FALSE)

  expect_identical(reorderClusterization(objCOTAN = obj)[1L:2L],
                   list("clusters" = mergedClusters, "coex" = mergedCoexDF))

  mergedLfcDF <- logFoldChangeOnClusters(objCOTAN = obj, clName = "merge")

  expect_setequal(colnames(mergedLfcDF), mergedClusters)
  expect_identical(rownames(mergedLfcDF), getGenes(objCOTAN = obj))
  # with 2 clusters the changes are symmetric
  expect_equal(mergedLfcDF[[1L]], -mergedLfcDF[[2L]], tolerance = 1.0e-12)

  # Test the low GDI (homogeneity) for each defined clusters
  simpleChecker <- checkers[[2L]]@thirdCheck
  for (cl in levels(mergedClusters)) {
    print(paste("Tested cluster:", cl))

    cellsToDrop <- names(clusters)[mergedClusters != cl]

    tmpObj <- dropGenesCells(objCOTAN = obj, cells = cellsToDrop)

    suppressWarnings({
      tmpObj <- proceedToCoex(objCOTAN = tmpObj, cores = 6L, saveObj = FALSE)
    })

    GDIData <- calculateGDI(objCOTAN = tmpObj)

    expect_lte(sum(GDIData[["GDI"]] >= simpleChecker@GDIThreshold),
               simpleChecker@maxRatioBeyond * nrow(GDIData))
  }

  gc()
})
