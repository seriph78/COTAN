tm <- tempdir()
stopifnot(file.exists(tm))

library(zeallot)
library(stats)

test_that("Merge Uniform Cells Clusters", {

  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, calcCoex = FALSE, cores = 12L, saveObj = FALSE)

  clusters <- factor(readRDS(file.path(getwd(), "clusters1.RDS")))
  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  coexDF <- DEAOnClusters(obj)

  obj <- addClusterizationCoex(obj, clName = "clusters", coexDF = coexDF)

  expect_setequal(colnames(coexDF), levels(clusters))
  expect_identical(rownames(coexDF), getGenes(obj))

  lfcDF <- logFoldChangeOnClusters(obj, clusters = clusters)

  expect_setequal(colnames(lfcDF), clusters)
  expect_identical(rownames(lfcDF), getGenes(obj))
  expect_gte(min(colSums(lfcDF > 0.0)), 280L)
  expect_lte(max(colSums(lfcDF > 0.0)), 302L)
  expect_lt(max(abs(colMeans(lfcDF))), 0.06)

  method <- "bonferroni"

  pValDF <- pValueFromDEA(coexDF, numCells = getNumCells(obj), method = "none")
  adjPValDF <- pValueFromDEA(coexDF, numCells = getNumCells(obj),
                             method = method)

  expect_setequal(colnames(pValDF), levels(clusters))
  expect_identical(rownames(pValDF), getGenes(obj))
  expect_setequal(colnames(adjPValDF), levels(clusters))
  expect_identical(rownames(adjPValDF), getGenes(obj))

  coexDF_exp <- readRDS(file.path(getwd(), "coex.test.cluster1.RDS"))
  pValDF_exp <- readRDS(file.path(getwd(), "pval.test.cluster1.RDS"))

  expect_equal(coexDF[genes.names.test, ], coexDF_exp, tolerance = 1.0e-12)
  expect_equal(pValDF[genes.names.test, ], pValDF_exp, tolerance = 1.0e-12)

  GDIThreshold <- 1.5

  deltaExpression <- clustersDeltaExpression(obj)

  expect_identical(rownames(deltaExpression), getGenes(obj))
  expect_setequal(colnames(deltaExpression), levels(clusters))

  #primaryMarkers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
                       G2 = c("g-000300", "g-000330"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  e.df <- geneSetEnrichment(clustersCoex = coexDF, groupMarkers = groupMarkers)

  expect_identical(nrow(e.df), length(groupMarkers))
  expect_identical(ncol(e.df), ncol(coexDF) + 2L)
  expect_lte(max(e.df[, 1L:(ncol(e.df) - 2L)]), 1L)
  expect_gte(min(e.df[, 1L:(ncol(e.df) - 2L)]), 0L)
  expect_gte(min(e.df[["N. total"]] - e.df[["N. detected"]]), 0L)
  expect_equal(e.df[["N. total"]], vapply(groupMarkers, length, integer(1L)),
               ignore_attr = TRUE)

  c(mergedClusters, mergedCoexDF) %<-%
    mergeUniformCellsClusters(objCOTAN = obj, clusters = clusters, cores = 12L,
                              GDIThreshold = GDIThreshold,
                              distance = "cosine", hclustMethod = "ward.D2",
                              saveObj = TRUE, outDir = tm)

  expect_lt(nlevels(mergedClusters), nlevels(clusters))
  expect_setequal(colnames(mergedCoexDF), mergedClusters)
  expect_identical(rownames(mergedCoexDF), getGenes(obj))

  obj <- addClusterization(obj, clName = "merge", clusters = mergedClusters,
                           coexDF = mergedCoexDF, override = FALSE)

  expect_identical(reorderClusterization(obj),
                   list("clusters" = mergedClusters, "coex" = mergedCoexDF))

  #cluster_data <- readRDS(file.path(getwd(), "cluster_data_merged.RDS"))
  #expect_identical(mergedClusters[genes.names.test], cluster_data)

  mergedLfcDF <- logFoldChangeOnClusters(obj, clName = "merge")

  expect_setequal(colnames(mergedLfcDF), mergedClusters)
  expect_identical(rownames(mergedLfcDF), getGenes(obj))
  # with 2 clusters the changes are symmetric
  expect_identical(mergedLfcDF[[1]], -mergedLfcDF[[2]])

  for (cl in levels(mergedClusters)) {
    cellsToDrop <- names(clusters)[mergedClusters != cl]

    clObj <- dropGenesCells(obj, cells = cellsToDrop)

    clObj <- proceedToCoex(clObj, cores = 12L, saveObj = FALSE)

    GDIData <- calculateGDI(clObj)

    expect_lte(sum(GDIData[["GDI"]] >= GDIThreshold), 0.01 * nrow(GDIData))

    gc()
  }
})
