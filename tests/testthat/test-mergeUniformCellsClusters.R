
test_that("Merge Uniform Cells Clusters", {

  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- proceedToCoex(obj, cores = 12, saveObj = FALSE)

  clusters <- readRDS(file.path(getwd(), "clusters1.RDS"))
  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  c(coexDF, pValDF) %<-% DEAOnClusters(obj)

  obj <- addClusterizationCoex(obj, clName = "clusters", coexDF = coexDF)

  expect_equal(colnames(coexDF), unique(clusters))
  expect_equal(rownames(coexDF), getGenes(obj))

  expect_equal(colnames(pValDF), unique(clusters))
  expect_equal(rownames(pValDF), getGenes(obj))

  coexDF_exp <- readRDS(file.path(getwd(), "coex.test.cluster1.RDS"))
  pValDF_exp <- readRDS(file.path(getwd(), "pval.test.cluster1.RDS"))

  expect_equal(coexDF[genes.names.test, ], coexDF_exp)
  expect_equal(pValDF[genes.names.test, ], pValDF_exp)


  deltaExpression <- clustersDeltaExpression(obj)

  expect_equal(rownames(deltaExpression), getGenes(obj))
  expect_equal(colnames(deltaExpression), unique(clusters))

  #primaryMarkers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
  groupMarkers <- list(G1 = c("Pcbp2", "Snrpe", "Nfyb"), G2 = c("Prpf40a", "Ergic2"),
                       G3 = c("Ncl", "Cd47", "Macrod2", "Fth1", "Supt16"))

  e.df <- geneSetEnrichment(clustersCoex = coexDF, groupMarkers = groupMarkers)

  expect_equal(nrow(e.df), length(groupMarkers))
  expect_equal(ncol(e.df), ncol(coexDF) + 2)
  expect_lte(max(e.df[,1:(ncol(e.df)-2)]), 1)
  expect_gte(min(e.df[,1:(ncol(e.df)-2)]), 0)
  expect_gte(min(e.df[["N. total"]] - e.df[["N. detected"]]), 0)
  expect_equal(e.df[["N. total"]], sapply(groupMarkers, length), ignore_attr = TRUE)

  c(mergedClusters, mergedCoexDF, mergedPValueDF) %<-%
    mergeUniformCellsClusters(objCOTAN = obj, clusters = clusters, cores = 12,
                              distance = "cosine", hclustMethod = "ward.D2",
                              saveObj = FALSE)

  expect_lt(length(unique(mergedClusters)), length(unique(clusters)))
  expect_setequal(mergedClusters, colnames(mergedCoexDF))
  expect_setequal(colnames(mergedCoexDF), colnames(mergedPValueDF))

  #cluster_data <- readRDS(file.path(getwd(), "cluster_data_merged.RDS"))
  #expect_equal(mergedClusters[genes.names.test], cluster_data)

  for (cl in unique(mergedClusters)) {
    cellsToDrop <- names(clusters)[mergedClusters != cl]

    clObj <- dropGenesCells(obj, cells = cellsToDrop)

    clObj <- proceedToCoex(clObj, cores = 12, saveObj = FALSE)

    GDIData <- calculateGDI(clObj)

    expect_lte(sum(GDIData[["GDI"]] >= 1.5), 0.01 * nrow(GDIData))

    gc()
  }
})
