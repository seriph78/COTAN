
test_that("Establish genes clusters", {
  data("test.dataset")
  objCOTAN <- COTAN(raw = test.dataset)
  objCOTAN <- proceedToCoex(objCOTAN, cores = 12L, saveObj = FALSE)

  #primaryMarkers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
                       G2 = c("g-000300", "g-000330"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  c(secondaryMarkers, GCS, rankGenes) %<-%
    genesCoexSpace(objCOTAN = objCOTAN,
                   primaryMarkers = unlist(groupMarkers),
                   numGenesPerMarker = 11L)

  expect_gt(length(secondaryMarkers), length(unlist(groupMarkers)))
  expect_equal(colnames(rankGenes), unlist(groupMarkers), ignore_attr = TRUE)
  expect_identical(rownames(rankGenes), secondaryMarkers)

  expect_identical(colnames(GCS), secondaryMarkers)
  expect_lte(max(abs(GCS)), 1L)

  if (TRUE) {
    # saveRDS(GCS, file = "genes.coex.space.RDS")
    GCS_old <- readRDS(file.path(getwd(), "genes.coex.space.RDS"))
    expect_identical(GCS, GCS_old)
  }

  c(gSpace, plotEigen, pcaClusters, treePlot) %<-%
    establishGenesClusters(objCOTAN = objCOTAN,
                           groupMarkers = groupMarkers,
                           numGenesPerMarker = 11L,
                           kCuts = 6L,
                           distance = "cosine",
                           hclustMethod = "ward.D2")

  pcaExtraCols <- c("highlight", "hclust", "sec_markers",
                     "colors", "col_branches", "groupLabels")

  expect_s3_class(pcaClusters, "data.frame")
  expect_identical(ncol(pcaClusters), 16L)
  expect_identical(colnames(pcaClusters),
                   c(paste0("PC", (1L:10L)), pcaExtraCols))

  if (FALSE) {
    # saveRDS(pcaClusters, "pca.clusters.RDS")
    pcaClustersExp <- readRDS(file.path(getwd(), "pca.clusters.RDS"))
    expect_identical(pcaClusters, pcaClustersExp)
  }

  expect_identical(nrow(gSpace), nrow(pcaClusters))
  expect_identical(GCS, gSpace)

  expect_s3_class(plotEigen, "ggplot")
  expect_identical(dim(plotEigen[["data"]]), c(10L, 2L))

  expect_s3_class(treePlot, "dendrogram")
})
