
test_that("Establish genes clusters", {
  data("test.dataset")
  objCOTAN <- COTAN(raw = test.dataset)
  objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)

  #primaryMarkers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
                       G2 = c("g-000300", "g-000330"),
                       G3 = c("g-000510", "g-000530", "g-000550", "g-000570", "g-000590"))

  c(secondaryMarkers, GCS, rankGenes) %<-%
    genesCoexSpace(objCOTAN = objCOTAN,
                   primaryMarkers = unlist(groupMarkers),
                   numGenesPerMarker = 11)

  expect_gt(length(secondaryMarkers), length(unlist(groupMarkers)))
  expect_equal(colnames(rankGenes), unlist(groupMarkers), ignore_attr = TRUE)
  expect_equal(rownames(rankGenes), secondaryMarkers)

  expect_equal(colnames(GCS), secondaryMarkers)
  expect_lte(max(abs(GCS)), 1)

  if (TRUE) {
    # saveRDS(GCS, file = "genes.coex.space.RDS")
    GCS_old <- readRDS(file.path(getwd(), "genes.coex.space.RDS"))
    expect_equal(GCS, GCS_old)
  }

  c(gSpace, plotEigen, pcaClusters, treePlot) %<-%
    establishGenesClusters(objCOTAN = objCOTAN,
                           groupMarkers = groupMarkers,
                           numGenesPerMarker = 11,
                           kCuts = 6,
                           distance = "cosine",
                           hclustMethod = "ward.D2")

  pcaExtraCols <- c("highlight", "hclust", "sec_markers",
                     "colors", "col_branches", "groupLabels")

  expect_s3_class(pcaClusters, "data.frame")
  expect_equal(ncol(pcaClusters), 16)
  expect_equal(colnames(pcaClusters), c(paste0("PC", c(1:10)), pcaExtraCols))

  if (FALSE) {
    # saveRDS(pcaClusters, "pca.clusters.RDS")
    pcaClustersExp <- readRDS(file.path(getwd(), "pca.clusters.RDS"))
    expect_equal(pcaClusters, pcaClustersExp)
  }

  expect_equal(nrow(gSpace), nrow(pcaClusters))
  expect_equal(GCS, gSpace)

  expect_s3_class(plotEigen, "ggplot")
  expect_equal(dim(plotEigen[["data"]]), c(10, 2))

  expect_s3_class(treePlot, "dendrogram")
})
