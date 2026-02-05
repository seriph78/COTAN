library(zeallot)

prevOptState <- options(parallelly.fork.enable = TRUE)

test_that("Establish genes clusters", {
  data("test.dataset")
  objCOTAN <- COTAN(raw = test.dataset)
  objCOTAN <- proceedToCoex(objCOTAN, cores = 6L,
                            optimizeForSpeed = TRUE, saveObj = FALSE)

  c(secondaryMarkers, gCS, rankGenes) %<-%
    genesCoexSpace(objCOTAN = objCOTAN,
                   primaryMarkers = "g-000300",
                   numGenesPerMarker = 5L)

  expect_gt(length(secondaryMarkers), 1L)
  expect_equal(colnames(rankGenes), "g-000300", ignore_attr = TRUE)
  expect_identical(rownames(rankGenes), secondaryMarkers)

  expect_identical(colnames(gCS), secondaryMarkers)
  expect_lte(max(abs(gCS)), 1L)

  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000138"),
                       G2 = c("g-000300", "g-000330", "g-000660"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  expect_warning(
    c(secondaryMarkers, gCS, rankGenes) %<-%
    genesCoexSpace(objCOTAN = objCOTAN,
                   primaryMarkers = unlist(groupMarkers),
                   numGenesPerMarker = 11L))

  expect_gt(length(secondaryMarkers), length(unlist(groupMarkers)))
  expect_equal(colnames(rankGenes), unlist(groupMarkers)[-6L],
               ignore_attr = TRUE)
  expect_identical(rownames(rankGenes), secondaryMarkers)

  expect_identical(colnames(gCS), secondaryMarkers)
  expect_lte(max(abs(gCS)), 1L)

  if (TRUE) {
    gcsOld <-
      as.matrix(readRDS(file.path(getwd(), "genes.coex.space.test.RDS")))
    expect_equal(gCS, gcsOld, tolerance = 5.0e-6)
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

  if (TRUE) {
    pcaClustersExp <- readRDS(file.path(getwd(), "pca.genes.clusters.test.RDS"))
    expect_identical(dim(pcaClusters), dim(pcaClustersExp))
    expect_identical(colnames(pcaClusters), colnames(pcaClustersExp))
    expect_identical(rownames(pcaClusters), rownames(pcaClustersExp))
    expect_identical(ncol(pcaClusters), 16L)
    expect_equal(abs(pcaClusters[, 1L:8L]),
                 abs(pcaClustersExp[, 1L:8L]), tolerance = 5.0e-5)
    expect_equal(abs(pcaClusters[, 9L:10L]),
                 abs(pcaClustersExp[, 9L:10L]), tolerance = 5.0e-4)
    expect_identical(pcaClusters[, 11L:13L], pcaClustersExp[, 11L:13L])
    expect_identical(pcaClusters[, 16L],     pcaClustersExp[, 16L])
  }

  expect_identical(nrow(gSpace), nrow(pcaClusters))
  expect_identical(gCS, gSpace)

  expect_s3_class(plotEigen, "ggplot")
  expect_identical(dim(plotEigen[["data"]]), c(10L, 2L))

  expect_s3_class(treePlot, "dendrogram")
})

options(prevOptState)
