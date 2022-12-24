
test_that("Establish genes clusters", {

  data("raw.dataset")
  objCOTAN <- COTAN(raw = raw.dataset)
  objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)

  #primaryMarkers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
  groupMarkers <- list(G1 = c("Pcbp2", "Snrpe", "Nfyb"), G2 = c("Prpf40a", "Ergic2"),
                       G3 = c("Ncl", "Cd47", "Macrod2", "Fth1", "Supt16"))

  list[secondaryMarkers, GCS, rankGenes] <-
    genesCoexSpace(objCOTAN = objCOTAN,
                   primaryMarkers = unlist(groupMarkers),
                   numGenesPerMarker = 11)

  expect_gt(length(secondaryMarkers), length(unlist(groupMarkers)))
  expect_equal(colnames(rankGenes), unlist(groupMarkers), ignore_attr = TRUE)
  expect_equal(rownames(rankGenes), secondaryMarkers)

  expect_equal(colnames(GCS), secondaryMarkers)
  expect_lte(max(abs(GCS)), 1)

  if (FALSE) {
    GCS_old <- readRDS(file.path(getwd(), "genes.coex.space.RDS"))
    expect_equal(GCS, GCS_old)
  }

  list[g.space, plot.eig, pca_clusters, tree_plot] <-
    establishGenesClusters(objCOTAN = objCOTAN,
                           groupMarkers = groupMarkers,
                           numGenesPerMarker = 11,
                           kCuts = 6,
                           distance = "cosine",
                           hclustMethod = "ward.D2")

  pcaExtraCols <- c("highlight", "hclust", "sec_markers",
                     "colors", "col_branches", "groupLabels")

  expect_s3_class(pca_clusters, "data.frame")
  expect_equal(ncol(pca_clusters), 16)
  expect_equal(colnames(pca_clusters), c(paste0("PC", c(1:10)), pcaExtraCols))

  if (FALSE) {
    pca_clusters_old <- readRDS(file.path(getwd(), "pca.clusters.RDS"))
    expect_equal(pca_clusters, pca_clusters_old)
  }

  expect_equal(nrow(g.space), nrow(pca_clusters))
  expect_equal(GCS, g.space)

  expect_s3_class(plot.eig, "ggplot")
  expect_equal(dim(plot.eig$data), c(10, 2))

  expect_s3_class(tree_plot, "dendrogram")
})
