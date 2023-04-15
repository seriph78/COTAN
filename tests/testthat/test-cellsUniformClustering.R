tm = tempdir()
stopifnot(file.exists(tm))

test_that("Cell Uniform Clustering", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, cores = 12, saveObj = FALSE)

  GDIThreshold <- 1.5

  clusters <- cellsUniformClustering(obj, cores = 12,
                                     GDIThreshold = GDIThreshold,
                                     saveObj = FALSE, outDir = tm)

  gc()

  expect_equal(length(levels(factor(clusters))), 4)

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  expect_equal(getClusterizationData(obj)[["clusters"]], clusters,
               ignore_attr = TRUE)

  firstCl <- clusters[[1L]]
  expect_true(
    checkClusterUniformity(obj, GDIThreshold = GDIThreshold,
                           cluster = paste0("Cluster_", firstCl),
                           cells = names(clusters)[clusters == firstCl],
                           saveObj = TRUE, outDir = tm))

  #clusters_exp <- readRDS(file.path(getwd(),"clusters1.RDS"))

  #expect_equal(clusters, clusters)

  clMarkersDF <- findClustersMarkers(obj)

  expect_equal(colnames(clMarkersDF), c("CL", "Gene", "Score", "pVal",
                                        "adjPVal", "DEA", "IsMarker"))
  expect_equal(nrow(clMarkersDF), 10L * 2L * length(unique(clusters)))
  expect_type(clMarkersDF[["Gene"]],     "character")
  expect_type(clMarkersDF[["IsMarker"]], "integer")
  expect_equal(sum(clMarkersDF[["IsMarker"]]), 0)

  topGenesNum <- as.integer(substring(clMarkersDF[["Gene"]], 6))
  highPos <- c(1:80) %in% c(1:10, 31:40, 51:70)
  expect_gt(min(topGenesNum[ highPos]), 480)
  expect_lt(max(topGenesNum[!highPos]), 241)

  primaryMarkers <-
    c("g-000010", "g-000020", "g-000030", "g-000300", "g-000330",
      "g-000510", "g-000530", "g-000550", "g-000570", "g-000590")
  clMarkersDF2 <- findClustersMarkers(obj, markers = primaryMarkers)

  expect_equal(colnames(clMarkersDF2), colnames(clMarkersDF))
  expect_equal(clMarkersDF2[, -7], clMarkersDF[, -7])
  expect_gt(sum(clMarkersDF2[["IsMarker"]]), 0)

  clMarkersDF3 <- findClustersMarkers(obj, clusters = clusters)

  expect_equal(clMarkersDF3, clMarkersDF)

  ####################################

  # Test the low GDI (homogeneity) for each defined clusters

  ####################################

  for (cl in sample(unique(clusters[!is.na(clusters)]), size = 2)) {
    print(cl)

    cellsToDrop <- which(clusters != cl)

    temp.obj <- dropGenesCells(objCOTAN = obj,
                               cells = getCells(obj)[cellsToDrop])

    temp.obj <- proceedToCoex(temp.obj, cores = 12, saveObj = FALSE)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_lt(nrow(GDI_data[GDI_data[["GDI"]] >= GDIThreshold, ]),
              0.01 * nrow(GDI_data))
  }
})
