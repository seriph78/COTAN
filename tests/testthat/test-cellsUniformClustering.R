tm <- tempdir()
stopifnot(file.exists(tm))

library(zeallot)

test_that("Cell Uniform Clustering", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, calcCoex = FALSE,
                       cores = 12L, saveObj = TRUE, outDir = tm)

  GDIThreshold <- 1.5

  suppressWarnings({
    clusters <- cellsUniformClustering(obj, GDIThreshold = GDIThreshold,
                                       cores = 12L, saveObj = TRUE,
                                       outDir = tm)[["clusters"]]
  })

  gc()

  expect_identical(nlevels(clusters), 4L)

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  expect_equal(getClusterizationData(obj)[["clusters"]], clusters,
               ignore_attr = TRUE)

  firstCl <- clusters[[1L]]
  c(isUniform, fracAbove, lastPerc) %<-%
    checkClusterUniformity(obj, GDIThreshold = GDIThreshold,
                           cluster = paste0("Cluster_", firstCl),
                           cells = names(clusters)[clusters == firstCl],
                           saveObj = TRUE, outDir = tm)
  expect_true(isUniform)
  expect_lte(fracAbove, 0.01)
  expect_lte(lastPerc, GDIThreshold)

  #clusters_exp <- readRDS(file.path(getwd(), "clusters1.RDS"))

  #expect_identical(clusters, clusters_exp)

  clMarkersDF <- findClustersMarkers(obj)

  expect_identical(colnames(clMarkersDF), c("CL", "Gene", "Score", "pVal",
                                            "adjPVal", "DEA", "IsMarker"))
  expect_identical(nrow(clMarkersDF), 10L * 2L * length(unique(clusters)))
  expect_type(clMarkersDF[["Gene"]],     "character")
  expect_type(clMarkersDF[["IsMarker"]], "integer")
  expect_identical(sum(clMarkersDF[["IsMarker"]]), 0L)

  topGenesNum <- as.integer(substring(clMarkersDF[["Gene"]], 6L))
  highPos <- (1L:80L) %in% c(1L:10L, 21L:30L, 51L:60L, 71L:80L)
  expect_gt(min(topGenesNum[ highPos]), 480L)
  expect_lt(max(topGenesNum[!highPos]), 241L)

  primaryMarkers <-
    c("g-000010", "g-000020", "g-000030", "g-000300", "g-000330",
      "g-000510", "g-000530", "g-000550", "g-000570", "g-000590")
  clMarkersDF2 <- findClustersMarkers(obj, markers = primaryMarkers)

  expect_identical(colnames(clMarkersDF2), colnames(clMarkersDF))
  expect_identical(clMarkersDF2[, -7L], clMarkersDF[, -7L])
  expect_gt(sum(clMarkersDF2[["IsMarker"]]), 0L)

  clMarkersDF3 <- findClustersMarkers(obj, clusters = clusters)

  expect_identical(clMarkersDF3, clMarkersDF)

  ####################################

  # Test the low GDI (homogeneity) for each defined clusters

  ####################################

  for (cl in sample(unique(clusters[!is.na(clusters)]), size = 2L)) {
    print(cl)

    cellsToDrop <- which(clusters != cl)

    temp.obj <- dropGenesCells(objCOTAN = obj,
                               cells = getCells(obj)[cellsToDrop])

    temp.obj <- proceedToCoex(temp.obj, cores = 12L, saveObj = FALSE)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_lt(nrow(GDI_data[GDI_data[["GDI"]] >= GDIThreshold, ]),
              0.01 * nrow(GDI_data))
  }
})
