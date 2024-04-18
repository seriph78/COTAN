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
                       cores = 6L, saveObj = TRUE, outDir = tm)

  GDIThreshold <- 1.46
  initialResolution <- 0.8
  suppressWarnings({
    clusters <- cellsUniformClustering(obj, GDIThreshold = GDIThreshold,
                                       initialResolution = initialResolution,
                                       cores = 6L, optimizeForSpeed = TRUE,
                                       deviceStr = "cuda", saveObj = TRUE,
                                       outDir = tm)[["clusters"]]
  })

  expect_true(file.exists(file.path(tm, "test", "reclustering",
                                    "partial_clusterization_1.csv")))
  expect_true(file.exists(file.path(tm, "test", "reclustering",
                                    "all_check_results_1.csv")))
  expect_true(file.exists(file.path(tm, "test", "split_check_results.csv")))

  gc()

  expect_identical(nlevels(clusters), 4L)

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  expect_equal(getClusters(obj), clusters, ignore_attr = TRUE)
  expect_identical(reorderClusterization(obj)[["clusters"]], clusters)

  firstCl <- clusters[[1L]]
  c(isUniform, fracAbove, firstPerc, clSize) %<-%
    checkClusterUniformity(obj, GDIThreshold = GDIThreshold,
                           clusterName = paste0("Cluster_", firstCl),
                           cells = names(clusters)[clusters == firstCl],
                           optimizeForSpeed = TRUE, deviceStr = "cpu",
                           saveObj = TRUE, outDir = tm)
  expect_true(isUniform)
  expect_lte(fracAbove, 0.01)
  expect_lte(firstPerc, GDIThreshold)
  expect_identical(clSize, sum(clusters == firstCl))

  clusters2 <- factor(clusters, levels = c(levels(clusters), "-1"))
  clusters2[1L:50L] <- "-1"
  coexDF2 <- DEAOnClusters(obj, clusters = clusters2)
  c(rClusters2, rCoexDF2) %<-%
    reorderClusterization(obj, reverse = TRUE, keepMinusOne = TRUE,
                          clusters = clusters2, coexDF = coexDF2)

  expect_identical(levels(rClusters2)[rClusters2[1L:50L]],
                   levels(clusters2)[clusters2[1L:50L]])
  expect_identical(rClusters2 == "-1", clusters2 == "-1")
  # this is an happenstance
  expect_identical(colnames(rCoexDF2), rev(levels(rClusters2)))

  c(clusters3, coexDF3) %<-%
    reorderClusterization(obj, useDEA = FALSE,
                          reverse = FALSE, keepMinusOne = TRUE,
                          clusters = clusters2, coexDF = coexDF2)

  expect_identical(levels(clusters3)[clusters3[1L:50L]],
                   levels(clusters2)[clusters2[1L:50L]])
  expect_identical(clusters3 == "-1", clusters2 == "-1")
  # this is an happenstance
  expect_identical(colnames(coexDF3)[-5L], levels(clusters3)[-1L])

  exactClusters <- set_names(rep(1L:2L, each = 600L), nm = getCells(obj))

  suppressWarnings({
    splitList <- cellsUniformClustering(obj, GDIThreshold = GDIThreshold,
                                        initialResolution = initialResolution,
                                        initialClusters = exactClusters,
                                        cores = 6L, optimizeForSpeed = FALSE,
                                        deviceStr = "cpu", saveObj = TRUE,
                                        outDir = tm)
  })

  expect_identical(splitList[["clusters"]], factor(exactClusters))

  clMarkersDF <- findClustersMarkers(obj)

  expect_identical(colnames(clMarkersDF), c("CL", "Gene", "Score", "adjPVal",
                                            "DEA", "IsMarker", "logFoldCh"))
  expect_identical(nrow(clMarkersDF), 10L * 2L * length(unique(clusters)))
  expect_type(clMarkersDF[["Gene"]],     "character")
  expect_type(clMarkersDF[["IsMarker"]], "integer")
  expect_identical(sum(clMarkersDF[["IsMarker"]]), 0L)
  expect_gt(min(clMarkersDF[["Score"]] * clMarkersDF[["DEA"]]), 0.0)
  expect_gt(min(clMarkersDF[["Score"]] * clMarkersDF[["logFoldCh"]]), 0.0)

  topGenesNum <- as.integer(substring(clMarkersDF[["Gene"]], 6L))
  highPos <- (1L:80L) %in% c(11L:20L, 31L:40L, 41L:50L, 61L:70L)
  expect_gt(min(topGenesNum[ highPos]), 480L)
  expect_lt(max(topGenesNum[!highPos]), 241L)

  primaryMarkers <-
    c("g-000010", "g-000020", "g-000030", "g-000300", "g-000330",
      "g-000510", "g-000530", "g-000550", "g-000570", "g-000590")
  clMarkersDF2 <- findClustersMarkers(obj, markers = primaryMarkers)

  expect_identical(colnames(clMarkersDF2), colnames(clMarkersDF))
  expect_identical(clMarkersDF2[, -6L], clMarkersDF[, -6L])
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

    temp.obj <- proceedToCoex(temp.obj, cores = 6L, saveObj = FALSE)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_lt(nrow(GDI_data[GDI_data[["GDI"]] >= GDIThreshold, ]),
              0.01 * nrow(GDI_data))
  }
})
