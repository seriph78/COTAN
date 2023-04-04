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

  expect_equal(getClusterizationData(obj)[["clusters"]], clusters, ignore_attr = TRUE)

  #clusters_exp <- readRDS(file.path(getwd(),"clusters1.RDS"))

  #expect_equal(clusters, clusters)

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

    expect_lt(nrow(GDI_data[GDI_data[["GDI"]] >= GDIThreshold, ]), 0.01 * nrow(GDI_data))
  }
})
