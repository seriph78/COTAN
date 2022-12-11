tm = tempdir()
stopifnot(file.exists(tm))

test_that("Cell Uniform Clustering", {
  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- clean(obj, calcExtraData = FALSE)[["objCOTAN"]]

  obj <- estimateDispersionBisection(obj, cores = 12)

  obj <- calculateCoex(obj)

  currentLevel <- setLoggingLevel(2)

  clusters <- cellsUniformClustering(obj, cores = 12,
                                     saveObj = FALSE, outDir = tm)

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  expect_equal(getClusterizationData(obj)[["clusters"]], clusters, ignore_attr = TRUE)

  #clusters_exp <- readRDS(file.path(getwd(),"clusters1.RDS"))

  #expect_equal(clusters, clusters)

  ####################################

  # Test the low GDI (homogeneity) for each defined clusters

  ####################################

  setLoggingLevel(currentLevel)

  for (cl in sample(unique(clusters[!is.na(clusters)]), size = 5)) {
    print(cl)

    cellsToDrop <- which(clusters != cl)

    temp.obj <- dropGenesCells(objCOTAN = obj,
                               cells = getCells(obj)[cellsToDrop])

    temp.obj <- clean(temp.obj, calcExtraData = FALSE)[["objCOTAN"]]

    temp.obj <- estimateDispersionBisection(temp.obj, cores = 12)
    temp.obj <- calculateCoex(temp.obj)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_lt(nrow(GDI_data[GDI_data[["GDI"]] >= 1.5, ]), 0.01 * nrow(GDI_data))
  }
})
