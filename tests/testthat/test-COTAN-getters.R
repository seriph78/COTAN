
test_that("COTAN getters", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")
  obj <- clean(obj)
  obj <- estimateDispersionBisection(obj)
  obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = FALSE)
  obj <- calculateCoex(obj, actOnCells = TRUE,  optimizeForSpeed = TRUE)

  obj <- addClusterization(obj, clName = "Test",
                           clusters = rep(c(1, 2), 10))
  obj <- addClusterization(obj, clName = "Test2",
                           clusters = rep(c(2, 1), 10))

  metaInfo <- c("V", "10X", "Test", "20", "TRUE", "TRUE",
                paste0(10/55), paste0(0))

  expect_equal(getRawData(obj), as(as(raw, "dMatrix"), "sparseMatrix"))
  expect_equal(getNumGenes(obj), 10)
  expect_equal(getNumCells(obj), 20)
  expect_equal(getGenes(obj), head(LETTERS, getNumGenes(obj)))
  expect_equal(getCells(obj), head(letters, getNumCells(obj)))
  expect_equal(getZeroOneProj(obj), sign(getRawData(obj)))
  expect_equal(getCellsSize(obj), colSums(getRawData(obj)))
  expect_equal(getNormalizedData(obj), t(t(getRawData(obj)) * (1/getNu(obj))))
  expect_equal(getMetadataDataset(obj)[[1]], datasetTags()[c(1:8)],
               ignore_attr = TRUE)
  expect_equal(getMetadataDataset(obj)[[2]], metaInfo)
  expect_setequal(colnames(getMetadataGenes(obj)),
                  c("lambda", "hkGenes", "dispersion"))
  expect_equal(rownames(getMetadataGenes(obj)), getGenes(obj))
  expect_setequal(colnames(getMetadataCells(obj)),
                  c("nu", "feCells", names(getClustersCoex(obj))))
  expect_equal(rownames(getMetadataCells(obj)), getCells(obj))
  expect_equal(length(getClustersCoex(obj)), 2)
  expect_equal(names(getClustersCoex(obj)),
               paste0("CL_", getClusterizations(obj)))
  expect_equal(getClustersCoex(obj)[["CL_Test"]], data.frame())
  expect_equal(getNu(obj), getMetadataCells(obj)[["nu"]],
               ignore_attr = TRUE)
  expect_equal(getLambda(obj), getMetadataGenes(obj)[["lambda"]],
               ignore_attr = TRUE)
  expect_equal(getDispersion(obj), getMetadataGenes(obj)[["dispersion"]],
               ignore_attr = TRUE)
  expect_equal(flagNotHousekeepingGenes(obj),
               c(FALSE, rep(TRUE, getNumGenes(obj) - 1)))
  expect_equal(getHousekeepingGenes(obj), c(LETTERS[1]))
  expect_equal(flagNotFullyExpressedCells(obj), rep(TRUE, getNumCells(obj)))
  expect_equal(getFullyExpressedCells(obj), vector(mode = "character"))
  expect_equal(dim(getGenesCoex(obj)),
               as.integer(c(getNumGenes(obj), getNumGenes(obj))))
  expect_equal(dim(getCellsCoex(obj)),
               as.integer(c(getNumCells(obj), getNumCells(obj))))
  expect_equal(getClusterizations(obj), c("Test", "Test2"))
  expect_equal(getClusterizationName(obj), "Test2")
  expect_equal(getClusterizationName(obj, clName = "Test", keepPrefix = TRUE),
               "CL_Test")
  expect_setequal(names(getClusterizationData(obj)), c("coex","clusters"))
  expect_equal(getClusterizationData(obj)[["clusters"]],
               getMetadataCells(obj)[["CL_Test2"]], ignore_attr = TRUE)
  expect_equal(getClusterizationData(obj)[["coex"]], data.frame())
  expect_equal(length(getDims(obj)), 7)
  expect_equal(getDims(obj)[["raw"]], dim(getRawData(obj)))
  expect_equal(getDims(obj)[["metaCells"]], ncol(getMetadataCells(obj)))
  expect_equal(getDims(obj)[["clusterCoex"]], NULL)
})
