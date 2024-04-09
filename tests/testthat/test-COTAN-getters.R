
test_that("COTAN getters", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")
  obj <- clean(obj)
  obj <- estimateDispersionBisection(obj)
  obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = FALSE)
  obj <- calculateCoex(obj, actOnCells = TRUE,  optimizeForSpeed = TRUE)

  obj <- addClusterization(obj, clName = "Test",
                           clusters = set_names(rep(c(1L, 2L), 10L),
                                                colnames(raw)))
  obj <- addClusterization(obj, clName = "Test2",
                           clusters = set_names(rep(c(2L, 1L), 10L),
                                                colnames(raw)))

  obj <- addCondition(obj, condName = "Test",
                      conditions = set_names(rep(c("F", "M"), 10L),
                                             colnames(raw)))
  obj <- addCondition(obj, condName = "Test", override = TRUE,
                      conditions = set_names(c(rep(c("F", "M"), 9L), "M", "F"),
                                             colnames(raw)))

  metaInfo <- c("V", "10X", "Test", "20", "TRUE", "TRUE",
                paste0(10.0 / 55.0), paste0(0L))

  expect_identical(getRawData(obj), as(as(raw, "dMatrix"), "sparseMatrix"))
  expect_identical(getNumGenes(obj), 10L)
  expect_identical(getNumCells(obj), 20L)
  expect_identical(getGenes(obj), head(LETTERS, getNumGenes(obj)))
  expect_identical(getCells(obj), head(letters, getNumCells(obj)))
  expect_identical(getZeroOneProj(obj), sign(getRawData(obj)))
  expect_identical(getCellsSize(obj), colSums(getRawData(obj)))
  expect_identical(getNumExpressedGenes(obj), colSums(getRawData(obj) != 0L))
  expect_identical(getGenesSize(obj), rowSums(getRawData(obj)))
  expect_identical(getNumOfExpressingCells(obj), rowSums(getRawData(obj) != 0L))
  expect_identical(getNormalizedData(obj),
                   t(t(getRawData(obj)) * (1.0 / getNu(obj))))
  expect_equal(getMetadataDataset(obj)[[1L]], datasetTags()[1L:8L],
               ignore_attr = TRUE)
  expect_identical(getMetadataDataset(obj)[[2L]], metaInfo)
  expect_setequal(colnames(getMetadataGenes(obj)),
                  c("lambda", "feGenes", "dispersion"))
  expect_identical(rownames(getMetadataGenes(obj)), getGenes(obj))
  expect_setequal(colnames(getMetadataCells(obj)),
                  c("nu", "feCells", names(getClustersCoex(obj)),
                    getAllConditions(obj, keepPrefix = TRUE)))
  expect_identical(rownames(getMetadataCells(obj)), getCells(obj))
  expect_length(getClustersCoex(obj), 2L)
  expect_named(getClustersCoex(obj), paste0("CL_", getClusterizations(obj)))
  expect_identical(getClustersCoex(obj)[["CL_Test"]], data.frame())
  expect_equal(getNu(obj), getMetadataCells(obj)[["nu"]],
               ignore_attr = TRUE)
  expect_equal(getLambda(obj), getMetadataGenes(obj)[["lambda"]],
               ignore_attr = TRUE)
  expect_equal(getDispersion(obj), getMetadataGenes(obj)[["dispersion"]],
               ignore_attr = TRUE)
  expect_identical(flagNotFullyExpressedGenes(obj),
                   c(FALSE, rep(TRUE, getNumGenes(obj) - 1L)))
  expect_identical(getFullyExpressedGenes(obj), LETTERS[[1L]])
  expect_identical(flagNotFullyExpressingCells(obj),
                   rep(TRUE, getNumCells(obj)))
  expect_identical(getFullyExpressingCells(obj), vector(mode = "character"))
  expect_identical(dim(getGenesCoex(obj)),
                   as.integer(c(getNumGenes(obj), getNumGenes(obj))))
  expect_identical(dim(getCellsCoex(obj)),
                   as.integer(c(getNumCells(obj), getNumCells(obj))))
  expect_identical(getClusterizations(obj), c("Test", "Test2"))
  expect_identical(getClusterizationName(obj), "Test2")
  expect_identical(getClusterizationName(obj, clName = "Test",
                                         keepPrefix = TRUE), "CL_Test")
  expect_setequal(names(getClusterizationData(obj)), c("coex", "clusters"))
  expect_equal(getClusterizationData(obj)[["clusters"]],
               getMetadataCells(obj)[["CL_Test2"]], ignore_attr = TRUE)
  expect_identical(getClusters(obj), getClusterizationData(obj)[["clusters"]])
  expect_identical(getClusterizationData(obj)[["coex"]], data.frame())
  expect_identical(getAllConditions(obj), "Test")
  expect_identical(getConditionName(obj), "")
  expect_identical(getConditionName(obj, condName = "Test",
                                    keepPrefix = TRUE), "COND_Test")
  expect_identical(levels(getCondition(obj)), "NoCond")
  expect_identical(factorToVector(getCondition(obj,
                                               condName = "Test"))[19L:20L],
                   set_names(c("M", "F"), letters[19L:20L]))
  expect_length(getDims(obj), 7L)
  expect_identical(getDims(obj)[["raw"]], dim(getRawData(obj)))
  expect_identical(getDims(obj)[["metaCells"]], ncol(getMetadataCells(obj)))
  expect_null(getDims(obj)[["clusterCoex"]])
})
