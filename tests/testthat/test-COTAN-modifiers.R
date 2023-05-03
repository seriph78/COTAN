
test_that("metaDataset", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0),
                nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")

  genesCoexInSync <- getMetadataElement(obj, datasetTags()[["gsync"]])
  cellsCoexInSync <- getMetadataElement(obj, datasetTags()[["csync"]])

  expect_equal(c(genesCoexInSync, genesCoexInSync), c("FALSE", "FALSE"))

  meta <- getMetadataDataset(obj)

  expect_equal(meta[[1]], head(datasetTags(), 6), ignore_attr = TRUE)
  expect_equal(meta[[2]], c("V", "10X", "Test", "20",
                            genesCoexInSync, cellsCoexInSync))

  obj <- addElementToMetaDataset(obj, tag = "Tag_1", value = 1)
  obj <- addElementToMetaDataset(obj, tag = "Tag_2", value = "Test")
  obj <- addElementToMetaDataset(obj, tag = "Tag_3",
                                 value = c("Array", "of", "strings"))

  meta <- getMetadataDataset(obj)

  expect_equal(meta[nrow(meta) - 2, 2], as.character(1))
  expect_equal(meta[nrow(meta) - 1, 2], "Test")
  expect_equal(getMetadataElement(obj, "Tag_3"),
               list("Array", "of", "strings"), ignore_attr = TRUE)
})


test_that("Fully-expressed Genes and Fully-expressing Cells", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0),
                nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)

  obj <- findFullyExpressedGenes(obj)

  expect_equal(flagNotFullyExpressedGenes(obj),
               c(FALSE, rep(TRUE, getNumGenes(obj) - 1)))
  expect_equal(getFullyExpressedGenes(obj), c(LETTERS[1]))

  obj <- findFullyExpressingCells(obj)

  expect_equal(flagNotFullyExpressingCells(obj), rep(TRUE, getNumCells(obj)))
  expect_equal(getFullyExpressingCells(obj), vector(mode = "character"))
})


test_that("dropGenesCells", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0),
                nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  tags <- c("GEO:", "scRNAseq method:",
            "starting n. of cells:", "Condition sample:")

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")
  obj <- clean(obj)
  obj <- estimateDispersionNuBisection(obj, cores = 4,
                                       enforceNuAverageToOne = TRUE)
  obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = TRUE)
  obj <- calculateCoex(obj, actOnCells = TRUE,  optimizeForSpeed = FALSE)

  obj <- dropGenesCoex(obj)

  expect_error(getGenesCoex(obj))
  expect_false(is_empty(getCellsCoex(obj)))

  obj <- dropCellsCoex(obj)

  expect_error(getCellsCoex(obj))

  obj <- dropGenesCells(obj, genes = LETTERS[8:9])

  expect_equal(getNumGenes(obj), 8)
  expect_true(all(!getGenes(obj) %in% LETTERS[8:9]))
  expect_equal(getCells(obj), letters[1:20])

  obj <- dropGenesCells(obj, cells = c(letters[8:9], letters[18:19]))

  expect_equal(getNumCells(obj), 16)
  expect_true(all(!getCells(obj) %in% c(letters[8:9], letters[18:19])))
  expect_true(all(!getGenes(obj) %in% LETTERS[8:9]))

  obj <- dropGenesCells(obj, genes = LETTERS[6:7], cells = letters[6:7])

  expect_equal(getNumGenes(obj), 6)
  expect_true(all(!getGenes(obj) %in% LETTERS[6:9]))
  expect_equal(getNumCells(obj), 14)
  expect_true(all(!getCells(obj) %in% c(letters[6:9], letters[18:19])))
})


test_that("Managed clusterizations", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0),
                nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)

  clusters <- set_names(rep(c(1, 2), 10), getCells(obj))
  obj <- addClusterization(obj, clName = "Test",
                           clusters = clusters)

  expect_equal(getClusterizations(obj), "Test")
  expect_length(getClusterizations(obj, dropNoCoex = TRUE), 0)
  expect_equal(colnames(getMetadataCells(obj)), "CL_Test")
  expect_equal(rownames(getMetadataCells(obj)), getCells(obj))
  expect_equal(getClusterizationData(obj)[["clusters"]], clusters)

  obj <- estimateNuLinear(obj)

  expect_equal(colnames(getMetadataCells(obj)), c("CL_Test", "nu"))

  coexDF <- set_names(
    as.data.frame(atan(getNormalizedData(obj)[,1:2]-0.5)/pi*2), c(1, 2))

  obj <- addClusterizationCoex(obj, clName = "Test", coexDF)

  expect_equal(getClusterizations(obj, dropNoCoex = TRUE, keepPrefix = TRUE),
               "CL_Test")
  expect_equal(getClusterizationData(obj)[["coex"]], coexDF)

  clusters2 = set_names(rep(c("2", "1"), 10), getCells(obj))
  coexDF2 <- set_names(
    as.data.frame(atan(getNormalizedData(obj)[,1:2]-0.7)/pi*2), c("1", "2"))

  obj <- addClusterization(obj, clName = "Test2",
                           clusters = clusters2,
                           coexDF = coexDF2)

  expect_equal(names(getClustersCoex(obj)), c("CL_Test", "CL_Test2"))
  expect_equal(getClusterizations(obj), c("Test", "Test2"))
  expect_equal(colnames(getMetadataCells(obj)), c("CL_Test", "nu", "CL_Test2"))
  expect_equal(getClusterizationData(obj),
               list("clusters" = clusters2, "coex" = coexDF2))
  expect_equal(getClusterizationData(obj, clName = "Test"),
               list("clusters" = clusters, "coex" = coexDF))

  obj <- dropClusterization(obj, "Test")

  expect_equal(names(getClustersCoex(obj)), "CL_Test2")
  expect_equal(getClusterizations(obj), "Test2")
  expect_equal(colnames(getMetadataCells(obj)), c("nu", "CL_Test2"))
  expect_equal(getClusterizationData(obj),
               list("clusters" = clusters2, "coex" = coexDF2))

  # no such clusterization
  expect_error(getClusterizationData(obj, clName = "Test"))

  # empyt name clusterization
  expect_error(addClusterization(obj, clName = "", clusters = rep(0, 20)))

  # already existing clusterization
  expect_error(addClusterization(obj, clName = "Test2", clusters = rep(0, 20)))

  # wrong clusters size
  expect_error(addClusterization(obj, clName = "Test", clusters = rep(0, 17)))

  # wrong coex data.frame size
  expect_error(addClusterization(obj, clName = "Test",
                                 clusters = clusters,
                                 coexDF = coexDF[1:8, ]))

  # wrong coex data.frame column names
  expect_error(addClusterization(obj, clName = "Test",
                                 clusters = rep(c("3", "4"), 20),
                                 coexDF = coexDF))
})
