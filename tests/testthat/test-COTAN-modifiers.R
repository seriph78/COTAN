
test_that("metaDataset", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")

  genesCoexInSync <- getMetadataElement(obj, datasetTags()[["gsync"]])
  cellsCoexInSync <- getMetadataElement(obj, datasetTags()[["csync"]])

  expect_identical(c(genesCoexInSync, genesCoexInSync), c("FALSE", "FALSE"))

  meta <- getMetadataDataset(obj)

  expect_equal(meta[[1L]], head(datasetTags(), 6L), ignore_attr = TRUE)
  expect_identical(meta[[2L]], c("V", "10X", "Test", "20",
                                 genesCoexInSync, cellsCoexInSync))

  obj <- addElementToMetaDataset(obj, tag = "Tag_1", value = 1L)
  obj <- addElementToMetaDataset(obj, tag = "Tag_2", value = "Test")
  obj <- addElementToMetaDataset(obj, tag = "Tag_3",
                                 value = c("Array", "of", "strings"))

  meta <- getMetadataDataset(obj)

  expect_identical(meta[nrow(meta) - 2L, 2L], "1")
  expect_identical(meta[nrow(meta) - 1L, 2L], "Test")
  expect_equal(getMetadataElement(obj, "Tag_3"),
               list("Array", "of", "strings"), ignore_attr = TRUE)
})


test_that("Fully-expressed Genes and Fully-expressing Cells", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)

  obj <- findFullyExpressedGenes(obj)

  expect_identical(flagNotFullyExpressedGenes(obj),
                   c(FALSE, rep(TRUE, getNumGenes(obj) - 1L)))
  expect_identical(getFullyExpressedGenes(obj), LETTERS[[1L]])

  obj <- findFullyExpressingCells(obj)

  expect_identical(flagNotFullyExpressingCells(obj),
                   rep(TRUE, getNumCells(obj)))
  expect_identical(getFullyExpressingCells(obj), vector(mode = "character"))
})


test_that("dropGenesCells", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  tags <- c("GEO:", "scRNAseq method:",
            "starting n. of cells:", "Condition sample:")

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")
  obj <- clean(obj)

  obj <- estimateLambdaLinear(obj)
  obj <- estimateDispersionNuBisection(obj, cores = 4L,
                                       enforceNuAverageToOne = TRUE)
  suppressWarnings({
    obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = TRUE)
  })
  expect_no_warning({
    obj <- calculateCoex(obj, actOnCells = TRUE, optimizeForSpeed = FALSE)
  })

  gdiDF <- calculateGDI(obj)
  obj <- storeGDI(obj, getColumnFromDF(gdiDF, "GDI"))
  expect_equal(gdiDF[, "GDI"], getGDI(obj), ignore_attr = TRUE)

  obj <- dropGenesCoex(obj)

  expect_error(getGenesCoex(obj))
  genesSel <- c("F", "D", "A", "C")
  expect_identical(dim(suppressWarnings(getGenesCoex(obj, genes = genesSel))),
                   c(getNumGenes(obj), 4L))
  expect_false(is_empty(getCellsCoex(obj)))

  obj <- dropCellsCoex(obj)

  expect_error(getCellsCoex(obj))
  cellsSel <- c("f", "d", "a", "c")
  expect_identical(dim(suppressWarnings(getCellsCoex(obj, cells = cellsSel))),
                   c(getNumCells(obj), 4L))

  obj <- dropGenesCells(obj, genes = LETTERS[8L:9L])

  expect_identical(getNumGenes(obj), 8L)
  expect_false(any(getGenes(obj) %in% LETTERS[8L:9L]))
  expect_identical(getCells(obj), letters[1L:20L])

  obj <- dropGenesCells(obj, cells = c(letters[8L:9L], letters[18L:19L]))

  expect_identical(getNumCells(obj), 16L)
  expect_false(any(getCells(obj) %in% c(letters[8L:9L], letters[18L:19L])))
  expect_false(any(getGenes(obj) %in% LETTERS[8L:9L]))

  obj <- dropGenesCells(obj, genes = LETTERS[6L:7L], cells = letters[6L:7L])

  expect_identical(getNumGenes(obj), 6L)
  expect_false(any(getGenes(obj) %in% LETTERS[6L:9L]))
  expect_identical(getNumCells(obj), 14L)
  expect_false(any(getCells(obj) %in% c(letters[6L:9L], letters[18L:19L])))
})


test_that("Managed clusterizations and conditions", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)

  clusters <- factor(set_names(rep(1L:2L, 10L), getCells(obj)))
  obj <- addClusterization(obj, clName = "Test",
                           clusters = clusters)

  expect_identical(getClusterizations(obj), "Test")
  expect_length(getClusterizations(obj, dropNoCoex = TRUE), 0L)
  expect_identical(colnames(getMetadataCells(obj)), "CL_Test")
  expect_identical(rownames(getMetadataCells(obj)), getCells(obj))
  expect_identical(getClusters(obj), clusters)

  obj <- estimateNuLinear(obj)

  expect_identical(colnames(getMetadataCells(obj)), c("CL_Test", "nu"))

  coexDF <-
    set_names(as.data.frame(atan(getNuNormData(obj)[, 1L:2L] - 0.5) / pi * 2.0),
              1L:2L)

  obj <- addClusterizationCoex(obj, clName = "Test", coexDF)

  expect_identical(getClusterizations(obj, dropNoCoex = TRUE,
                                      keepPrefix = TRUE), "CL_Test")
  expect_identical(getClusterizationData(obj)[["coex"]], coexDF)

  clusters2 <- factor(set_names(rep(c("2", "1"), 10L), getCells(obj)))
  coexDF2 <-
    set_names(as.data.frame(atan(getNuNormData(obj)[, 1L:2L] - 0.7) / pi * 2.0),
              c("1", "2"))

  obj <- addClusterization(obj, clName = "Test2",
                           clusters = clusters2,
                           coexDF = coexDF2)

  expect_named(getClustersCoex(obj), c("CL_Test", "CL_Test2"))
  expect_identical(getClusterizations(obj), c("Test", "Test2"))
  expect_identical(colnames(getMetadataCells(obj)),
                   c("CL_Test", "nu", "CL_Test2"))
  expect_identical(getClusterizationData(obj),
                   list("clusters" = clusters2, "coex" = coexDF2))
  expect_identical(getClusterizationData(obj, clName = "Test"),
                   list("clusters" = clusters, "coex" = coexDF))

  obj <- dropClusterization(obj, "Test")

  expect_named(getClustersCoex(obj), "CL_Test2")
  expect_identical(getClusterizations(obj), "Test2")
  expect_identical(colnames(getMetadataCells(obj)), c("nu", "CL_Test2"))
  expect_identical(getClusterizationData(obj),
                   list("clusters" = clusters2, "coex" = coexDF2))

  genre <- factor(set_names(rep(c("F", "M"), 10L), getCells(obj)))

  obj <- addCondition(obj, condName = "Test", conditions = genre)

  expect_identical(getAllConditions(obj), "Test")
  expect_identical(colnames(getMetadataCells(obj)),
                   c("nu", "CL_Test2", "COND_Test"))
  expect_identical(levels(getCondition(obj)), "NoCond")
  expect_identical(getCondition(obj, condName = "Test"), genre)

  obj <- dropCondition(obj, condName = "Test")

  expect_identical(getAllConditions(obj), vector(mode = "character"))
  expect_identical(colnames(getMetadataCells(obj)), c("nu", "CL_Test2"))

  # no such clusterization/condition
  expect_error(getClusterizationData(obj, clName = "Test"))
  expect_error(getClusters(obj, clName = "Test"))
  expect_error(getCondition(obj, condName = "Test"))

  # empty name clusterization/condition
  expect_error(addClusterization(obj, clName = "", clusters = rep(0L, 20L)))
  expect_error(addCondition(obj, condName = "", conditions = rep(0L, 20L)))

  # no names in clusterization
  expect_error(addClusterization(obj, clName = "Test", clusters = rep(0L, 20L)))
  expect_error(addCondition(obj, condName = "Cond", conditions = rep(0L, 20L)))

  # already existing clusterization
  expect_error(addClusterization(obj, clName = "Test2",
                                 clusters = clusters2))

  # wrong clusters/conditions size
  expect_error(addClusterization(obj, clName = "Test", clusters = rep(0L, 17L)))
  expect_error(addCondition(obj, condName = "Test", conditions = rep("A", 17L)))

  # wrong COEX data.frame size
  expect_error(addClusterization(obj, clName = "Test",
                                 clusters = clusters,
                                 coexDF = coexDF[1L:8L, ]))

  # wrong COEX data.frame column names
  expect_error(addClusterization(obj, clName = "Test",
                                 clusters = rep(c("3", "4"), 20L),
                                 coexDF = coexDF))
})
