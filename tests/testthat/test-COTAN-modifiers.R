
test_that("metaDataset", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  tags <- c("GEO:", "scRNAseq method:", "starting n. of cells:", "Condition sample:")

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X", sampleCondition = "Test")

  meta <- getMetadataDataset(obj)

  expect_equal(meta[[1]], tags)
  expect_equal(meta[[2]], c("V", "10X", "20", "Test"))

  obj <- addElementToMetaDataset(obj, tag = "Tag_1", value = 1)
  obj <- addElementToMetaDataset(obj, tag = "Tag_2", value = "Test")
  obj <- addElementToMetaDataset(obj, tag = "Tag_3", value = c("Array", "of", "strings"))

  meta <- getMetadataDataset(obj)

  expect_equal(meta[nrow(meta) - 2, 2], as.character(1))
  expect_equal(meta[nrow(meta) - 1, 2], "Test")
  expect_equal(meta[nrow(meta), ], list("Tag_3", "Array", "of", "strings"), ignore_attr = TRUE)
})
