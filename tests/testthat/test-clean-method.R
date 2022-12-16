library(gsubfn)

test_that("clean COTAN object", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)
  obj <- clean(obj)

  expect_true(validObject(obj))

  expect_equal(colnames(getMetadataCells(obj)), c("nu", "feCells"))
  expect_equal(colnames(getMetadataGenes(obj)), c("lambda","hkGenes"))
})

test_that("Clean on test dataset", {
  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- clean(obj)

  obj <- estimateDispersionBisection(obj, cores = 12)

  raw.norm <- readRDS(file.path(getwd(), "raw.norm.test.RDS"))
  lambda <- readRDS(file.path(getwd(), "lambda.test.RDS"))
  dispersion <- readRDS(file.path(getwd(), "a.test.RDS"))
  nu <- readRDS(file.path(getwd(), "nu.test.RDS"))

  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))
  cells.names.test <- readRDS(file.path(getwd(), "cell.names.test.RDS"))

  expect_equal(getNormalizedData(obj)[genes.names.test, cells.names.test], raw.norm)
  expect_equal(getLambda(obj)[genes.names.test], lambda)
  expect_equal(getNu(obj)[cells.names.test] , nu)
  expect_equal(getDispersion(obj)[genes.names.test], dispersion)
})
