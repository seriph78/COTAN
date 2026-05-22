
prevOptState <- options(parallelly.fork.enable = TRUE)

test_that("clean COTAN object", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- clean(obj)

  expect_true(validObject(obj))

  expect_identical(colnames(getMetadataCells(obj)), c("feCells", "nu"))
  expect_identical(colnames(getMetadataGenes(obj)), "feGenes")
})

test_that("Clean on test dataset", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- clean(obj)

  obj <- estimateLambdaLinear(obj)
  obj <- estimateDispersionViaSolver(obj, cores = 6L)

  rawNorm <- readRDS(file.path(getwd(), "raw.norm.test.RDS"))
  lambda <- readRDS(file.path(getwd(), "lambda.test.RDS"))
  dispersion <- readRDS(file.path(getwd(), "dispersion.test.RDS"))
  nu <- readRDS(file.path(getwd(), "nu.test.RDS"))

  genesNamesTest <- readRDS(file.path(getwd(), "genes.names.test.RDS"))
  cellsNamesTest <- readRDS(file.path(getwd(), "cells.names.test.RDS"))

  expect_equal(getNuNormData(obj)[genesNamesTest, cellsNamesTest],
               rawNorm, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getLambda(obj)[genesNamesTest],
               lambda, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getNu(obj)[cellsNamesTest],
               nu, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getDispersion(obj)[genesNamesTest],
               dispersion, tolerance = 1.0e-10, ignore_attr = FALSE)
})


test_that("clean accepts CleaningOptions and rejects mixed legacy arguments", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)

  expect_true(validObject(clean(obj, cleaningOptions = CleaningOptions())))

  expect_error(
    clean(
      obj,
      cellsCutoff = 0.005,
      cleaningOptions = CleaningOptions()
    ),
    regexp = "Do not mix `cleaningOptions`"
  )

  expect_error(
    clean(
      obj,
      cleaningOptions = ReductionOptions()
    ),
    regexp = "unable to find an inherited method"
  )
})


test_that("proceedToCoex accepts independent execution and cleaning option packs", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  exec <- ExecutionOptions(cores = 1L, optimizeForSpeed = FALSE)
  cleanOpt <- CleaningOptions()

  expect_true(validObject(proceedToCoex(
    obj,
    calcCoex = FALSE,
    cleaningOptions = cleanOpt,
    saveObj = FALSE
  )))

  expect_true(validObject(proceedToCoex(
    obj,
    calcCoex = FALSE,
    executionOptions = exec,
    saveObj = FALSE
  )))

  expect_true(validObject(proceedToCoex(
    obj,
    calcCoex = FALSE,
    executionOptions = exec,
    cleaningOptions = cleanOpt,
    saveObj = FALSE
  )))

  expect_error(
    proceedToCoex(
      obj,
      calcCoex = FALSE,
      cellsCutoff = 0.005,
      cleaningOptions = cleanOpt,
      saveObj = FALSE
    ),
    regexp = "Do not mix `cleaningOptions`"
  )

  expect_error(
    proceedToCoex(
      obj,
      calcCoex = FALSE,
      cores = 2L,
      executionOptions = exec,
      saveObj = FALSE
    ),
    regexp = "Do not mix `executionOptions`"
  )
})

test_that("proceedToCoex accepts execution and cleaning option packs", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  exec <- ExecutionOptions(cores = 1L, optimizeForSpeed = FALSE)
  cleanOpt <- CleaningOptions()

  expect_true(validObject(proceedToCoex(
    obj,
    calcCoex = FALSE,
    executionOptions = exec,
    saveObj = FALSE
  )))

  expect_true(validObject(proceedToCoex(
    obj,
    calcCoex = FALSE,
    cleaningOptions = cleanOpt,
    saveObj = FALSE
  )))

  expect_true(validObject(proceedToCoex(
    obj,
    calcCoex = FALSE,
    executionOptions = exec,
    cleaningOptions = cleanOpt,
    saveObj = FALSE
  )))

  expect_error(
    proceedToCoex(
      obj,
      calcCoex = FALSE,
      cellsCutoff = 0.005,
      cleaningOptions = cleanOpt,
      saveObj = FALSE
    ),
    regexp = "Do not mix `cleaningOptions`"
  )

  expect_error(
    proceedToCoex(
      obj,
      calcCoex = FALSE,
      cores = 2L,
      executionOptions = exec,
      saveObj = FALSE
    ),
    regexp = "Do not mix `executionOptions`"
  )
})

options(prevOptState)
