
test_that("Linear estimates", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- addCondition(obj, condName = "batch",
                      conditions = rlang::rep_named(getCells(obj), c(1L, 2L)))

  obj <- estimateLambdaLinear(obj)

  expect_identical(getLambda(obj), rowMeans(getRawData(obj), dims = 1L))
  expect_equal(getMetadataDataset(obj)[[1L]], datasetTags()[5L:6L],
               ignore_attr = TRUE)
  expect_identical(getMetadataDataset(obj)[[2L]], c("FALSE", "FALSE"))

  obj <- estimateNuLinear(obj)

  expect_identical(getNu(obj),
                   Matrix::colMeans(getRawData(obj)) /
                     mean(Matrix::colMeans(getRawData(obj))))

  obj <- resetBatches(obj, conditionToUse = "batch")

  obj <- proceedToCoex(obj, calcCoex = FALSE, cores = 6L, saveObj = FALSE)

  expect_identical(getNu(obj), rlang::rep_named(getCells(obj), 1.0))

  expect_identical(getLambda(obj, batchName = "All"),
                   as.matrix(getRawData(obj)))
  expect_identical(getLambda(obj, batchName = "1"), getRawData(obj)[, 1L])
  expect_identical(getLambda(obj, batchName = "2"), getRawData(obj)[, 2L])
})


test_that("Bisection estimates", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- clean(obj)

  obj <- estimateLambdaLinear(obj)
  obj <- estimateDispersionBisection(obj, cores = 3L, chunkSize = 2L)

  expect_length(getDispersion(obj), getNumGenes(obj))
  expect_equal(getDispersion(obj)[[1L]], -Inf, ignore_attr = TRUE)

  expect_equal(rowSums(getZeroOneProj(obj) + getProbabilityOfZero(obj)),
               rep(getNumCells(obj), getNumGenes(obj)),
               tolerance = 0.001, ignore_attr = TRUE)

  obj <- estimateNuBisection(obj, cores = 6L, chunkSize = 3L)

  expect_length(getNu(obj), getNumCells(obj))

  expect_equal(colSums(getZeroOneProj(obj) + getProbabilityOfZero(obj)),
               rep(getNumGenes(obj), getNumCells(obj)),
               tolerance = 0.001, ignore_attr = TRUE)

  obj <- estimateDispersionNuBisection(obj, enforceNuAverageToOne = TRUE)

  expect_length(getDispersion(obj), getNumGenes(obj))
  expect_equal(getDispersion(obj)[[1L]], -Inf, ignore_attr = TRUE)
  expect_length(getNu(obj), getNumCells(obj))
  expect_equal(mean(getNu(obj)), 1.0, tolerance = 1.0e-12)

  expect_equal(rowSums(getZeroOneProj(obj) + getProbabilityOfZero(obj)),
               rep(getNumCells(obj), getNumGenes(obj)),
               tolerance = 0.001, ignore_attr = TRUE)

  expect_equal(colSums(getZeroOneProj(obj) + getProbabilityOfZero(obj)),
               rep(getNumGenes(obj), getNumCells(obj)),
               tolerance = 0.001, ignore_attr = TRUE)

  expect_error(estimateDispersionNuNlminb(obj))
})
