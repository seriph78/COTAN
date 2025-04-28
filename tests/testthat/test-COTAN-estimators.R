
test_that("Linear estimates", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)

  obj <- estimateLambdaLinear(obj)

  expect_identical(getLambda(obj), rowMeans(getRawData(obj), dims = 1L))
  expect_equal(getMetadataDataset(obj)[[1L]], datasetTags()[5L:6L],
               ignore_attr = TRUE)
  expect_identical(getMetadataDataset(obj)[[2L]], c("FALSE", "FALSE"))

  obj <- estimateNuLinear(obj)

  expect_identical(getNu(obj),
                   colMeans(getRawData(obj), dims = 1L) /
                     mean(colMeans(getRawData(obj), dims = 1L)))

  clusters <- set_names(rep(1L:2L, times = 10L), getCells(obj))
  obj <- estimateNuLinearByCluster(obj, clusters = clusters)

  expect_identical(getNu(obj), set_names(rep_len(1.0, getNumCells(obj)),
                                         getCells(obj)))
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
