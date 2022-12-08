
test_that("Calculations on genes", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)
  obj <- clean(obj, calcExtraData = FALSE)[[1]]

  mu <- calculateMu(obj)

  expect_equal(dim(mu), dim(getRawData(obj)))
  expect_equal(mu[ 1,  1], getLambda(obj)[ 1] * getNu(obj)[ 1], ignore_attr = TRUE)
  expect_equal(mu[10,  1], getLambda(obj)[10] * getNu(obj)[ 1], ignore_attr = TRUE)
  expect_equal(mu[ 1, 20], getLambda(obj)[ 1] * getNu(obj)[20], ignore_attr = TRUE)
  expect_equal(mu[10, 10], getLambda(obj)[10] * getNu(obj)[10], ignore_attr = TRUE)

  list[observedYY, observedY] <-
    observedContingencyTablesYY(obj, actOnCells = FALSE, asDspMatrices = FALSE)

  expect_s4_class(observedYY, "dsyMatrix")
  expect_equal(dim(observedYY), rep(getNumGenes(obj), 2))
  expect_equal(diag(as.matrix(observedYY)), observedY)
  expect_equal(length(observedY), getNumGenes(obj))
  expect_equal(observedY, c(20, rep(10 ,9)), ignore_attr = TRUE)

  list[observedNN, observedNY, observedYN, ] <-
    observedContingencyTables(obj, actOnCells = FALSE, asDspMatrices = FALSE)

  expect_equal(unlist(lapply(c(observedNN, observedNY, observedYN), dim)), rep(dim(observedYY), 3))
  expect_equal(diag(as.matrix(observedNN)), c(0, rep(10, 9)), ignore_attr = TRUE)
  expect_equal(observedNY, t(observedYN))
  expect_equal(diag(as.matrix(observedNY)), rep(0, 10), ignore_attr = TRUE)
  expect_equal(as.matrix(observedNN + observedNY + observedYN + observedYY),
               matrix(getNumCells(obj), nrow = getNumGenes(obj), ncol = getNumGenes(obj)),
               ignore_attr = TRUE)

  obj <- estimateDispersionBisection(obj, step = 4, cores = 4)

  list[expectedNN, expectedN] <-
    expectedContingencyTablesNN(obj, actOnCells = FALSE, asDspMatrices = TRUE)

  expect_s4_class(expectedNN, "dspMatrix")
  expect_equal(dim(expectedNN), rep(getNumGenes(obj), 2))
  expect_equal(length(expectedN), getNumGenes(obj))
  expect_equal(expectedN, c(0, rep(10, 9)), ignore_attr = TRUE, tolerance = 1e-4)

  list[, expectedNY, expectedYN, expectedYY] <-
    expectedContingencyTables(obj, actOnCells = FALSE, asDspMatrices = TRUE)

  expect_equal(unlist(lapply(c(expectedYY, expectedNY, expectedYN), dim)), rep(dim(expectedNN), 3))
  expect_equal(as.matrix(expectedNN + expectedNY + expectedYN + expectedYY),
               matrix(getNumCells(obj), nrow = getNumGenes(obj), ncol = getNumGenes(obj)),
               ignore_attr = TRUE)
})


test_that("Calculations on cells", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)
  obj <- clean(obj, calcExtraData = FALSE)[[1]]

  list[observedYY, observedY] <-
    observedContingencyTablesYY(obj, actOnCells = TRUE, asDspMatrices = FALSE)

  expect_s4_class(observedYY, "dsyMatrix")
  expect_equal(dim(observedYY), rep(getNumCells(obj), 2))
  expect_equal(diag(as.matrix(observedYY)), observedY)
  expect_equal(length(observedY), getNumCells(obj))
  expect_equal(observedY, rep(c(7,4), 10), ignore_attr = TRUE)

  list[observedNN, observedNY, observedYN, ] <-
    observedContingencyTables(obj, actOnCells = TRUE, asDspMatrices = FALSE)

  expect_equal(unlist(lapply(c(observedNN, observedNY, observedYN), dim)), rep(dim(observedYY), 3))
  expect_equal(diag(as.matrix(observedNN)), rep(c(3, 6), 10), ignore_attr = TRUE)
  expect_equal(observedNY, t(observedYN))
  expect_equal(diag(as.matrix(observedNY)), rep(0, 20), ignore_attr = TRUE)
  expect_equal(as.matrix(observedNN + observedNY + observedYN + observedYY),
               matrix(getNumGenes(obj), nrow = getNumCells(obj), ncol = getNumCells(obj)),
               ignore_attr = TRUE)

  obj <- estimateDispersionBisection(obj, step = 4, cores = 4)
  obj <- estimateNuBisection(obj, step = 4, cores = 4)

  list[expectedNN, expectedN] <-
    expectedContingencyTablesNN(obj, actOnCells = TRUE, asDspMatrices = TRUE)

  expect_s4_class(expectedNN, "dspMatrix")
  expect_equal(dim(expectedNN), rep(getNumCells(obj), 2))
  expect_equal(length(expectedN), getNumCells(obj))
  expect_equal(expectedN, rep(c(3, 6), 10), ignore_attr = TRUE, tolerance = 1e-3)

  list[, expectedNY, expectedYN, expectedYY] <-
    expectedContingencyTables(obj, actOnCells = TRUE, asDspMatrices = TRUE)

  expect_equal(unlist(lapply(c(expectedYY, expectedNY, expectedYN), dim)), rep(dim(expectedNN), 3))
  expect_equal(as.matrix(expectedNN + expectedNY + expectedYN + expectedYY),
               matrix(getNumGenes(obj), nrow = getNumCells(obj), ncol = getNumCells(obj)),
               ignore_attr = TRUE)
})


test_that("Coex", {
  skip("Test not ready yet")

  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  obj <- COTAN(raw = raw)
  obj <- clean(obj, calcExtraData = FALSE)[[1]]

  obj <- estimateDispersionBisection(obj, step = 4, cores = 4)
  obj <- estimateNuBisection(obj, step = 4, cores = 4)

  obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = FALSE)

  expect_equal(dim(getGenesCoex(obj)), rep(getNumGenes(obj), 2))
  #expect_equal

  obj <- calculateCoex(obj, actOnCells = TRUE, optimizeForSpeed = TRUE)

  expect_equal(dim(getCellsCoex(obj)), rep(getNumCells(obj), 2))
  #expect_equal

  calculateS(obj)
  calculateG(obj)
  calculatePValue(obj, statType = "S")
  calculatePValue(obj, statType = "G",
                  geneSubsetCol = getGenes(obj)[5:15],
                  geneSubsetRow = getGenes(obj)[1:10])
  calculateGDI(obj,statType = "G")
  calculateGDI(obj,statType = "S")
})


test_that("Coex vs saved results", {
  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- clean(obj, calcExtraData = FALSE)[["objCOTAN"]]

  obj <- estimateDispersionBisection(obj, cores = 12)

  obj <- calculateCoex(obj)

  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))

  coex_test <- readRDS(file.path(getwd(), "coex.test.RDS"))

  expect_equal(as.matrix(getGenesCoex(obj, genes = genes.names.test)), coex_test)

  pval <- calculatePValue(obj, geneSubsetCol = genes.names.test,
                          geneSubsetRow = genes.names.test)

  pval_exp  <- readRDS(file.path(getwd(), "pval.test.RDS"))

  expect_equal(pval, pval_exp)

  GDI <- calculateGDI(obj)[genes.names.test, ]

  GDI_exp <- readRDS(file.path(getwd(), "GDI.test.RDS"))
  expect_equal(GDI, GDI_exp)
})
