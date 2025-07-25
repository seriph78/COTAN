tm <- tempdir()
stopifnot(file.exists(tm))

library(zeallot)
library(rlang)

options(COTAN.TorchWarning = NULL)
options(parallelly.fork.enable = TRUE)

crossEntrVector <- function(zeroOne, probZero) {
  crossEntr <- rep_len(0.0, nrow(zeroOne))
  for (r in seq_len(nrow(zeroOne))) {
    for (c in seq_len(ncol(zeroOne))) {
      if (probZero[r, c] != 0.0) {
        crossEntr[r] <- crossEntr[r] -
          (1.0 - zeroOne[r, c]) * log(probZero[r, c]) -
          zeroOne[r, c] * log(1.0 - probZero[r, c])
      }
    }
    crossEntr[r] <- crossEntr[r] / ncol(zeroOne)
  }

  return(crossEntr)
}

coexPoint <- function(o, e, n) {
  num <- ( ((o[[1L]] - e[[1L]]) / max(1.0, e[[1L]])) -
           ((o[[2L]] - e[[2L]]) / max(1.0, e[[2L]])) -
           ((o[[3L]] - e[[3L]]) / max(1.0, e[[3L]])) +
           ((o[[4L]] - e[[4L]]) / max(1.0, e[[4L]])) )
  den <- sqrt(n * ( (1.0 / max(1.0, e[[1L]])) +
                    (1.0 / max(1.0, e[[2L]])) +
                    (1.0 / max(1.0, e[[3L]])) +
                    (1.0 / max(1.0, e[[4L]])) ))
  return(num / den)
}

coexMatrix <- function(obs, exp, n, s) {
  coex <- matrix(NA, s, s)
  for (i in (1L:s)) for (j in (i:s)) {
    o <- c(obs[[1L]][i, j], obs[[2L]][i, j], obs[[3L]][i, j], obs[[4L]][i, j])
    e <- c(exp[[1L]][i, j], exp[[2L]][i, j], exp[[3L]][i, j], exp[[4L]][i, j])
    coex[i, j] <- coexPoint(o, e, n)
  }
  return(as.matrix(forceSymmetric(coex)))
}

test_that("Calculations on genes", {
  set.seed(137L)

  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- clean(obj)

  obj <- estimateLambdaLinear(obj)

  lambda <- getLambda(obj)
  mu <- getMu(obj)

  expect_identical(dim(mu), dim(getRawData(obj)))
  expect_equal(mu[ 1L,  1L], lambda[[ 1L]] * getNu(obj)[[ 1L]],
               ignore_attr = TRUE)
  expect_equal(mu[10L,  1L], lambda[[10L]] * getNu(obj)[[ 1L]],
               ignore_attr = TRUE)
  expect_equal(mu[ 1L, 20L], lambda[[ 1L]] * getNu(obj)[[20L]],
               ignore_attr = TRUE)
  expect_equal(mu[10L, 10L], lambda[[10L]] * getNu(obj)[[10L]],
               ignore_attr = TRUE)

  c(observedYY, observedY) %<-%
    observedContingencyTablesYY(obj, actOnCells = FALSE, asDspMatrices = FALSE)

  expect_s4_class(observedYY, "dsyMatrix")
  expect_identical(dim(observedYY), rep(getNumGenes(obj), 2L))
  expect_identical(diag(as.matrix(observedYY)), observedY)
  expect_length(observedY, getNumGenes(obj))
  expect_equal(observedY, c(20L, rep(10L, 9L)), ignore_attr = TRUE)

  observed <- observedContingencyTables(obj, actOnCells = FALSE,
                                        asDspMatrices = FALSE)
  c(observedNN, observedNY, observedYN, .) %<-% observed

  expect_identical(unlist(lapply(c(observedNN, observedNY, observedYN), dim)),
                   rep(dim(observedYY), 3L))
  expect_equal(diag(as.matrix(observedNN)), c(0L, rep(10L, 9L)),
               ignore_attr = TRUE)
  expect_identical(observedNY, t(observedYN))
  expect_equal(diag(as.matrix(observedNY)), rep(0L, 10L), ignore_attr = TRUE)
  expect_equal(as.matrix(observedNN + observedNY + observedYN + observedYY),
               matrix(getNumCells(obj), nrow = getNumGenes(obj),
                      ncol = getNumGenes(obj)),
               ignore_attr = TRUE)

  obj <- estimateDispersionBisection(obj, cores = 4L, chunkSize = 4L)

  c(expectedNN, expectedN) %<-%
    expectedContingencyTablesNN(obj, actOnCells = FALSE, asDspMatrices = TRUE)

  expect_s4_class(expectedNN, "dspMatrix")
  expect_identical(dim(expectedNN), rep(getNumGenes(obj), 2L))
  expect_length(expectedN, getNumGenes(obj))
  expect_equal(expectedN, c(0L, rep(10L, 9L)),
               ignore_attr = TRUE, tolerance = 1e-4)

  expected <- expectedContingencyTables(obj, actOnCells = FALSE,
                                        asDspMatrices = TRUE)
  c(., expectedNY, expectedYN, expectedYY) %<-% expected

  expect_identical(substring(names(observed), 9L),
                   substring(names(expected), 9L))
  expect_identical(unlist(lapply(c(expectedYY, expectedNY, expectedYN), dim)),
                   rep(dim(expectedNN), 3L))
  expect_equal(as.matrix(expectedNN + expectedNY + expectedYN + expectedYY),
               matrix(getNumCells(obj), nrow = getNumGenes(obj),
                      ncol = getNumGenes(obj)),
               ignore_attr = TRUE)

  # take a gene pair ensuring to poll only the upper triangle side of the
  # matrices as the flag 'asDspMatrices = TRUE' makes them incorrect on the
  # other side
  e1 <- sample(getNumGenes(obj), 1L)
  e2 <- sample(getNumGenes(obj), 1L)
  g1 <- getGenes(obj)[[min(e1, e2)]]
  g2 <- getGenes(obj)[[max(e1, e2)]]
  c(gpObs, gpExp) %<-% contingencyTables(obj, g1, g2)

  expect_identical(as.vector(gpObs), c(observedYY[g1, g2], observedYN[g1, g2],
                                       observedNY[g1, g2], observedNN[g1, g2]))
  expect_equal(as.vector(gpExp), c(expectedYY[g1, g2], expectedYN[g1, g2],
                                   expectedNY[g1, g2], expectedNN[g1, g2]),
               tolerance = 1.0e-12)

  zeroOne <- getZeroOneProj(obj)
  probZero <- getProbabilityOfZero(obj)

  gce <- calculateGenesCE(obj)

  expect_named(gce, getGenes(obj))
  expect_identical(gce[[1L]], 0.0)
  expect_equal(gce, crossEntrVector(zeroOne, probZero), ignore_attr = TRUE)

  expect_no_warning({
    obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = FALSE,
                         returnPPFract = TRUE)
  })
  expect_null(getOption("COTAN.TorchWarning"))

  legacyCoex <- getGenesCoex(obj, zeroDiagonal = FALSE)

  expect_true(isCoexAvailable(obj, actOnCells = FALSE, ignoreSync = FALSE))
  expect_identical(dim(legacyCoex), rep(getNumGenes(obj), 2L))
  expect_identical(getGenesCoex(obj)[1L, 1L], 0.0)
  expect_equal(abs(as.vector(legacyCoex[-1L, -1L])),
               rep(1.0, 81L), tolerance = 0.01)

  expect_equal(as.matrix(legacyCoex),
               coexMatrix(observed, expected,
                          getNumCells(obj), getNumGenes(obj)),
               tolerance = 0.001, ignore_attr = TRUE)

  expect_equal(getMetadataDataset(obj)[[1L]], datasetTags()[c(5L, 6L, 7L)],
               ignore_attr = TRUE)
  expect_identical(getMetadataElement(obj, datasetTags()[["gbad"]]),
                   paste0(10.0 / 55.0))

  suppressWarnings({
    obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = TRUE)
  })

  torchCoex <- getGenesCoex(obj, zeroDiagonal = FALSE)

  expect_true(isCoexAvailable(obj, actOnCells = FALSE, ignoreSync = FALSE))
  expect_identical(dim(torchCoex), rep(getNumGenes(obj), 2L))
  expect_identical(getGenesCoex(obj)[1L, 1L], 0.0)
  expect_equal(abs(as.vector(torchCoex[-1L, -1L])),
               rep(1.0, 81L), tolerance = 0.01)

  expect_equal(as.matrix(torchCoex),
               coexMatrix(observed, expected,
                          getNumCells(obj), getNumGenes(obj)),
               tolerance = 0.001, ignore_attr = TRUE)

  expect_equal(getMetadataDataset(obj)[[1L]], datasetTags()[c(5L, 6L, 7L)],
               ignore_attr = TRUE)
  expect_identical(getMetadataElement(obj, datasetTags()[["gbad"]]),
                   paste0(NA))

  genesSample1 <- sample(getNumGenes(obj), 3L)
  partialCoex1 <- calculatePartialCoex(obj, genesSample1)

  # These need tolerance
  expect_equal(partialCoex1,
               legacyCoex[, sort(genesSample1)],
               tolerance = 1e-12)

  genesSample2 <- getGenes(obj)[sample(getNumGenes(obj), 3L)]
  partialCoex2 <- calculatePartialCoex(obj, genesSample2,
                                       optimizeForSpeed = FALSE)

  expect_equal(partialCoex2, legacyCoex[, sort(genesSample2), drop = FALSE],
               tolerance = 1e-12)
})


test_that("Calculations on cells", {
  set.seed(137L)

  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- clean(obj)

  c(observedYY, observedY) %<-%
    observedContingencyTablesYY(obj, actOnCells = TRUE, asDspMatrices = TRUE)

  expect_s4_class(observedYY, "dspMatrix")
  expect_identical(dim(observedYY), rep(getNumCells(obj), 2L))
  expect_identical(diag(as.matrix(observedYY)), observedY)
  expect_length(observedY, getNumCells(obj))
  expect_equal(observedY, rep(c(7L, 4L), 10L), ignore_attr = TRUE)

  observed <- observedContingencyTables(obj, actOnCells = TRUE,
                                        asDspMatrices = TRUE)
  c(observedNN, observedNY, observedYN, .) %<-% observed

  expect_identical(unlist(lapply(c(observedNN, observedNY, observedYN), dim)),
                   rep(dim(observedYY), 3L))
  expect_equal(diag(as.matrix(observedNN)), rep(c(3L, 6L), 10L),
               ignore_attr = TRUE)
  expect_equal(diag(as.matrix(observedNY)), rep(0L, 20L), ignore_attr = TRUE)
  expect_equal(as.matrix(observedNN + observedNY + observedYN + observedYY),
               matrix(getNumGenes(obj), nrow = getNumCells(obj),
                      ncol = getNumCells(obj)),
               ignore_attr = TRUE)

  obj <- estimateLambdaLinear(obj)
  obj <- estimateDispersionNuBisection(obj, cores = 4L, chunkSize = 4L,
                                       enforceNuAverageToOne = FALSE)

  c(expectedNN, expectedN) %<-%
    expectedContingencyTablesNN(obj, actOnCells = TRUE, asDspMatrices = FALSE)

  expect_s4_class(expectedNN, "dsyMatrix")
  expect_identical(dim(expectedNN), rep(getNumCells(obj), 2L))
  expect_length(expectedN, getNumCells(obj))
  expect_equal(expectedN, rep(c(3L, 6L), 10L),
               ignore_attr = TRUE, tolerance = 1e-3)

  expected <- expectedContingencyTables(obj, actOnCells = TRUE,
                                        asDspMatrices = FALSE)
  c(., expectedNY, expectedYN, expectedYY) %<-% expected

  expect_identical(unlist(lapply(c(expectedYY, expectedNY, expectedYN), dim)),
                   rep(dim(expectedNN), 3L))
  expect_identical(expectedNY, t(expectedYN))
  expect_equal(as.matrix(expectedNN + expectedNY + expectedYN + expectedYY),
               matrix(getNumGenes(obj), nrow = getNumCells(obj),
                      ncol = getNumCells(obj)),
               ignore_attr = TRUE)

  expect_warning({
    obj <- calculateCoex(obj, actOnCells = TRUE,
                         optimizeForSpeed = TRUE,
                         returnPPFract = TRUE)
  })

  genesCoexInSync <- getMetadataElement(obj, datasetTags()[["gsync"]])
  cellsCoexInSync <- getMetadataElement(obj, datasetTags()[["csync"]])

  expect_identical(c(genesCoexInSync, cellsCoexInSync), c("FALSE", "TRUE"))

  expect_true(isCoexAvailable(obj, actOnCells = TRUE))
  expect_identical(dim(getCellsCoex(obj)), rep(getNumCells(obj), 2L))

  # as all cells are repeated altenating
  expect_true(
    all(abs(getCellsCoex(obj, zeroDiagonal = FALSE)[, seq_len(getNumCells(obj))
                                                        %% 2L == 1L] -
              getCellsCoex(obj, zeroDiagonal = FALSE)[, 1L]) < 1e-12))
  expect_true(
    all(abs(getCellsCoex(obj, zeroDiagonal = FALSE)[, seq_len(getNumCells(obj))
                                                        %% 2L == 0L] -
            getCellsCoex(obj, zeroDiagonal = FALSE)[, 2L]) < 1e-12))

  expect_equal(as.matrix(getCellsCoex(obj, zeroDiagonal = FALSE)),
               coexMatrix(observed, expected, getNumGenes(obj),
                          getNumCells(obj)),
               tolerance = 0.001, ignore_attr = TRUE)

  expect_equal(getMetadataDataset(obj)[[1L]], datasetTags()[c(5L, 6L, 8L)],
               ignore_attr = TRUE)
  expect_identical(getMetadataElement(obj, datasetTags()[["cbad"]]), paste0(0L))

  cellsSample1 <- sample(getNumCells(obj), 3L)
  partialCoex1 <- calculatePartialCoex(obj, cellsSample1,
                                       actOnCells = TRUE)

  # These need tolerance
  expect_equal(partialCoex1,
               getCellsCoex(obj, zeroDiagonal = FALSE)[, sort(cellsSample1)],
               tolerance = 1e-12)

  cellsSample2 <- getCells(obj)[sample(getNumCells(obj), 3L)]
  partialCoex2 <- calculatePartialCoex(obj, cellsSample2, actOnCells = TRUE,
                                       optimizeForSpeed = FALSE)

  expect_equal(partialCoex2,
               getCellsCoex(obj, cellsSample2, zeroDiagonal = FALSE),
               tolerance = 1e-12)
})


test_that("Coex", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- clean(obj)

  obj <- estimateLambdaLinear(obj)
  obj <- estimateDispersionNuBisection(obj, cores = 4L, chunkSize = 4L,
                                       enforceNuAverageToOne = FALSE)

  expect_no_warning({
    obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = FALSE)
  })
  expect_true(isCoexAvailable(obj))
  expect_identical(dim(getGenesCoex(obj)), rep(getNumGenes(obj), 2L))

  S <- as.matrix(calculateS(obj))
  G <- as.matrix(calculateG(obj))

  expect_identical(dim(S), dim(getGenesCoex(obj)))
  expect_identical(dim(G), dim(getGenesCoex(obj)))
  expect_equal(diag(S), rep(0L, nrow(S)), ignore_attr = TRUE)
  expect_identical(diag(S), diag(G))
  diag(S) <- 1.0
  diag(G) <- 1.4
  expect_equal(S[-1L, 1L] / G[-1L, 1L], rep(11L, 9L), tolerance = 1e-3,
               ignore_attr = TRUE)
  expect_true(all(((1.4 * S[-1L, -1L]) / G[-1L, -1L]) < 1.2))
  expect_true(all((G[-1L, -1L] / (1.4 * S[-1L, -1L])) < 1.2))

  pVS <- calculatePValue(obj, statType = "S")[2L:5L, 6L:9L]
  pVG <- calculatePValue(obj, statType = "G",
                         geneSubsetCol = getGenes(obj)[6L:9L],
                         geneSubsetRow = getGenes(obj)[2L:5L])

  expect_identical(dim(pVS), dim(pVG))
  expect_true(all((pVS / (pVG * 50.0)) < 1.1))
  expect_true(all(((pVG * 50.0) / pVS) < 2.5))

  GDI_S <- calculateGDI(obj, statType = "S")
  GDI_G <- calculateGDI(obj, statType = "G")

  expect_identical(dim(GDI_S), as.integer(c(getNumGenes(obj), 3L)))
  expect_identical(dim(GDI_S), dim(GDI_G))
  expect_identical(colnames(GDI_S), c("sum.raw.norm", "GDI", "exp.cells"))
  expect_identical(colnames(GDI_S), colnames(GDI_G))
  expect_identical(GDI_S[[1L]], GDI_G[[1L]])
  expect_equal(GDI_S[[2L]], GDI_G[[2L]], tolerance = 0.1)
  expect_identical(GDI_S[[3L]], GDI_G[[3L]])
  expect_equal(GDI_S[[3L]],
               c(100L, rep(50L, getNumGenes(obj) - 1L)), ignore_attr = TRUE)

  GDI_S_2 <- calculateGDIGivenS(calculateS(obj))
  GDI_S_3 <- calculateGDIGivenCorr(getGenesCoex(obj),
                                   numDegreesOfFreedom = getNumCells(obj))

  expect_equal(GDI_S[["GDI"]], GDI_S_2, ignore_attr = TRUE)
  expect_identical(GDI_S_2, GDI_S_3)
})


test_that("Coex vs saved results", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj,
                       cores = 6L,
                       optimizeForSpeed = FALSE,
                       deviceStr = "cpu",
                       saveObj = FALSE,
                       outDir = tm)

  genesCoexInSync <- getMetadataElement(obj, datasetTags()[["gsync"]])
  cellsCoexInSync <- getMetadataElement(obj, datasetTags()[["csync"]])

  expect_identical(c(genesCoexInSync, cellsCoexInSync), c("TRUE", "FALSE"))

  suppressWarnings({
    obj2 <- automaticCOTANObjectCreation(raw = test.dataset,
                                         GEO = " ",
                                         sequencingMethod = "artificial",
                                         sampleCondition = "test",
                                         cores = 6L,
                                         optimizeForSpeed = FALSE,
                                         deviceStr = "cpu",
                                         saveObj = FALSE,
                                         outDir = tm)
  })

  expect_identical(obj2, obj)

  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))

  coex_test <- readRDS(file.path(getwd(), "coex.test.RDS"))

  expect_true(isCoexAvailable(obj))
  expect_equal(getGenesCoex(obj, genes = genes.names.test,
                            zeroDiagonal = FALSE),
               coex_test, tolerance = 1.0e-12)

  pval <- calculatePValue(obj, geneSubsetCol = genes.names.test)

  pval_exp <- readRDS(file.path(getwd(), "pval.test.RDS"))
  diag(pval_exp[genes.names.test, ]) <- 1L
  expect_equal(pval, pval_exp, tolerance = 1.0e-12)

  GDI <- calculateGDI(obj)[genes.names.test, ]

  GDI_exp <- readRDS(file.path(getwd(), "GDI.test.RDS"))

  expect_equal(GDI, GDI_exp, tolerance = 1.0e-12)

  # Torch CPU
  suppressWarnings({
    obj3 <- automaticCOTANObjectCreation(raw = test.dataset,
                                         GEO = " ",
                                         sequencingMethod = "artificial",
                                         sampleCondition = "test",
                                         cores = 6L,
                                         optimizeForSpeed = TRUE,
                                         deviceStr = "cpu",
                                         saveObj = FALSE,
                                         outDir = tm)
  })

  expect_equal(obj3@genesCoex, obj@genesCoex, tolerance = 5.0e-6)

  expect_true(isCoexAvailable(obj3))
  expect_equal(getGenesCoex(obj3, genes = genes.names.test,
                            zeroDiagonal = FALSE),
               coex_test, tolerance = 5.0e-6)

  pval <- calculatePValue(obj3, geneSubsetCol = genes.names.test)

  expect_equal(pval, pval_exp, tolerance = 5.0e-6)

  GDI <- calculateGDI(obj3)[genes.names.test, ]

  expect_equal(GDI, GDI_exp, tolerance = 5.0e-6)

  # Torch GPU
  suppressWarnings({
    obj4 <- automaticCOTANObjectCreation(raw = test.dataset,
                                         GEO = " ",
                                         sequencingMethod = "artificial",
                                         sampleCondition = "test",
                                         cores = 6L,
                                         optimizeForSpeed = TRUE,
                                         deviceStr = "cuda",
                                         saveObj = FALSE,
                                         outDir = tm)
  })

  expect_equal(obj4@genesCoex, obj@genesCoex, tolerance = 5.0e-6)

  expect_true(isCoexAvailable(obj4))
  expect_equal(getGenesCoex(obj4, genes = genes.names.test,
                            zeroDiagonal = FALSE),
               coex_test, tolerance = 5.0e-6)

  pval <- calculatePValue(obj4, geneSubsetCol = genes.names.test)

  expect_equal(pval, pval_exp, tolerance = 5.0e-6)

  GDI <- calculateGDI(obj4)[genes.names.test, ]

  expect_equal(GDI, GDI_exp, tolerance = 5.0e-6)
})


test_that("Coex with negative dispersion genes", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  cells.names.test <- readRDS(file.path(getwd(), "cells.names.test.RDS"))
  cellsToDrop <- getCells(obj)[!getCells(obj) %in% cells.names.test]
  obj <- dropGenesCells(obj, cells = cellsToDrop)

  obj <- proceedToCoex(obj, cores = 6L, calcCoex = FALSE, saveObj = FALSE)

  expect_true(any(getDispersion(obj) < 0.0))

  expect_no_warning({
    obj <- calculateCoex(obj, optimizeForSpeed = FALSE, deviceStr = "cpu")
  })
  coex1 <- getGenesCoex(obj, zeroDiagonal = FALSE)

  suppressWarnings({
    obj <- calculateCoex(obj, optimizeForSpeed = TRUE, deviceStr = "cpu")
  })
  coex2 <- getGenesCoex(obj, zeroDiagonal = FALSE)

  suppressWarnings({
    obj <- calculateCoex(obj, optimizeForSpeed = TRUE, deviceStr = "cuda")
  })
  coex3 <- getGenesCoex(obj, zeroDiagonal = FALSE)

  expect_equal(coex1, coex2, tolerance = 1e-7)
  expect_equal(coex1, coex3, tolerance = 1e-7)
  expect_equal(coex2, coex3, tolerance = 5e-7)

  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
                       G2 = c("g-000300", "g-000330", "g-000660"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  hmDF <- singleHeatmapDF(obj, genesLists = groupMarkers, sets = 2L:3L,
                          pValueThreshold = 0.05)

  numOtherGenes <- sum(lengths(groupMarkers)[2L:3L])
  expectedColNames <- c("g2", "g1", "coex", "pValue",
                        "cond", "type", "absent", "fe_genes")
  expect_identical(nrow(hmDF), numOtherGenes * length(groupMarkers[[1L]]))
  expect_setequal(hmDF[["g1"]], groupMarkers[[1L]])
  expect_setequal(hmDF[["g2"]], unlist(groupMarkers[2L:3L]))
  expect_setequal(hmDF[["type"]], names(groupMarkers)[2L:3L])
  expect_setequal(hmDF[["cond"]], "test")
  expect_identical(colnames(hmDF), expectedColNames)
  expect_identical(hmDF[["absent"]], hmDF[["g2"]] == "g-000660")
  expect_true(all(hmDF[["pValue"]] > 0.05))
  expect_true(all(hmDF[["coex"]] == 0.0))
})
