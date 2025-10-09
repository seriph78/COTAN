tm <- tempdir()
stopifnot(file.exists(tm))

library(Matrix)
library(rlang)

options(parallelly.fork.enable = TRUE)

test_that("Logging", {
  logPath <- file.path(tm, "COTAN_Test.log")

  suppressMessages({
    currentLevel <- setLoggingLevel(0L)
  })
  suppressMessages({
    currentFile  <- setLoggingFile(logPath)
  })

  expect_false(is.null(currentLevel[[1L]]))
  expect_null(currentFile[[1L]])
  expect_true(file.exists(logPath))

  expect_no_message(suppressMessages(logThis("This should not appear",
                                             logLevel = 0L, appendLF = FALSE)))

  expect_message(   logThis("This should appear",     logLevel = 0L))
  expect_no_message(logThis("This should not appear", logLevel = 1L))

  suppressMessages(setLoggingLevel(1L))

  expect_message(   logThis("This should appear",     logLevel = 1L))
  expect_no_message(logThis("This should not appear"))

  suppressMessages(setLoggingLevel(3L))
  suppressMessages(setLoggingFile(logPath))

  expect_message(logThis("This should appear",        logLevel = 3L))

  # restore logging status
  suppressMessages(setLoggingLevel(currentLevel))
  suppressMessages(setLoggingFile(""))

  expect_equal(R.utils::countLines(logPath), 5L, ignore_attr = TRUE)
  file.remove(logPath)
})


test_that("Clusterizations manipulations", {
  set.seed(1675787192L)
  elemValues <- paste0("", as.roman(sample(7L, 100L, replace = TRUE)))
  elemNames <- paste0("el_", 1L:100L)

  clusters <- as.data.frame(list("a" = elemNames, "b" = elemValues))
  clusters <- asClusterization(clusters, elemNames)

  expect_s3_class(clusters, "factor")
  expect_setequal(levels(clusters), paste0("", as.roman(1L:7L)))
  expect_named(clusters, elemNames)

  clustersList <- toClustersList(clusters)

  expect_length(clustersList, 7L)
  expect_setequal(lengths(clustersList), as.vector(table(clusters)))

  clusters2 <- fromClustersList(clustersList, elemNames)

  expect_identical(clusters2, clusters)

  clusters3 <- fromClustersList(clustersList, elemNames = NULL)
  expect_equal(table(clusters2), table(clusters3), ignore_attr = TRUE)

  positions <- groupByClusters(clusters)

  expect_identical(clusters2[positions], clusters3)

  expect_identical(groupByClusters(clusters2), positions)
  expect_identical(groupByClustersList(elemNames, clustersList), positions)

  # cause mismatches between the element names and the clusterization
  elemNames <- append(elemNames, paste0("el_", 201L:210L), after = 20L)[1L:100L]

  expect_setequal(fromClustersList(clustersList, elemNames)[21L:30L], "-1")

  expect_identical(groupByClustersList(elemNames, clustersList)[91L:100L],
                   (21L:30L))

  clusterM1 <- mergeClusters(clusters, names = as.roman(c(5L, 1L)),
                             mergedName = "I'V")

  expect_identical(levels(clusterM1)[[nlevels(clusterM1)]], "I'V")
  expect_identical(table(clusterM1)[["I'V"]], sum(table(clusters)[c(1L, 5L)]))

  clusterM2 <-
    multiMergeClusters(clusters3, namesList = list(as.roman(c(1L, 5L)),
                                                   as.roman(c(6L, 2L, 4L))))

  expect_setequal(levels(clusterM2),
                  c("I_V-merge", "II_IV_VI-merge", "III", "VII"))
  expect_identical(sum(clusterM2 == "I_V-merge"), sum(clusterM1 == "I'V"))

  niceClusters <- niceFactorLevels(clusters)
  expect_identical(max(nchar(factorToVector(niceClusters))), 3L)
  expect_identical(min(nchar(factorToVector(niceClusters))), 3L)
  expect_true(all(endsWith(factorToVector(niceClusters),
                           factorToVector(clusters))))

  levels(niceClusters) <- c(1L:3L, 11L:13L, 100L)
  niceClusters <- niceFactorLevels(niceClusters)
  expect_identical(max(nchar(factorToVector(niceClusters))), 3L)
  expect_identical(min(nchar(factorToVector(niceClusters))), 3L)
  expect_setequal(as.integer(levels(niceClusters)), c(1L:3L, 11L:13L, 100L))
})


test_that("Adding/extracting columns to/from data.frames", {
  df <- data.frame()

  df <- setColumnInDF(df, colName = "constant", colToSet = rep(1L, 10L))
  expect_identical(dim(df), c(10L, 1L))
  expect_setequal(df[["constant"]], 1L)

  df <- setColumnInDF(df, colToSet = (1L:10L),
                      colName = "sequence", rowNames = LETTERS[1L:10L])
  expect_identical(rownames(df), LETTERS[1L:10L])
  expect_identical(colnames(df), c("constant", "sequence"))
  expect_identical(getColumnFromDF(df, "constant"),
                   rep_named(LETTERS[1L:10L], 1L))
  expect_identical(getColumnFromDF(df, 2L),
                   set_names(seq_len(10L), LETTERS[1L:10L]))

  df <- setColumnInDF(df, colName = "constant", colToSet = rep(2L, 10L))
  expect_identical(colnames(df), c("constant", "sequence"))
  expect_setequal(df[["constant"]], 2L)
})


test_that("funProbZeroNegBin", {
  # Cases with mu = 0 are not actually in use
  expect_identical(funProbZeroNegBin(-Inf, 0.0), 0.0)
  expect_identical(funProbZeroNegBin(-1.0, 0.0), 1.0)
  expect_identical(funProbZeroNegBin( 0.0, 0.0), 1.0)
  expect_identical(funProbZeroNegBin(1e-6, 0.0), 1.0)
  expect_identical(funProbZeroNegBin( 1.0, 0.0), 1.0)
  expect_identical(funProbZeroNegBin(10.0, 0.0), 1.0)
  expect_identical(funProbZeroNegBin( Inf, 0.0), 1.0)

  # Cases with infinite disp can happen
  expect_identical(funProbZeroNegBin(-Inf, 1.0),            0.0)
  expect_identical(funProbZeroNegBin(-1.0, 1.0),      exp(-2.0))
  expect_identical(funProbZeroNegBin( 0.0, 1.0),      exp(-1.0))
  expect_identical(funProbZeroNegBin(1e-6, 1.0), exp(-0.999999))
  expect_identical(funProbZeroNegBin( 1.0, 1.0),      1.0 / 2.0)
  expect_identical(funProbZeroNegBin(10.0, 1.0),      11.0^-0.1)
  expect_identical(funProbZeroNegBin( Inf, 1.0),            1.0)

  # Cases with mu = Inf are not actually in use
  expect_identical(funProbZeroNegBin(-Inf, Inf), 0.0)
  expect_identical(funProbZeroNegBin(-1.0, Inf), 0.0)
  expect_identical(funProbZeroNegBin( 0.0, Inf), NaN)
  expect_identical(funProbZeroNegBin(1e-6, Inf), NaN)
  expect_identical(funProbZeroNegBin( 1.0, Inf), 0.0)
  expect_identical(funProbZeroNegBin(10.0, Inf), 0.0)
  expect_identical(funProbZeroNegBin( Inf, Inf), 1.0)
})


test_that("funProbZeroNegBin with matrices", {
  mu <- matrix((1L:25L) / 7.0, nrow = 10L, ncol = 10L)
  disp <- c(-Inf, -1.0/3.0, 1e-6, (1L:6L) / 3.0, Inf)

  p <- funProbZeroNegBin(disp, mu)

  expect_identical(dim(p), dim(mu))

  expect_identical(p[ 1L,  1L], funProbZeroNegBin(disp[[ 1L]], mu[ 1L,  1L]))
  expect_identical(p[ 1L, 10L], funProbZeroNegBin(disp[[ 1L]], mu[ 1L, 10L]))

  expect_identical(p[ 2L,  2L], funProbZeroNegBin(disp[[ 2L]], mu[ 2L,  2L]))
  expect_identical(p[ 3L,  7L], funProbZeroNegBin(disp[[ 3L]], mu[ 3L,  7L]))
  expect_identical(p[ 6L,  4L], funProbZeroNegBin(disp[[ 6L]], mu[ 6L,  4L]))

  expect_identical(p[10L,  1L], funProbZeroNegBin(disp[[10L]], mu[10L,  1L]))
  expect_identical(p[10L, 10L], funProbZeroNegBin(disp[[10L]], mu[10L, 10L]))
})


test_that("funProbZeroMixPoi", {
  # Cases with mu = 0 are not actually in use
  expect_identical(funProbZeroMixPoi(0.0, 1.0), 0.0 + 1.0 * exp(-1.0))
  expect_identical(funProbZeroMixPoi(1.0, 1.0), 1.0 + 0.0 * exp(-1.0))
  expect_identical(funProbZeroMixPoi(0.5, 2.0), 0.5 + 0.5 * exp(-2.0))
  expect_identical(funProbZeroMixPoi(0.2, 0.5), 0.2 + 0.8 * exp(-0.5))
  expect_identical(funProbZeroMixPoi(0.8, 4.0), 0.8 + 0.2 * exp(-4.0))
  expect_identical(funProbZeroMixPoi(0.6, 0.8), 0.6 + 0.4 * exp(-0.8))
})


test_that("funProbZeroMixPoi with matrices", {
  mu <- matrix((1L:25L) / 7.0, nrow = 10L, ncol = 10L)
  pi <- (0L:9L) / 9.0

  p <- funProbZeroMixPoi(pi, mu)

  expect_identical(dim(p), dim(mu))

  expect_identical(p[ 1L,  1L], funProbZeroMixPoi(pi[[ 1L]], mu[ 1L,  1L]))
  expect_identical(p[ 1L, 10L], funProbZeroMixPoi(pi[[ 1L]], mu[ 1L, 10L]))

  expect_identical(p[ 3L,  8L], funProbZeroMixPoi(pi[[ 3L]], mu[ 3L,  8L]))
  expect_identical(p[ 6L,  4L], funProbZeroMixPoi(pi[[ 6L]], mu[ 6L,  4L]))

  expect_identical(p[10L,  1L], funProbZeroMixPoi(pi[[10L]], mu[10L,  1L]))
  expect_identical(p[10L, 10L], funProbZeroMixPoi(pi[[10L]], mu[10L, 10L]))
})


test_that("dispersionSolvers", {
  lambda   <- c(3.0, 1.75, 2.0)
  nu       <- rep(c(0.5, 1.5), 5L)
  sumZeros <- c(0.0, 5.0, 1.0)

  d1 <- mapply(
    FUN       = dispersionBisection,
    lambda    = lambda,
    sumZeros  = sumZeros,
    MoreArgs  = list("nu" = nu),
    SIMPLIFY  = TRUE,    # collapse to atomic vector
    USE.NAMES = FALSE)

  expect_identical(d1, c(-Inf, 1.980468750, -0.646484375))

  expect_lt(max(abs(rowSums(funProbZero(d1, lambda %o% nu)) - sumZeros)),
            1.0e-3)

  d2 <- mapply(
    FUN       = dispersionNewton,
    lambda    = lambda,
    sumZeros  = sumZeros,
    MoreArgs  = list("nu" = nu),
    SIMPLIFY  = TRUE,    # collapse to atomic vector
    USE.NAMES = FALSE)

  expect_equal(d2, c(-Inf, 1.980664643877, -0.645937418860), tolerance = 1e-12)

  expect_lt(max(abs(rowSums(funProbZero(d2, lambda %o% nu)) - sumZeros)),
            2.0e-5)

  d1p <- parallelDispersionBisection(1L:3L, sumZeros, lambda, nu)

  expect_identical(d1p, d1)

  d2p <- parallelDispersionNewton(1L:3L, sumZeros, lambda, nu)

  expect_equal(d2p, d2, tolerance = 1e-12)
})

test_that("lambdaNewton", {
  nu <- rep(c(0.5, 1.5), 5L)
  avgNumNonZeros <- c(0.0, 0.1, 0.4, 0.6, 0.8,  1.0)
  avgCounts      <- c(0.0, 0.1, 1.0, 1.3, 3.0, 10.0)
  zeroPi <- (rowSums(exp(-(avgCounts %o% nu))) / 10.0 + avgNumNonZeros > 1.0)
  zeroPi <- (zeroPi | avgNumNonZeros == avgCounts)

  lambda <- mapply(lambdaNewton,
                   avgNumNonZeros = avgNumNonZeros,
                   avgCounts = avgCounts,
                   zeroPi = zeroPi,
                   MoreArgs = list(nu = nu))

  expectedRes <- c(0.0, 0.1067852, 1.9671601, 1.5695487, 3.3949319, +Inf)
  expect_equal(lambda, expectedRes,
               tolerance = 1e-6)

  pi <- ifelse(zeroPi, 0.0, 1.0 - avgCounts / lambda)

  piError <- avgNumNonZeros -
    (1.0 - rowSums(funProbZeroMixPoi(pi, lambda %o% nu)) / 10.0)

  expect_identical(piError[c(1L, 6L)], c(0.0, 0.0))
  expect_lt(max(abs(piError[2L:5L])), 5.0e-5)

  lambda2 <- parallelLambdaNewton(genes = rep_len(TRUE, 6L),
                                  avgNumNonZeros = avgNumNonZeros,
                                  avgCounts = avgCounts,
                                  nu = nu,
                                  zeroPi = zeroPi)

  expect_equal(lambda2, expectedRes, tolerance = 1e-6)
})


test_that("nuBisection", {
  lambda       <- c(5.5,  4.0, 2.0, 1.0, 5.5, 1.5, 3.0, 3.5, 1.0, 4.5)
  dispersion   <- c(-Inf, 4.0, 2.5, 0.9, 5.0, 1.9, 3.5, 4.0, 0.9, 4.5)
  sumZeros     <- c(3L, 6L)
  initialGuess <- c(1.0, 1.5)

  nu <- c(nuBisection(sumZeros = sumZeros[[1L]], lambda = lambda,
                      dispersion = dispersion,
                      initialGuess = initialGuess[[1L]]),
          nuBisection(sumZeros = sumZeros[[2L]], lambda = lambda,
                      dispersion = dispersion,
                      initialGuess = initialGuess[[2L]]))

  expect_identical(nu, c(3.1484375, 0.349365234375))
})


test_that("plotTheme", {
  expect_warning(plotTheme(plotKind = "Wrong", textSize = 12.0))
})


test_that("Raw data normalization", {
  utils::data("test.dataset", package = "COTAN")
  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))
  cells.names.test <- readRDS(file.path(getwd(), "cells.names.test.RDS"))

  raw <- test.dataset[genes.names.test, cells.names.test]

  nu <- readRDS(file.path(getwd(), "nu.test.RDS"))
  raw.norm <- as.matrix(readRDS(file.path(getwd(), "raw.norm.test.RDS")))

  expect_identical(t(t(raw) * (1.0 / nu)), raw.norm)
})


test_that("parallelDist - cosine dissimilarity", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  cd <- as.matrix(calcDist(t(raw), method = "cosine"))

  expect_equal(cd[(row(cd) + col(cd)) %% 2L == 1L], rep(cd[2L, 1L], 200L),
               tolerance = 1.0e-15)
  expect_equal(cd[(row(cd) + col(cd)) %% 2L == 0L], rep(cd[3L, 1L], 200L),
               tolerance = 1.0e-15)

  expect_error(expect_message(expect_message(calcDist(t(raw),
                                                      method = "dummy"))))
})


test_that("pca usage", {
  utils::data("test.dataset", package = "COTAN")

  pcaRaw <- runPCA(x = as.matrix(test.dataset), rank = 10L,
                   BSPARAM = IrlbaParam(), get.rotation = FALSE)[["x"]]

  expect_identical(rownames(pcaRaw), rownames(test.dataset))

  pcaExp <- readRDS(file.path(getwd(), "pca.raw.test.RDS"))

  pcaRaw <- pcaRaw[rownames(pcaExp), ]

  correlations <- c(cor(pcaRaw[, 1L], pcaExp[, 1L]),
                    cor(pcaRaw[, 2L], pcaExp[, 2L]))

  dists <- sqrt(c(sum((pcaRaw[, 1L] - correlations[[1L]] * pcaExp[, 1L])^2L),
                  sum((pcaRaw[, 2L] - correlations[[2L]] * pcaExp[, 2L])^2L)))

  expect_lt(max(dists), 10.0^(-4L))
})


# legacy
test_that("vec2mat_rfast", {
  mat <- matrix(0.0, nrow = 10L, ncol = 10L)
  mat <- Rfast::lower_tri.assign(mat, (1L:55L), diag = TRUE)
  mat <- Rfast::upper_tri.assign(mat,
                                 v = Rfast::upper_tri(Rfast::transpose(mat)))

  colnames(mat) <- paste0("row.", (1L:10L))
  rownames(mat) <- paste0("row.", (1L:10L))

  genes <- paste0("row.", c(1L, 2L, 9L, 10L))

  expect_identical(mat, vec2mat_rfast(mat2vec_rfast(mat)))
  expect_identical(mat[, genes],
                   vec2mat_rfast(mat2vec_rfast(mat), genes = genes))
})


test_that("mat2vec_rfast", {
  names.v <- paste0("raw", (1L:15L))

  vec <- list("genes" = names.v, "values" = 1L:120L)

  expect_equal(vec, mat2vec_rfast(vec2mat_rfast(vec)), ignore_attr = TRUE)
})
