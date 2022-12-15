library(Matrix)

test_that("Logging levels", {
  suppressMessages(currentLevel <- setLoggingLevel(0))

  expect_false(is.null(currentLevel))

  expect_no_message(suppressMessages(logThis("This should not appear", logLevel = 0)))

  expect_message(   logThis("This should appear",     logLevel = 0))
  expect_no_message(logThis("This should not appear", logLevel = 1))

  suppressMessages(setLoggingLevel(1))

  expect_message(   logThis("This should appear",     logLevel = 1))
  expect_no_message(logThis("This should not appear"))

  suppressMessages(setLoggingLevel(3))

  expect_message(logThis("This should appear", logLevel = 3))

  # restore logging level
  suppressMessages(setLoggingLevel(currentLevel))
})


test_that("Adding columns to data.frames", {
  df <- data.frame()

  df <- setColumnInDF(df, colName = "constant", colToSet = rep(1, 10))
  expect_equal(dim(df), as.integer(c(10,1)))
  expect_setequal(df[["constant"]], 1)

  df <- setColumnInDF(df, colToSet = c(1:10), colName = "sequence", rowNames = LETTERS[1:10])
  expect_equal(rownames(df), LETTERS[1:10])
  expect_equal(colnames(df), c("constant", "sequence"))

  df <- setColumnInDF(df, colName = "constant", colToSet = rep(2, 10))
  expect_equal(colnames(df), c("constant", "sequence"))
  expect_setequal(df[["constant"]], 2)
})


test_that("funProbZero", {
  # Cases with mu = 0 are not actually in use
  expect_equal(funProbZero(-Inf,   0), NaN)
  expect_equal(funProbZero(  -1,   0),   1)
  expect_equal(funProbZero(   0,   0),   1)
  expect_equal(funProbZero(   1,   0),   1)
  expect_equal(funProbZero(  10,   0),   1)
  expect_equal(funProbZero( Inf,   0), NaN)

  # Cases with infinite disp can happen
  expect_equal(funProbZero(-Inf,   1),          0)
  expect_equal(funProbZero(  -1,   1),    exp(-2))
  expect_equal(funProbZero(   0,   1),    exp(-1))
  expect_equal(funProbZero(   1,   1),        1/2)
  expect_equal(funProbZero(  10,   1), 11^(-1/10))
  expect_equal(funProbZero( Inf,   1),          1)

  # Cases with mu = Inf are not actually in use
  expect_equal(funProbZero(-Inf, Inf),   0)
  expect_equal(funProbZero(  -1, Inf),   0)
  expect_equal(funProbZero(   0, Inf), NaN)
  expect_equal(funProbZero(   1, Inf),   0)
  expect_equal(funProbZero(  10, Inf),   0)
  expect_equal(funProbZero( Inf, Inf),   1)
})


test_that("funProbZero with matrices", {
  mu <- matrix(c(1:25)/7, nrow = 10, ncol = 10)
  disp <- c(-1:8)/3

  p <- funProbZero(disp, mu)

  expect_equal(dim(p),dim(mu))

  expect_equal(p[ 1, 1], funProbZero(disp[ 1], mu[ 1, 1]))
  expect_equal(p[ 1,10], funProbZero(disp[ 1], mu[ 1,10]))

  expect_equal(p[ 3, 7], funProbZero(disp[ 3], mu[ 3, 7]))
  expect_equal(p[ 6, 4], funProbZero(disp[ 6], mu[ 6, 4]))

  expect_equal(p[10, 1], funProbZero(disp[10], mu[10, 1]))
  expect_equal(p[10,10], funProbZero(disp[10], mu[10,10]))
})


test_that("dispersionBisection", {
  lambda   <- c(3,  1.75)
  nu       <- rep(c(0.5, 1.5), 5)
  sumZeros <- c(0, 5)

  d <- c(dispersionBisection(sumZeros = sumZeros[1], lambda = lambda[1], nu = nu),
         dispersionBisection(sumZeros = sumZeros[2], lambda = lambda[2], nu = nu))

  expect_equal(d, c(-Inf, 1.98046875))
})


test_that("nuBisection", {
  lambda       <- c(5.5,  4,2,    1,   5.5, 1.5, 3,   3.5, 1,   4.5)
  dispersion   <- c(-Inf, 4,      2.5, 0.9, 5,   1.9, 3.5, 4,   0.9, 4.5)
  sumZeros     <- c(3, 6)
  initialGuess <- c(1, 1.5)

  nu <- c(nuBisection(sumZeros = sumZeros[1], lambda = lambda,
                      dispersion = dispersion, initialGuess = initialGuess[1]),
          nuBisection(sumZeros = sumZeros[2], lambda = lambda,
                      dispersion = dispersion, initialGuess = initialGuess[2]))

  expect_equal(nu, c(3.1484375, 0.349365234375))
})


test_that("plotTheme", {
  expect_warning(plotTheme(plotKind = "Wrong", textSize = 12))
})


test_that("Cosine dissimilarity", {
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  cd <- as.matrix(cosineDissimilarity(raw))

  expect_equal(cd[(row(cd) + col(cd)) %% 2 == 1], rep(cd[2,1], 200))
  expect_equal(cd[(row(cd) + col(cd)) %% 2 == 0], rep(cd[3,1], 200))
})

test_that("Raw data normalization", {
  utils::data("test.dataset.col", package = "COTAN")
  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))
  cell.names.test <- readRDS(file.path(getwd(), "cell.names.test.RDS"))

  raw <- test.dataset.col[genes.names.test, cell.names.test]

  nu <- readRDS(file.path(getwd(), "nu.test.RDS"))
  raw.norm <- as.matrix(readRDS(file.path(getwd(), "raw.norm.test.RDS")))

  expect_equal(t(t(raw) * (1/nu)), raw.norm)
})


test_that("prcomp_irlba usage", {
  utils::data("raw.dataset", package = "COTAN")
  pca.tb = readRDS(file.path(getwd(),"pca.tb.RDS"))

  # strip names column from raw.dataset if necessary
  if (nrow(raw.dataset) != nrow(pca.tb) &&
      rownames(raw.dataset)[[1]] == "HK") {
    raw.dataset <- raw.dataset[2:nrow(raw.dataset),]
  }

  pca <- irlba::prcomp_irlba(raw.dataset, n = 5)

  pca.raw <- pca$x
  rm(pca)

  rownames(pca.raw) <- rownames(raw.dataset)
  colnames(pca.raw) <- paste0("PC_", c(1:5))

  correlation1.value <- cor(pca.raw[,1], pca.tb[,1])
  correlation2.value <- cor(pca.raw[,2], pca.tb[,2])
  pca.tb[,1] <- correlation1.value*pca.tb[,1]
  pca.tb[,2] <- correlation2.value*pca.tb[,2]
  x1 <- pca.raw[rownames(pca.tb),1] - pca.tb[,1]
  x2 <- pca.raw[rownames(pca.tb),2] - pca.tb[,2]

  dist1 <- sqrt(sum(x1^2))
  dist2 <- sqrt(sum(x2^2))

  expect_true(dist1 < 10^(-4))
  expect_true(dist2 < 10^(-4))
})


# legacy
test_that("vec2mat_rfast", {
  mat <- matrix(0,nrow = 10, ncol = 10)
  mat <- Rfast::lower_tri.assign(mat, c(1:55), diag = T)
  mat <- Rfast::upper_tri.assign(mat, v = Rfast::upper_tri(Rfast::transpose(mat)))

  colnames(mat) <- paste0("row.", c(1:10))
  rownames(mat) <- paste0("row.", c(1:10))

  genes <- paste0("row.", c(1, 2, 9, 10))

  expect_equal(mat, vec2mat_rfast(mat2vec_rfast(mat)))
  expect_equal(mat[, genes], vec2mat_rfast(mat2vec_rfast(mat), genes = genes))
})


test_that("mat2vec_rfast", {
  names.v <- paste0("raw", c(1:15))

  vec <- list("genes" = names.v, "values" =  c(1:120))

  expect_equal(vec, mat2vec_rfast(vec2mat_rfast(vec)))
})
