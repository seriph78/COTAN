library(Matrix)

test_that("Logging levels", {
  currentLevel <- getOption("COTAN.LogLevel", default = 1)

  suppressMessages(setLoggingLevel(0))

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

  # clean raw.dataset if necessary
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
