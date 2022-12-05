
test_that("Logging levels work", {
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

test_that("Adding columns to data.frames is OK",{
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

test_that("funProbZero gives the correct values", {

})
