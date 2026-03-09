

test_that("Raw and Clean plots", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  expect_no_condition(plot(ECDPlot(obj)))
  expect_no_condition(plot(cellSizePlot(obj)))
  expect_no_condition(plot(genesSizePlot(obj)))
  expect_no_condition(plot(scatterPlot(obj)))

  gpp <- genesPercentagePlot(obj, genes = getGenes(obj)[1:30])
  expect_identical(dim(gpp[["sizes"]]), c(ncol(test.dataset), 5L))
  expect_identical(colnames(gpp[["sizes"]]),
                   c("sizes", "n", "sample", "sum", "percentage"))
  expect_no_condition(plot(gpp[["plot"]]))


  mpp <- mitochondrialPercentagePlot(obj, genePrefix = "g-0000")
  expect_identical(dim(mpp[["sizes"]]), c(ncol(test.dataset), 5L))
  expect_identical(colnames(mpp[["sizes"]]),
                   c("sizes", "n", "sample", "sum.mit", "mit.percentage"))
  expect_no_condition(plot(mpp[["plot"]]))

  obj <- clean(obj)
  clPlots <- cleanPlots(obj, includePCA = TRUE)

  expect_identical(names(clPlots),
                   c("pcaCells", "pcaCellsData", "genes",
                     "UDE", "nu", "zoomedNu"))
  expect_identical(dim(clPlots[["pcaCellsData"]]), c(ncol(test.dataset), 6L))
  expect_identical(colnames(clPlots[["pcaCellsData"]]),
                   c("PC1", "PC2", "PC3", "PC4", "PC5", "groups"))
  expect_no_condition(plot(clPlots[["pcaCells"]]))
  expect_no_condition(plot(clPlots[["pcaCellsData"]]))
  expect_no_condition(plot(clPlots[["genes"]]))
  expect_no_condition(plot(clPlots[["UDE"]]))
  expect_no_condition(plot(clPlots[["nu"]]))
  expect_no_condition(plot(clPlots[["zoomedNu"]]))

  batch <- factor(rep(c("L", "H"), each = getNumCells(obj) / 2L))
  names(batch) <- getCells(obj)
  obj <- addCondition(obj, condName = "H/L", conditions = batch)

  expect_no_condition(plot(ECDPlot(obj, condName = "H/L")))
  expect_no_condition(plot(cellSizePlot(obj, condName = "H/L")))
  expect_no_condition(plot(genesSizePlot(obj, condName = "H/L")))
  expect_no_condition(plot(scatterPlot(obj, condName = "H/L")))

  gpp <- genesPercentagePlot(obj, genes = getGenes(obj)[1:30], condName = "H/L")
  expect_identical(dim(gpp[["sizes"]]), c(ncol(test.dataset), 5L))
  expect_identical(colnames(gpp[["sizes"]]),
                   c("sizes", "n", "sample", "sum", "percentage"))
  expect_no_condition(plot(gpp[["plot"]]))

  mpp <-
    mitochondrialPercentagePlot(obj, genePrefix = "g-0000", condName = "H/L")
  expect_identical(dim(mpp[["sizes"]]), c(ncol(test.dataset), 5L))
  expect_identical(colnames(mpp[["sizes"]]),
                   c("sizes", "n", "sample", "sum.mit", "mit.percentage"))
  expect_no_condition(plot(mpp[["plot"]]))
})
