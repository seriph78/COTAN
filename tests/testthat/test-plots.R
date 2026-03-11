

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
  expect_identical(names(gpp), c("plot", "sizes"))
  expect_identical(dim(gpp[["sizes"]]), c(ncol(test.dataset), 5L))
  expect_identical(colnames(gpp[["sizes"]]),
                   c("sizes", "n", "sample", "sum", "percentage"))
  expect_no_condition(plot(gpp[["plot"]]))

  mpp <-
    mitochondrialPercentagePlot(obj, genePrefix = "g-0000", condName = "H/L")
  expect_identical(names(mpp), c("plot", "sizes"))
  expect_identical(dim(mpp[["sizes"]]), c(ncol(test.dataset), 5L))
  expect_identical(colnames(mpp[["sizes"]]),
                   c("sizes", "n", "sample", "sum.mit", "mit.percentage"))
  expect_no_condition(plot(mpp[["plot"]]))
})


test_that("Heatmap plots", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, calcCoex = TRUE, cores = 3L)

  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000138",
                              "g-000150", "g-000160", "g-000170"),
                       G2 = c("g-000300", "g-000330", "g-000450",
                              "g-000460", "g-000470"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  expect_no_warning(
    plot(heatmapPlot(obj, cores = 3L, genesLists = groupMarkers))
  )

  expect_no_warning(
    genesHeatmapPlot(obj, symmetric = TRUE, cores = 3L,
                     primaryMarkers = lapply(groupMarkers, \(g) g[[1L]]))
  )

  expect_no_warning(
    genesHeatmapPlot(obj, symmetric = FALSE, cores = 3L,
                     primaryMarkers = lapply(groupMarkers, \(g) g[[1L]]))
    )

  expect_no_warning(
    genesHeatmapPlot(obj, symmetric = TRUE, cores = 3L,
                     primaryMarkers = lapply(groupMarkers, \(g) g[[1L]]),
                     secondaryMarkers = lapply(groupMarkers, \(g) g[[2L]]))
  )

  expect_no_warning(
    genesHeatmapPlot(obj, symmetric = TRUE, cores = 3L,
                     primaryMarkers = lapply(groupMarkers, \(g) g[[1L]]),
                     secondaryMarkers = lapply(groupMarkers, \(g) g[[2L]]))
  )

  suppressWarnings(obj <- calculateCoex(obj, actOnCells = TRUE))
  expect_no_warning(
    plot(cellsHeatmapPlot(obj, getCells(obj)[seq.int(20, 1000, 20)]))
  )
})


test_that("Clusters plots", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, calcCoex = TRUE, cores = 3L)

  batch <- factor(rep(c("L", "H"), each = getNumCells(obj) / 2L))
  names(batch) <- getCells(obj)
  obj <- addCondition(obj, condName = "H/L", conditions = batch)
  obj <- addClusterization(obj, clName = "batch", batch)

  csdp1 <- clustersSummaryPlot(obj, clName = "batch",
                               plotTitle = "by batch")
  expect_identical(names(csdp1), c("data", "plot"))
  expect_identical(dim(csdp1[["data"]]), c(nlevels(batch), 8L))
  expect_identical(colnames(csdp1[["data"]]),
                   c("batch", "NoCond", "CellNumber", "CellPercentage",
                     "MeanUDE", "MedianUDE", "ExpGenes25", "ExpGenes"))
  expect_no_condition(plot(csdp1[["plot"]]))

  csdp2 <- clustersSummaryPlot(obj, clName = "batch", condName = "H/L",
                               plotTitle = "by batch")
  expect_identical(names(csdp2), c("data", "plot"))
  expect_identical(dim(csdp2[["data"]]), c(nlevels(batch) * nlevels(batch), 8L))
  expect_identical(colnames(csdp2[["data"]]),
                   c("batch", "H/L", "CellNumber", "CellPercentage",
                     "MeanUDE", "MedianUDE", "ExpGenes25", "ExpGenes"))
  expect_no_error(suppressWarnings(plot(csdp2[["plot"]])))

  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000138",
                              "g-000150", "g-000160", "g-000170"),
                       G2 = c("g-000300", "g-000330", "g-000450",
                              "g-000460", "g-000470"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  chpd1 <- clustersMarkersHeatmapPlot(obj, clName = "batch",
                                      groupMarkers = groupMarkers, kCuts = 2L)
  expect_identical(names(chpd1), c("heatmapPlot", "dataScore", "pValues"))
  expect_identical(dim(chpd1[["dataScore"]]),
                   c(length(unlist(groupMarkers)), 2L))
  expect_identical(dim(chpd1[["pValues"]]),
                   c(length(unlist(groupMarkers)), 2L))
  expect_no_warning(
    plot(chpd1[["heatmapPlot"]])
  )

  suppressWarnings(
    chpd2 <-
      clustersMarkersHeatmapPlot(obj, clName = "batch", condNameList = "H/L",
                                 groupMarkers = groupMarkers, kCuts = 2L)
  )
  expect_identical(names(chpd2), c("heatmapPlot", "dataScore", "pValues"))
  expect_identical(dim(chpd2[["dataScore"]]),
                   c(length(unlist(groupMarkers)), 2L))
  expect_identical(dim(chpd2[["pValues"]]),
                   c(length(unlist(groupMarkers)), 2L))
  expect_no_warning(
    plot(chpd2[["heatmapPlot"]])
  )

  cupd1 <- cellsUMAPPlot(obj, dataMethod = "LogLikelihood", clName = "batch",
                         useCoexEigen = TRUE, numComp = 5L)
  expect_identical(names(cupd1), c("plot", "cellsRDM"))
  expect_identical(dim(cupd1[["cellsRDM"]]), c(getNumCells(obj), 5L))
  expect_warning(
    plot(cupd1[["plot"]]),
    regexp = "No shared levels"
  )

  cupd2 <- cellsUMAPPlot(obj, dataMethod = "AdjBinarized",  clName = "batch",
                         useCoexEigen = FALSE, numComp = 5L,
                         genesSel = "HGDI", numGenes = 100)
  expect_identical(names(cupd2), c("plot", "cellsRDM"))
  expect_identical(dim(cupd2[["cellsRDM"]]), c(getNumCells(obj), 5L))
  expect_warning(
    plot(cupd2[["plot"]]),
    regexp = "No shared levels"
  )
})
