
# Creates the files to be reloaded by the tests for comparisons

outputTestDatasetCreation <-
  function(testsDir = file.path("tests", "testthat")) {
  #nolint start: object_name_linter
  utils::data("test.dataset", package = "COTAN")
  options(parallelly.fork.enable = TRUE)
  setLoggingLevel(3L)

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(objCOTAN = obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(objCOTAN = obj, cores = 6L, saveObj = FALSE)

  if (FALSE) {
    saveRDS(obj, file = file.path(testsDir, "test.COTAN.RDS"))
  }

  cells.names.test <-
    getCells(objCOTAN = obj)[c(1L:10L, 591L:610L, 991L:1000L)]
  genes.names.test <-
    getGenes(objCOTAN = obj)[c(131L:140L, 291L:310L, 591L: 600L)]
  saveRDS(cells.names.test, file.path(testsDir, "cells.names.test.RDS"))
  saveRDS(genes.names.test, file.path(testsDir, "genes.names.test.RDS"))

  pcaRaw <- runPCA(x = getRawData(objCOTAN = obj), rank = 10L,
                   BSPARAM = IrlbaParam(), get.rotation = FALSE)[["x"]]

  pca.raw.test <- pcaRaw[genes.names.test, ]
  saveRDS(pca.raw.test, file.path(testsDir, "pca.raw.test.RDS"))

  dispersion.test <- getDispersion(objCOTAN = obj)[genes.names.test]
  saveRDS(dispersion.test, file.path(testsDir, "dispersion.test.RDS"))

  raw.norm.test <-
    getNuNormData(objCOTAN = obj)[genes.names.test, cells.names.test]
  saveRDS(raw.norm.test, file.path(testsDir, "raw.norm.test.RDS"))

  coex.test <-
    getGenesCoex(objCOTAN = obj, genes = genes.names.test, zeroDiagonal = FALSE)
  saveRDS(coex.test, file.path(testsDir, "coex.test.RDS"))

  lambda.test <- getLambda(objCOTAN = obj)[genes.names.test]
  saveRDS(lambda.test, file.path(testsDir, "lambda.test.RDS"))

  GDI.test <- calculateGDI(objCOTAN = obj)
  GDI.test <- GDI.test[genes.names.test, ]
  saveRDS(GDI.test, file.path(testsDir, "GDI.test.RDS"))

  nu.test <- getNu(objCOTAN = obj)[cells.names.test]
  saveRDS(nu.test, file.path(testsDir, "nu.test.RDS"))

  pvalues.test <- calculatePValue(objCOTAN = obj,
                                  geneSubsetCol = genes.names.test,
                                  geneSubsetRow = genes.names.test)
  saveRDS(pvalues.test, file.path(testsDir, "pvalues.test.RDS"))

  groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000138"),
                       G2 = c("g-000300", "g-000330", "g-000660"),
                       G3 = c("g-000510", "g-000530", "g-000550",
                              "g-000570", "g-000590"))

  gcsData <- genesCoexSpace(objCOTAN = obj,
                            primaryMarkers = unlist(groupMarkers),
                            numGenesPerMarker = 11L)

  genes.coex.space.test <- gcsData[["GCS"]]
  saveRDS(genes.coex.space.test,
          file.path(testsDir, "genes.coex.space.test.RDS"))

  genesClustersData <-
    establishGenesClusters(objCOTAN = obj,
                           groupMarkers = groupMarkers,
                           numGenesPerMarker = 11L,
                           kCuts = 6L, distance = "cosine",
                           hclustMethod = "ward.D2")

  pca.genes.clusters.test <- genesClustersData[["pca_clusters"]]
  saveRDS(pca.genes.clusters.test,
          file.path(testsDir, "pca.genes.clusters.test.RDS"))

  # Make it a less strict check as it is only for testing
  checker <- new("AdvancedGDIUniformityCheck")
  checker <- shiftCheckerThresholds(checker, 0.1)

  initialResolution <- 1.3
  splitData <- cellsUniformClustering(objCOTAN = obj,
                                      checker = checker,
                                      initialResolution = initialResolution,
                                      useCoexEigen = TRUE,
                                      dataMethod = "LL",
                                      numReducedComp = 50L,
                                      cores = 6L, optimizeForSpeed = TRUE,
                                      deviceStr = "cuda", saveObj = FALSE)

  split.clusters.test <- splitData[["clusters"]]
  saveRDS(split.clusters.test,
          file = file.path(testsDir, "split.clusters.test.RDS"))

  test.dataset.clusters1 <- split.clusters.test
  save(test.dataset.clusters1, compress = TRUE,
       file = file.path("data", "test.dataset.clusters1.rda"))

  obj <- addClusterization(objCOTAN = obj,
                           clName = "split",
                           clusters = splitData[["clusters"]],
                           coexDF = splitData[["coex"]])

  coex.clusters.test <- splitData[["coex"]][genes.names.test, ]
  saveRDS(coex.clusters.test, file.path(testsDir, "coex.clusters.test.RDS"))

  pvalDF <- pValueFromDEA(splitData[["coex"]],
                          getNumCells(objCOTAN = obj),
                          adjustmentMethod = "none")

  pvalues.clusters.test <- pvalDF[genes.names.test, ]
  saveRDS(pvalues.clusters.test,
          file.path(testsDir, "pvalues.clusters.test.RDS"))

  mergedData <- mergeUniformCellsClusters(objCOTAN = obj,
                                          clusters = splitData[["clusters"]],
                                          checkers = checker,
                                          batchSize = 1L,
                                          cores = 6L,
                                          distance = "cosine",
                                          hclustMethod = "ward.D2",
                                          saveObj = FALSE)

  merge.clusters.test <- mergedData[["clusters"]]
  saveRDS(merge.clusters.test,
          file = file.path(testsDir, "merge.clusters.test.RDS"))

  test.dataset.clusters2 <- merge.clusters.test
  save(test.dataset.clusters2, compress = TRUE,
       file = file.path("data", "test.dataset.clusters2.rda"))

  obj <- addClusterization(objCOTAN = obj,
                           clName = "merge",
                           clusters = mergedData[["clusters"]],
                           coexDF = mergedData[["coex"]])
  # nolint end
}
