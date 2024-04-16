
# Creates the files to be reloaded by the tests for comparisons
library(zeallot)

outputTestDatasetCreation <- function(testsDir = file.path("tests",
                                                           "testthat")) {
  utils::data("test.dataset", package = "COTAN")
  options(parallelly.fork.enable = TRUE)

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, cores = 6L, saveObj = FALSE)
  #saveRDS(obj, file = file.path(testsDir,"temp.RDS"))

  cells.names.test <- getCells(obj)[c(1L:10L, 591L:610L, 991L:1000L)]
  genes.names.test <- getGenes(obj)[c(1L:10L, 291L:310L, 591L: 600L)]
  saveRDS(cells.names.test, file.path(testsDir, "cells.names.test.RDS"))
  saveRDS(genes.names.test, file.path(testsDir, "genes.names.test.RDS"))

  dispersion.test <- getDispersion(obj)[genes.names.test]
  saveRDS(dispersion.test, file.path(testsDir, "dispersion.test.RDS"))

  raw.norm.test <- getNormalizedData(obj)[genes.names.test, cells.names.test]
  saveRDS(raw.norm.test, file.path(testsDir, "raw.norm.test.RDS"))

  coex.test <- getGenesCoex(obj, genes = genes.names.test, zeroDiagonal = FALSE)
  saveRDS(coex.test, file.path(testsDir, "coex.test.RDS"))

  lambda.test <- getLambda(obj)[genes.names.test]
  saveRDS(lambda.test, file.path(testsDir, "lambda.test.RDS"))

  GDI.test <- calculateGDI(obj)
  GDI.test <- GDI.test[genes.names.test, ]
  saveRDS(GDI.test, file.path(testsDir, "GDI.test.RDS"))

  nu.test <- getNu(obj)[cells.names.test]
  saveRDS(nu.test, file.path(testsDir, "nu.test.RDS"))

  pval.test <- calculatePValue(obj, geneSubsetCol = genes.names.test)
  saveRDS(pval.test, file.path(testsDir, "pval.test.RDS"))

  GDIThreshold <- 1.46
  initialResolution <- 0.8

  clusters <- cellsUniformClustering(obj, GDIThreshold = GDIThreshold,
                                     initialResolution =   initialResolution,
                                     cores = 6L, saveObj = FALSE)[["clusters"]]
  saveRDS(clusters, file.path(testsDir, "clusters1.RDS"))

  coexDF <- DEAOnClusters(obj, clusters = clusters)
  obj <- addClusterization(obj, clName = "clusters",
                           clusters = clusters, coexDF = coexDF)

  saveRDS(coexDF[genes.names.test, ],
          file.path(testsDir, "coex.test.cluster1.RDS"))

  pvalDF <- pValueFromDEA(coexDF, getNumCells(obj), method = "none")

  saveRDS(pvalDF[genes.names.test, ],
          file.path(testsDir, "pval.test.cluster1.RDS"))

  c(mergedClusters, mCoexDF) %<-%
    mergeUniformCellsClusters(objCOTAN = obj,
                              clusters = NULL,
                              GDIThreshold = GDIThreshold,
                              cores = 6L,
                              distance = "cosine",
                              hclustMethod = "ward.D2",
                              saveObj = FALSE)

  saveRDS(mergedClusters[genes.names.test],
          file.path(testsDir, "cluster_data_merged.RDS"))
}
