
# Creates the files to be reloaded by the tests for comparisons

outputTestDatasetCreation <- function(testsDir = "tests/testthat"){
  utils::data("test.dataset.col", package = "COTAN")
  options(parallelly.fork.enable = TRUE)

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- proceedToCoex(obj, cores = 12, saveObj = FALSE)
  #saveRDS(obj, file = file.path(testsDir,"temp.RDS"))

  cell.names.test  <- getCells(obj)[c(1:10,2474:2484,4990:5000)]
  genes.names.test <- getGenes(obj)[c(1:10,990:1000,1971:1981)]
  saveRDS(cell.names.test, file.path(testsDir, "cell.names.test.RDS"))
  saveRDS(genes.names.test, file.path(testsDir, "genes.names.test.RDS"))

  dispersion.test <- getDispersion(obj)[genes.names.test]
  saveRDS(dispersion.test, file.path(testsDir, "dispersion.test.RDS"))

  raw.norm.test <- getNormalizedData(obj)[genes.names.test, cell.names.test]
  saveRDS(raw.norm.test, file.path(testsDir, "raw.norm.test.RDS"))

  coex.test <- getGenesCoex(obj, genes = genes.names.test, zeroDiagonal = FALSE)
  saveRDS(coex.test, file.path(testsDir, "coex.test.RDS"))

  lambda.test <- getLambda(obj)[genes.names.test]
  saveRDS(lambda.test, file.path(testsDir, "lambda.test.RDS"))

  GDI.test <- calculateGDI(obj)
  GDI.test <- GDI.test[genes.names.test, ]
  saveRDS(GDI.test, file.path(testsDir, "GDI.test.RDS"))

  nu.test <- getNu(obj)[cell.names.test]
  saveRDS(nu.test, file.path(testsDir, "nu.test.RDS"))

  pval.test <- calculatePValue(obj, geneSubsetCol = genes.names.test)
  saveRDS(pval.test, file.path(testsDir, "pval.test.RDS"))

  clusters <- cellsUniformClustering(obj, cores = 12, saveObj = FALSE)
  saveRDS(clusters, file.path(testsDir, "clusters1.RDS"))

  c(coexDF, pvalDF) %<-% DEAOnClusters(obj, clusters = clusters)
  obj <- addClusterization(obj, clName = "clusters",
                           clusters = clusters, coexDF = coexDF)

  saveRDS(coexDF[genes.names.test, ], file.path(testsDir, "coex.test.cluster1.RDS"))
  saveRDS(pvalDF[genes.names.test, ], file.path(testsDir, "pval.test.cluster1.RDS"))

  c(mergedClusters, mCoexDF, mPValueDf) %<-%
    mergeUniformCellsClusters(objCOTAN = obj,
                              clusters = NULL,
                              cores = 12,
                              distance = "cosine",
                              hclustMethod = "ward.D2",
                              saveObj = FALSE)

  saveRDS(mergedClusters[genes.names.test],
          file.path(testsDir, "cluster_data_merged.RDS"))
}
