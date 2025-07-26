
options(parallelly.fork.enable = TRUE)

test_that("Convert COTAN to and from SCE on test dataset", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, calcCoex = FALSE, cores = 6L, saveObj = FALSE)

  coex_test <- readRDS(file.path(getwd(), "coex.test.RDS"))

  fold15 <- function(m) {
    mm <- cbind(m, m, m, m, m)
    return(cbind(mm, mm, mm))
  }
  coex <- cbind(fold15(coex_test[,  1L:10L]), fold15(coex_test[, 11L:20L]),
                fold15(coex_test[, 21L:30L]), fold15(coex_test[, 31L:40L]))
  colnames(coex) <- rownames(coex)
  coex <- pack(forceSymmetric(coex, uplo = "U"))

  obj@genesCoex <- coex
  obj@metaDataset <- updateMetaInfo(obj@metaDataset,
                                    datasetTags()[["gsync"]], TRUE)

  objSCE <- convertToSingleCellExperiment(objCOTAN = obj)

  expect_identical(counts(objSCE), getRawData(obj))
  expect_identical(rowData(objSCE), S4Vectors::DataFrame(getMetadataGenes(obj)))
  expect_identical(colData(objSCE), S4Vectors::DataFrame(getMetadataCells(obj)))
  expect_named(S4Vectors::metadata(objSCE),
               c("genesCoex", "cellsCoex", "datasetMetadata"))
  expect_identical(S4Vectors::metadata(objSCE)[["datasetMetadata"]],
                   getMetadataDataset(obj))
  expect_identical(S4Vectors::metadata(objSCE)[["genesCoex"]], coex)
  expect_identical(S4Vectors::metadata(objSCE)[["cellsCoex"]],
                   emptySymmetricMatrix())

  newObj <- convertFromSingleCellExperiment(objSCE = objSCE)

  expect_identical(getRawData(newObj), getRawData(obj))
  expect_identical(getMetadataGenes(newObj), getMetadataGenes(obj))
  expect_identical(getMetadataCells(newObj), getMetadataCells(obj))
  expect_identical(getMetadataDataset(newObj), getMetadataDataset(obj))
  expect_identical(getMetadataElement(obj, datasetTags()[["gsync"]]),
                   getMetadataElement(newObj, datasetTags()[["gsync"]]))
  expect_identical(getGenesCoex(newObj), getGenesCoex(obj))

  raw.norm <- readRDS(file.path(getwd(), "raw.norm.test.RDS"))
  lambda <- readRDS(file.path(getwd(), "lambda.test.RDS"))
  dispersion <- readRDS(file.path(getwd(), "dispersion.test.RDS"))
  nu <- readRDS(file.path(getwd(), "nu.test.RDS"))

  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))
  cells.names.test <- readRDS(file.path(getwd(), "cells.names.test.RDS"))

  expect_equal(getNuNormData(newObj)[genes.names.test, cells.names.test],
               raw.norm, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getLambda(newObj)[genes.names.test],
               lambda, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getNu(newObj)[cells.names.test],
               nu, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getDispersion(newObj)[genes.names.test],
               dispersion, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_identical(getGenesCoex(newObj, zeroDiagonal = FALSE), coex)
})


test_that("Convert COTAN to and from Seurat via SCE on test dataset", {
  utils::data("test.dataset", package = "COTAN")

  rawData <- as(as.matrix(test.dataset), "dgCMatrix")

  srat <- Seurat::CreateSeuratObject(counts = rawData,
                                     project = "conversion_test")
  srat <- Seurat::NormalizeData(srat)
  srat <- Seurat::FindVariableFeatures(srat,
                                       nfeatures = min(2000L, nrow(rawData)),
                                       selection.method = "vst")
  srat <- Seurat::ScaleData(srat, features = VariableFeatures(object = srat))
  srat <- Seurat::RunPCA(srat, npcs = 25L,
                         features = VariableFeatures(object = srat))
  srat <- Seurat::FindNeighbors(srat, dims = seq_len(25L))

  srat <- Seurat::FindClusters(srat, resolution = 0.8, algorithm = 2L)

  sce <- suppressWarnings(Seurat::as.SingleCellExperiment(srat))
  obj <- convertFromSingleCellExperiment(sce)

  allDims <- set_names(c(600L, 1200L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 6L, 2L),
    c("raw1", "raw2", "numBatches", "genesCoex1", "genesCoex2", "cellsCoex1",
      "cellsCoex2", "metaDataset", "metaGenes", "metaCells", "clustersCoex"))
  expect_identical(unlist(getDims(obj)), allDims)
  expect_identical(getRawData(obj), rawData)
  expect_identical(getClusterizations(obj),
                   c("RNA_snn_res.0.8", "seurat_clusters"))
  expect_identical(getAllConditions(obj), "orig.ident")
  expect_identical(getClusters(obj, clName = "seurat_clusters"),
                   Seurat::Idents(srat))
  expect_identical(getCondition(obj, condName = "orig.ident"),
                   asClusterization(Seurat::FetchData(srat, var = "orig.ident"),
                                    getCells(obj)))

  obj <- proceedToCoex(obj, calcCoex = FALSE, cores = 6L, saveObj = FALSE)

  expect_identical(colnames(getMetadataGenes(obj)),
                   c("feGenes", "lambda", "dispersion"))
  expect_identical(colnames(getMetadataCells(obj)),
                   c("nCount_RNA", "nFeature_RNA", "ident",
                     "CL_RNA_snn_res.0.8", "CL_seurat_clusters",
                     "COND_orig.ident", "feCells", "nu"))
})
