
options(parallelly.fork.enable = TRUE)

test_that("Convert COTAN to and from SCE on test dataset", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "artificial",
                               sampleCondition = "test")

  obj <- proceedToCoex(obj, calcCoex = FALSE, cores = 1L, saveObj = FALSE)

  coexTest <- readRDS(test_path("coex.test.RDS"))

  fold15 <- function(m) {
    mm <- cbind(m, m, m, m, m)
    return(cbind(mm, mm, mm))
  }
  coex <- cbind(fold15(coexTest[,  1L:10L]), fold15(coexTest[, 11L:20L]),
                fold15(coexTest[, 21L:30L]), fold15(coexTest[, 31L:40L]))
  colnames(coex) <- rownames(coex)
  coex <- pack(forceSymmetric(coex, uplo = "U"))

  obj@genesCoex <- coex
  obj@metaDataset <- updateMetaInfo(obj@metaDataset,
                                    datasetTags()[["gsync"]], TRUE)

  objSCE <- convertToSingleCellExperiment(objCOTAN = obj)

  genesMeta <- cbind(GenesNames = getGenes(obj), getMetadataGenes(obj))
  cellsMeta <- cbind(CellsIDs   = getCells(obj), getMetadataCells(obj))

  expect_identical(counts(objSCE), getRawData(obj))
  expect_identical(rowData(objSCE), S4Vectors::DataFrame(genesMeta))
  expect_identical(colData(objSCE), S4Vectors::DataFrame(cellsMeta))
  expect_named(S4Vectors::metadata(objSCE),
               c("genesCoex", "cellsCoex", "datasetMetadata"))
  expect_identical(S4Vectors::metadata(objSCE)[["datasetMetadata"]],
                   getMetadataDataset(obj))
  expect_identical(S4Vectors::metadata(objSCE)[["genesCoex"]], coex)
  expect_identical(S4Vectors::metadata(objSCE)[["cellsCoex"]],
                   emptySymmetricMatrix())

  newObj <- convertFromSingleCellExperiment(objSCE = objSCE)

  expect_identical(getRawData(newObj), getRawData(obj))
  expect_identical(getMetadataGenes(newObj), genesMeta)
  expect_identical(getMetadataCells(newObj), cellsMeta)
  expect_identical(getMetadataDataset(newObj), getMetadataDataset(obj))
  expect_identical(getMetadataElement(obj, datasetTags()[["gsync"]]),
                   getMetadataElement(newObj, datasetTags()[["gsync"]]))
  expect_identical(getGenesCoex(newObj), getGenesCoex(obj))

  rawNorm <- readRDS(file.path(getwd(), "raw.norm.test.RDS"))
  lambda <- readRDS(file.path(getwd(), "lambda.test.RDS"))
  dispersion <- readRDS(file.path(getwd(), "dispersion.test.RDS"))
  nu <- readRDS(file.path(getwd(), "nu.test.RDS"))

  genesNamesTest <- readRDS(file.path(getwd(), "genes.names.test.RDS"))
  cellsNamesTest <- readRDS(file.path(getwd(), "cells.names.test.RDS"))

  expect_equal(getNuNormData(newObj)[genesNamesTest, cellsNamesTest],
               rawNorm, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getLambda(newObj)[genesNamesTest],
               lambda, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getNu(newObj)[cellsNamesTest],
               nu, tolerance = 1.0e-14, ignore_attr = FALSE)
  expect_equal(getDispersion(newObj)[genesNamesTest],
               dispersion, tolerance = 1.0e-10, ignore_attr = FALSE)
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

  sce <- Seurat::as.SingleCellExperiment(srat)
  suppressWarnings({
    expect_warning(obj <- convertFromSingleCellExperiment(sce))
  })

  allDims <- set_names(c(600L, 1000L, 0L, 0L, 0L, 0L, 0L, 0L, 6L, 2L),
    c("raw1", "raw2", "genesCoex1", "genesCoex2", "cellsCoex1", "cellsCoex2",
      "metaDataset", "metaGenes", "metaCells", "clustersCoex"))
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

  obj <- proceedToCoex(obj, calcCoex = FALSE, cores = 1L, saveObj = FALSE)

  expect_identical(colnames(getMetadataGenes(obj)),
                   c("feGenes", "lambda", "dispersion"))
  expect_identical(colnames(getMetadataCells(obj)),
                   c("nCount_RNA", "nFeature_RNA", "ident",
                     "CL_RNA_snn_res.0.8", "CL_seurat_clusters",
                     "COND_orig.ident", "feCells", "nu"))
})
