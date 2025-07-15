options(COTAN.TorchWarning = NULL)
options(parallelly.fork.enable = TRUE)


test_that("COTAN getters", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L,  0L, 3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")
  obj <- clean(obj)

  obj <- estimateLambdaLinear(obj)
  obj <- estimateDispersionBisection(obj, cores = 2L)

  # set tag label as legacy value
  row <- getMetaInfoRow(getMetadataDataset(obj), datasetTags()[["gsync"]])
  obj@metaDataset[row, 1L] <- "genes' coex is in sync:"

  suppressWarnings({
    obj <- calculateCoex(obj, actOnCells = FALSE, returnPPFract = TRUE,
                         optimizeForSpeed = TRUE, deviceStr = "cpu")
  })
  expect_no_warning({
    obj <- calculateCoex(obj, actOnCells = TRUE, returnPPFract = TRUE,
                         optimizeForSpeed = FALSE)
  })

  obj <- addClusterization(obj, clName = "Test",
                           clusters = set_names(rep(c(1L, 2L), 10L),
                                                colnames(raw)))
  obj <- addClusterization(obj, clName = "Test2",
                           clusters = set_names(rep(c(2L, 1L), 10L),
                                                colnames(raw)))

  obj <- addCondition(obj, condName = "Test",
                      conditions = set_names(rep(c("F", "M"), 10L),
                                             colnames(raw)))
  obj <- addCondition(obj, condName = "Test", override = TRUE,
                      conditions = set_names(c(rep(c("F", "M"), 9L), "M", "F"),
                                             colnames(raw)))

  metaInfo <- c("V", "10X", "Test", "20", "TRUE", "TRUE",
                paste0(10.0 / 55.0), paste0(0L))

  expect_identical(getRawData(obj), as(as(raw, "dMatrix"), "sparseMatrix"))
  expect_identical(getNumGenes(obj), 10L)
  expect_identical(getNumCells(obj), 20L)
  expect_identical(getGenes(obj), head(LETTERS, getNumGenes(obj)))
  expect_identical(getCells(obj), head(letters, getNumCells(obj)))
  expect_identical(getZeroOneProj(obj), sign(getRawData(obj)))
  expect_identical(getCellsSize(obj), colSums(getRawData(obj)))
  expect_identical(getNumExpressedGenes(obj), colSums(getRawData(obj) != 0L))
  expect_identical(getGenesSize(obj), rowSums(getRawData(obj)))
  expect_identical(getNumOfExpressingCells(obj), rowSums(getRawData(obj) != 0L))
  expect_identical(getNuNormData(obj),
                   t(t(getRawData(obj)) * (1.0 / getNu(obj))))
  expect_identical(getLogNormData(obj),
                   log1p(t(t(getRawData(obj)) * (1.0e4 / getCellsSize(obj)))) /
                     log(10.0))
  expect_equal(getMetadataDataset(obj)[[1L]], datasetTags()[1L:8L],
               ignore_attr = TRUE)
  expect_identical(getMetadataDataset(obj)[[2L]], metaInfo)
  expect_identical(getMetaInfoRow(getMetadataDataset(obj),
                                  "genes' COEX is in sync:"), 5L)
  expect_identical(getMetadataElement(obj, tag = "genes' coex is in sync:"),
                   "TRUE")
  expect_true(isCoexAvailable(obj, actOnCells = FALSE))
  expect_true(isCoexAvailable(obj, actOnCells = FALSE, ignoreSync = TRUE))
  expect_true(isCoexAvailable(obj, actOnCells = TRUE))
  expect_true(isCoexAvailable(obj, actOnCells = TRUE, ignoreSync = TRUE))
  expect_setequal(colnames(getMetadataGenes(obj)),
                  c("lambda", "feGenes", "dispersion"))
  expect_identical(rownames(getMetadataGenes(obj)), getGenes(obj))
  expect_setequal(colnames(getMetadataCells(obj)),
                  c("nu", "feCells", names(getClustersCoex(obj)),
                    getAllConditions(obj, keepPrefix = TRUE)))
  expect_identical(rownames(getMetadataCells(obj)), getCells(obj))
  expect_length(getClustersCoex(obj), 2L)
  expect_named(getClustersCoex(obj), paste0("CL_", getClusterizations(obj)))
  expect_identical(getClustersCoex(obj)[["CL_Test"]], data.frame())
  expect_equal(getNu(obj), getMetadataCells(obj)[["nu"]],
               ignore_attr = TRUE)
  expect_equal(getLambda(obj), getMetadataGenes(obj)[["lambda"]],
               ignore_attr = TRUE)
  expect_equal(getDispersion(obj), getMetadataGenes(obj)[["dispersion"]],
               ignore_attr = TRUE)
  expect_identical(flagNotFullyExpressedGenes(obj),
                   c(FALSE, rep(TRUE, getNumGenes(obj) - 1L)))
  expect_identical(getFullyExpressedGenes(obj), LETTERS[[1L]])
  expect_identical(flagNotFullyExpressingCells(obj),
                   rep(TRUE, getNumCells(obj)))
  expect_identical(getFullyExpressingCells(obj), vector(mode = "character"))
  expect_identical(dim(getGenesCoex(obj)),
                   as.integer(c(getNumGenes(obj), getNumGenes(obj))))
  expect_identical(dim(getCellsCoex(obj)),
                   as.integer(c(getNumCells(obj), getNumCells(obj))))
  expect_identical(getClusterizations(obj), c("Test", "Test2"))
  expect_identical(getClusterizationName(obj), "Test2")
  expect_identical(getClusterizationName(obj, clName = "Test",
                                         keepPrefix = TRUE), "CL_Test")
  expect_setequal(names(getClusterizationData(obj)), c("coex", "clusters"))
  expect_equal(getClusterizationData(obj)[["clusters"]],
               getMetadataCells(obj)[["CL_Test2"]], ignore_attr = TRUE)
  expect_identical(getClusters(obj), getClusterizationData(obj)[["clusters"]])
  expect_identical(getClusterizationData(obj)[["coex"]], data.frame())
  expect_identical(getAllConditions(obj), "Test")
  expect_identical(getConditionName(obj), "")
  expect_identical(getConditionName(obj, condName = "Test",
                                    keepPrefix = TRUE), "COND_Test")
  expect_identical(levels(getCondition(obj)), "NoCond")
  expect_identical(factorToVector(getCondition(obj,
                                               condName = "Test"))[19L:20L],
                   set_names(c("M", "F"), letters[19L:20L]))
  expect_length(getDims(obj), 7L)
  expect_identical(getDims(obj)[["raw"]], dim(getRawData(obj)))
  expect_identical(getDims(obj)[["metaCells"]], ncol(getMetadataCells(obj)))
  expect_null(getDims(obj)[["clusterCoex"]])

  # getting data matrices
  expect_identical(getZeroOneProj(obj), zeroOne)
  expect_identical(getProbabilityOfZero(obj), probZero)

  ## bound probablitities to avoid extremes
  bProbZero <- pmax(1.0e-8, pmin(1.0 - 1.0e-8, probZero))


  getLH <- function(objCOTAN, formula) {
    return(calculateLikelihoodOfObserved(objCOTAN, formula))
  }

  asDataMatrix <- function(objCOTAN, array) {
    matrix(array, nrow = getNumGenes(objCOTAN),
           dimnames = list(getGenes(objCOTAN), getCells(objCOTAN)))
  }

  ## check default
  expect_identical(calculateLikelihoodOfObserved(obj),
                   getLH(obj, formula = "raw"))

  expect_identical(
    getLH(obj, formula = "raw"),
    asDataMatrix(obj, ifelse(zeroOne, (1.0 - bProbZero), bProbZero)))
  expect_identical(
    getLH(obj, formula = "log"),
    asDataMatrix(obj, ifelse(zeroOne, log1p(-bProbZero), log(bProbZero))))
  expect_identical(
    getLH(obj, formula = "der"),
    asDataMatrix(obj,  ifelse(zeroOne, -1.0/(1.0 - bProbZero), 1.0/bProbZero)))
  expect_identical(
    getLH(obj, formula = "sLog"),
    asDataMatrix(obj, ifelse(zeroOne, -log1p(-bProbZero), log(bProbZero))))

  ## strings are case sensitive
  expect_error(getLH(obj, formula = "SLog"))

  # general data retrieval

  ## check default
  expect_identical(getDataMatrix(obj), as.matrix(getLogNormData(obj)))

  expect_identical(getDataMatrix(obj, dataMethod = "RW"),                      as.matrix(getRawData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "Raw"),                     as.matrix(getRawData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "RawData"),                 as.matrix(getRawData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "NN"),                      as.matrix(getNuNormData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "NuNorm"),                  as.matrix(getNuNormData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "Normalized"),              as.matrix(getNuNormData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "LN"),                      as.matrix(getLogNormData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "LogNorm"),                 as.matrix(getLogNormData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "LogNormalized"),           as.matrix(getLogNormData(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "BI"),                      as.matrix(zeroOne))
  expect_identical(getDataMatrix(obj, dataMethod = "Bin"),                     as.matrix(zeroOne))
  expect_identical(getDataMatrix(obj, dataMethod = "Binarized"),               as.matrix(zeroOne))
  expect_identical(getDataMatrix(obj, dataMethod = "BD"),                      as.matrix(    zeroOne + probZero - 1.0))
  expect_identical(getDataMatrix(obj, dataMethod = "BinDiscr"),                as.matrix(    zeroOne + probZero - 1.0))
  expect_identical(getDataMatrix(obj, dataMethod = "BinarizedDiscrepancy"),    as.matrix(    zeroOne + probZero - 1.0))
  expect_identical(getDataMatrix(obj, dataMethod = "AB"),                      as.matrix(abs(zeroOne + probZero - 1.0)))
  expect_identical(getDataMatrix(obj, dataMethod = "AdjBin"),                  as.matrix(abs(zeroOne + probZero - 1.0)))
  expect_identical(getDataMatrix(obj, dataMethod = "AdjBinarized"),            as.matrix(abs(zeroOne + probZero - 1.0)))
  expect_identical(getDataMatrix(obj, dataMethod = "LH"),                      as.matrix(getLH(obj, "raw")))
  expect_identical(getDataMatrix(obj, dataMethod = "Like"),                    as.matrix(getLH(obj, "raw")))
  expect_identical(getDataMatrix(obj, dataMethod = "Likelihood"),              as.matrix(getLH(obj, "raw")))
  expect_identical(getDataMatrix(obj, dataMethod = "LL"),                      as.matrix(getLH(obj, "log")))
  expect_identical(getDataMatrix(obj, dataMethod = "LogLike"),                 as.matrix(getLH(obj, "log")))
  expect_identical(getDataMatrix(obj, dataMethod = "LogLikelihood"),           as.matrix(getLH(obj, "log")))
  expect_identical(getDataMatrix(obj, dataMethod = "DL"),                      as.matrix(getLH(obj, "der")))
  expect_identical(getDataMatrix(obj, dataMethod = "DerLogL"),                 as.matrix(getLH(obj, "der")))
  expect_identical(getDataMatrix(obj, dataMethod = "DerivativeLogLikelihood"), as.matrix(getLH(obj, "der")))
  expect_identical(getDataMatrix(obj, dataMethod = "SL"),                      as.matrix(getLH(obj, "sLog")))
  expect_identical(getDataMatrix(obj, dataMethod = "SignLogL"),                as.matrix(getLH(obj, "sLog")))
  expect_identical(getDataMatrix(obj, dataMethod = "SignedLogLikelihood"),     as.matrix(getLH(obj, "sLog")))

  ## strings are case sensitive
  expect_error(getDataMatrix(obj, dataMethod = "signLogl"))

  expect_message(getSelectedGenes(obj))
  expect_identical(getSelectedGenes(obj),
                   getSelectedGenes(obj, genesSel = "HGDI", numGenes = 2000L))

  expect_identical(getSelectedGenes(obj, genesSel = "HGDI", numGenes = 5L),
                   c("C", "D", "F", "G", "I"))
  expect_identical(suppressWarnings(
                     getSelectedGenes(obj, genesSel = "HVG_Seurat",
                                      numGenes = 5L)[1L:4L]),
                   c("B", "C", "D", "E"))
  expect_error(getSelectedGenes(obj, genesSel = "HVG_Scanpy", numGenes = 5L))
  expect_identical(getSelectedGenes(obj, genesSel = c("C", "A", "D", "E", "B")),
                   LETTERS[1L:5L])
})
