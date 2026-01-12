
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
  obj <- estimateDispersionViaSolver(obj, cores = 2L)

  # add second column to global meta-data
  obj <- addElementToMetaDataset(
    obj, tag = datasetTags()[["sc.m"]], value = c("10X", "V3"))

  # set tag label as legacy value
  row <- getMetaInfoRow(getMetadataDataset(obj), datasetTags()[["gsync"]])
  obj@metaDataset[row, 1L:2L] <- c("genes' coex is in sync:", FALSE)

  suppressWarnings({
    obj <- calculateCoex(obj, actOnCells = FALSE, returnPPFract = TRUE,
                         optimizeForSpeed = TRUE, deviceStr = "cpu")
  })
  expect_no_warning({
    obj <- calculateCoex(obj, actOnCells = TRUE, returnPPFract = TRUE,
                         optimizeForSpeed = FALSE)
  })

  obj <- storeGDI(obj, genesGDI = calculateGDI(obj))

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
                "NegativeBinomial", paste0(10.0 / 55.0), paste0(0L))

  zeroOne <- sign(getRawData(obj))
  probZero <- funProbZeroNegBin(getDispersion(obj),
                                getLambda(obj) %o% getNu(obj))

  expect_identical(getRawData(obj), as(as(raw, "dMatrix"), "sparseMatrix"))
  expect_identical(getNumGenes(obj), 10L)
  expect_identical(getNumCells(obj), 20L)
  expect_identical(getGenes(obj), head(LETTERS, getNumGenes(obj)))
  expect_identical(getCells(obj), head(letters, getNumCells(obj)))
  expect_identical(getCellsSize(obj), colSums(getRawData(obj)))
  expect_identical(getNumExpressedGenes(obj), colSums(getRawData(obj) != 0L))
  expect_identical(getGenesSize(obj), rowSums(getRawData(obj)))
  expect_identical(getNumOfExpressingCells(obj), rowSums(getRawData(obj) != 0L))
  expect_identical(
    getNuNormData(obj),
    t(t(getRawData(obj)) * (1.0 / getNu(obj)))
  )
  expect_identical(
    getLogNormData(obj),
    log1p(t(t(getRawData(obj)) * (1.0e4 / getCellsSize(obj)))) / log(10.0)
  )
  expect_equal(
    getMetadataDataset(obj)[[1L]],
    datasetTags()[c(1L:6L, 9L, 7L, 8L)],
    ignore_attr = TRUE
  )
  expect_identical(getMetadataDataset(obj)[[2L]], metaInfo)
  expect_identical(
    getMetaInfoRow(getMetadataDataset(obj), "genes' COEX is in sync:"), 5L)
  expect_identical(
    as.list(getMetadataElement(obj, tag = "genes' coex is in sync:")),
    list(V2 = "TRUE", V3 = NA_character_)
  )
  expect_true(isCoexAvailable(obj, actOnCells = FALSE))
  expect_true(isCoexAvailable(obj, actOnCells = FALSE, ignoreSync = TRUE))
  expect_true(isCoexAvailable(obj, actOnCells = TRUE))
  expect_true(isCoexAvailable(obj, actOnCells = TRUE, ignoreSync = TRUE))

  expect_setequal(colnames(getMetadataGenes(obj)),
                  c("lambda", "feGenes", "dispersion", "GDI"))
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
  expect_identical(getGDI(obj), getColumnFromDF(calculateGDI(obj), "GDI"))
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

  binDiscr <- function(objCOTAN) {
    rawData <- getRawData(objCOTAN)
    probZero <- getProbabilityOfZero(objCOTAN)
    return(asDataMatrix(objCOTAN,
                        ifelse(rawData != 0L, probZero, probZero - 1.0)))
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
  expect_identical(getDataMatrix(obj, dataMethod = "BD"),                      as.matrix(    binDiscr(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "BinDiscr"),                as.matrix(    binDiscr(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "BinarizedDiscrepancy"),    as.matrix(    binDiscr(obj)))
  expect_identical(getDataMatrix(obj, dataMethod = "AB"),                      as.matrix(abs(binDiscr(obj))))
  expect_identical(getDataMatrix(obj, dataMethod = "AdjBin"),                  as.matrix(abs(binDiscr(obj))))
  expect_identical(getDataMatrix(obj, dataMethod = "AdjBinarized"),            as.matrix(abs(binDiscr(obj))))
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
  suppressWarnings(expect_warning(
      getSelectedGenes(obj, genesSel = "HVG_Seurat", numGenes = 5L)
  ))
  expect_error(getSelectedGenes(obj, genesSel = "HVG_Scanpy", numGenes = 5L))
  expect_identical(getSelectedGenes(obj, genesSel = c("C", "A", "D", "E", "B")),
                   LETTERS[1L:5L])

  calcRDM <- function(objCOTAN, useCoexEigen, dataMethod,
                      numComp, genesSel = "", numGenes = 2000L) {
    return(suppressWarnings(
      calculateReducedDataMatrix(objCOTAN = objCOTAN,
                                 useCoexEigen = useCoexEigen,
                                 dataMethod = dataMethod,
                                 numComp = numComp,
                                 genesSel = genesSel,
                                 numGenes = numGenes)))
  }

  expect_equal(suppressWarnings(abs(calculateReducedDataMatrix(obj))),
               abs(calcRDM(obj, useCoexEigen = FALSE,
                           dataMethod = "LogNormalized", numComp = 50L,
                           genesSel = "HGDI", numGenes = 2000L)),
               tolerance = 1.0e-12)

  m1 <- as.matrix(cbind(rep(2.756809750418045, times = 20L),
                        rep(0.0, times = 20L), rep(0.0, times = 20L),
                        rep(0.0, times = 20L), rep(0.0, times = 20L)))
  colnames(m1) <- paste0("PC", c(1L:5L))
  rownames(m1) <- letters[1L:20L]

  expect_equal(abs(calcRDM(obj, useCoexEigen = FALSE,
                           dataMethod = "AdjBinarized", numComp = 5L,
                           genesSel = "HGDI", numGenes = 8L)),
               m1, tolerance = 1.0e-12)

  expect_equal(abs(calcRDM(obj, useCoexEigen = FALSE,
                           dataMethod = "LogLikelihood", numComp = 5L,
                           genesSel = "HVG_Seurat", numGenes = 8L)),
               m1, tolerance = 1.0e-12)

  m2 <- cbind(
    rep(c(1.991509612121359,     -1.991510865753923    ), times = 10L),
    rep(c(4.556799230950100e-05, -4.815389335445746e-05), times = 10L),
    rep(c(5.701730287707164e-05, -8.375641609457851e-05), times = 10L))
  colnames(m2) <- paste0("EC_", c(1L:3L))
  rownames(m2) <- letters[1L:20L]

  expect_equal(calcRDM(obj, useCoexEigen = TRUE,
                       dataMethod = "BinDiscr", numComp = 5L)[, 1L:3L],
               m2, tolerance = 5.0e-5)

  m3 <- cbind(
    rep(c(-1.96537974590724,    1.96547022027574     ), times = 10L),
    rep(c(-5.14674535931154e-4, 2.4124745673443485e-2), times = 10L),
    rep(c( 0.227956910102598,   0.226238243931497    ), times = 10L))
  colnames(m3) <- paste0("EC_", c(1L:3L))
  rownames(m3) <- letters[1L:20L]
  if (FALSE) {
    attr(m3, "scaled:scale") <-
      set_names(rep(c(3.020580325047483, 3.020462711979378), times = 10L),
                letters[1L:20L])
  }
  expect_equal(calcRDM(obj, useCoexEigen = TRUE,
                       dataMethod = "DerLogL", numComp = 5L)[, 1L:3L],
               m3, tolerance = 5.0e-4)
})

