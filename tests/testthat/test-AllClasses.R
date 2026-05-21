library(Matrix)

setLoggingLevel(4L)

options(COTAN.TorchWarning = NULL)

test_that("Empty matrices", {
  expect_identical(dim(emptySparseMatrix()), c(0L, 0L))
  expect_s4_class(emptySparseMatrix(), "dgCMatrix")

  expect_identical(dim(emptySymmetricMatrix()), c(0L, 0L))
  expect_s4_class(emptySymmetricMatrix(), "dspMatrix")
})

test_that("'COTAN' constructor", {
  utils::data("test.dataset", package = "COTAN")

  obj <- COTAN(raw = test.dataset)

  expect_s4_class(obj, "COTAN")

  expect_identical(obj@raw, as(as.matrix(test.dataset), "sparseMatrix"))
})

test_that("'scCOTAN' converters", {
  raw <- matrix(c(1L,  0L, 4L, 2L, 11L, 0L, 6L, 7L, 0L, 9L,
                  10L, 8L, 0L, 0L, 0L,  3L, 0L, 0L, 2L, 0L),
                nrow = 10L, ncol = 20L)
  rownames(raw) <- LETTERS[1L:10L]
  colnames(raw) <- letters[1L:20L]

  tags <- c("GEO:", "scRNAseq method:",
            "starting n. of cells:", "Condition sample:")

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X",
                               sampleCondition = "Test")
  obj <- proceedToCoex(
    obj,
    executionOptions = ExecutionOptions(cores = 1L, optimizeForSpeed = FALSE)
  )

  coexDF <-
    set_names(as.data.frame(atan(getNuNormData(obj)[, 1L:2L] - 0.5) / pi * 2.0),
              c(1L, 2L))

  obj <- addClusterization(obj, clName = "clusters",
                           clusters = set_names(rep(colnames(coexDF), 10L),
                                                colnames(raw)),
                           coexDF = coexDF)

  # coerce 'COTAN' -> 'scCOTAN'
  objSc <- as(obj, "scCOTAN")

  expect_identical(objSc@raw,      obj@raw)
  expect_identical(objSc@raw.norm, getNuNormData(obj))
  expect_identical(objSc@coex,     obj@genesCoex)
  expect_identical(objSc@nu,       getNu(obj))
  expect_identical(objSc@lambda,   getLambda(obj))
  expect_identical(objSc@nu,       getNu(obj))
  expect_identical(objSc@a,        getDispersion(obj))
  expect_identical(objSc@hk,       getFullyExpressedGenes(obj))
  expect_identical(objSc@meta,     obj@metaDataset)
  expect_null(objSc@yes_yes)
  expect_length(objSc@clusters, ncol(objSc@raw))
  if (!all(is.na(objSc@clusters))) {
    expect_identical(objSc@clusters, factorToVector(getClusters(obj)))
    expect_identical(objSc@cluster_data, getClusterizationData(obj)[["coex"]])
  } else {
    expect_length(objSc@cluster_data, 0L)
  }

  # coerce 'scCOTAN' -> 'COTAN'
  obj2 <- as(objSc, "COTAN")

  if (all(flagNotFullyExpressedGenes(obj))) {
    #drop the fe column as it won't appear in the obj2 in this case!
    obj@metaGenes[["feGenes"]] <- NULL
  }
  if (all(flagNotFullyExpressingCells(obj))) {
    #drop the fe column as it won't appear in the obj2 anyway!
    obj@metaCells[["feCells"]] <- NULL
  }

  expect_identical(obj2, obj)
})


test_that("ExecutionOptions stores defaults and explicit values", {
  exec_default <- ExecutionOptions()

  expect_s4_class(exec_default, "ExecutionOptions")
  expect_identical(exec_default@cores, 1L)
  expect_identical(exec_default@optimizeForSpeed, TRUE)
  expect_identical(exec_default@deviceStr, "cuda")
  expect_identical(exec_default@chunkSize, 1024L)

  exec <- ExecutionOptions(
    cores = 4.0,
    optimizeForSpeed = FALSE,
    deviceStr = "cpu",
    chunkSize = 256.0
  )

  expect_s4_class(exec, "ExecutionOptions")
  expect_type(exec@cores, "integer")
  expect_identical(exec@cores, 4L)
  expect_type(exec@optimizeForSpeed, "logical")
  expect_identical(exec@optimizeForSpeed, FALSE)
  expect_type(exec@deviceStr, "character")
  expect_identical(exec@deviceStr, "cpu")
  expect_type(exec@chunkSize, "integer")
  expect_identical(exec@chunkSize, 256L)
})

test_that("ExecutionOptions rejects invalid values", {
  expect_error(ExecutionOptions(cores = -1L))
  expect_error(ExecutionOptions(cores = NA_integer_))

  expect_error(ExecutionOptions(optimizeForSpeed = NA))
  expect_error(ExecutionOptions(optimizeForSpeed = c(TRUE, FALSE)))

  expect_error(ExecutionOptions(deviceStr = NA_character_))
  expect_error(ExecutionOptions(deviceStr = c("cpu", "cuda")))

  expect_error(ExecutionOptions(chunkSize = 0L))
  expect_error(ExecutionOptions(chunkSize = -1L))
  expect_error(ExecutionOptions(chunkSize = NA_integer_))
})

test_that("ReductionOptions stores defaults and explicit values", {
  red_default <- ReductionOptions()

  expect_s4_class(red_default, "ReductionOptions")
  expect_identical(red_default@useCoexEigen, TRUE)
  expect_identical(red_default@dataMethod, "LogLikelihood")
  expect_identical(red_default@numComp, 25L)
  expect_identical(red_default@genesSel, "HGDI")
  expect_identical(red_default@numGenes, 2000L)

  red <- ReductionOptions(
    useCoexEigen = FALSE,
    dataMethod = "LogNormalized",
    numComp = 12.0,
    genesSel = c("GeneA", "GeneB"),
    numGenes = 500.0
  )

  expect_s4_class(red, "ReductionOptions")
  expect_type(red@useCoexEigen, "logical")
  expect_identical(red@useCoexEigen, FALSE)
  expect_identical(red@dataMethod, "LogNormalized")
  expect_type(red@numComp, "integer")
  expect_identical(red@numComp, 12L)
  expect_identical(red@genesSel, c("GeneA", "GeneB"))
  expect_type(red@numGenes, "integer")
  expect_identical(red@numGenes, 500L)
})

test_that("ReductionOptions rejects invalid values", {
  expect_error(ReductionOptions(useCoexEigen = NA))
  expect_error(ReductionOptions(useCoexEigen = c(TRUE, FALSE)))

  expect_error(ReductionOptions(dataMethod = NA_character_))
  expect_error(ReductionOptions(dataMethod = c("Raw", "LogNormalized")))

  expect_error(ReductionOptions(numComp = 0L))
  expect_error(ReductionOptions(numComp = -1L))
  expect_error(ReductionOptions(numComp = NA_integer_))

  expect_error(ReductionOptions(genesSel = character()))
  expect_error(ReductionOptions(genesSel = NA_character_))
  expect_error(ReductionOptions(genesSel = c("GeneA", "")))

  expect_error(ReductionOptions(useCoexEigen = FALSE, numGenes = 0L))
  expect_error(ReductionOptions(useCoexEigen = FALSE, numGenes = NA_integer_))
})


test_that("ClusterDistanceOptions stores defaults and explicit values", {
  dist_default <- ClusterDistanceOptions()

  expect_s4_class(dist_default, "ClusterDistanceOptions")
  expect_identical(dist_default@useDEA, TRUE)
  expect_identical(dist_default@distance, "")

  dist_opt <- ClusterDistanceOptions(
    useDEA = FALSE,
    distance = "euclidean"
  )

  expect_s4_class(dist_opt, "ClusterDistanceOptions")
  expect_type(dist_opt@useDEA, "logical")
  expect_identical(dist_opt@useDEA, FALSE)
  expect_type(dist_opt@distance, "character")
  expect_identical(dist_opt@distance, "euclidean")
})

test_that("ClusterDistanceOptions rejects invalid values", {
  expect_error(ClusterDistanceOptions(useDEA = NA))
  expect_error(ClusterDistanceOptions(useDEA = c(TRUE, FALSE)))

  expect_error(ClusterDistanceOptions(distance = NA_character_))
  expect_error(ClusterDistanceOptions(distance = c("cosine", "euclidean")))
})


gc()
