library(Matrix)

setLoggingLevel(3)

tm = tempdir()
stopifnot(file.exists(tm))

test_that("Empty matrices", {
  expect_equal(dim(emptySparseMatrix()), c(0,0))
  expect_s4_class(emptySparseMatrix(), "dgCMatrix")

  expect_equal(dim(emptySymmetricMatrix()), c(0,0))
  expect_s4_class(emptySymmetricMatrix(), "dspMatrix")
})

test_that("'COTAN' constructor", {
  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)

  expect_s4_class(obj, "COTAN")

  expect_equal(obj@raw, as(as.matrix(test.dataset.col), "sparseMatrix"))
})

test_that("'scCOTAN' converters",{
  raw <- matrix(c(1,0,4,2,11,0,6,7,0,9,10,8,0,0,0,3,0,0,2,0), nrow = 10, ncol = 20)
  rownames(raw) = LETTERS[1:10]
  colnames(raw) = letters[1:20]

  tags <- c("GEO:", "scRNAseq method:", "starting n. of cells:", "Condition sample:")

  obj <- COTAN(raw = raw)
  obj <- initializeMetaDataset(obj, GEO = "V", sequencingMethod = "10X", sampleCondition = "Test")
  obj <- clean(obj, calcExtraData = FALSE)[[1]]
  obj <- estimateDispersion(obj)
  obj <- calculateCoex(obj, actOnCells = FALSE, optimizeForSpeed = FALSE)

  coexDF <- as.data.frame(atan(getNormalizedData(obj)[,1:2]-0.5)/pi*2)
  colnames(coexDF) <- c(1, 2)

  obj <- addClusterization(obj, clusterizationName = "clusters",
                           clusters = rep(colnames(coexDF), 10),
                           coexDF = coexDF)

  # coerce 'COTAN' -> 'scCOTAN'
  obj_sc <- as(obj, "scCOTAN")

  expect_equal(obj_sc@raw,      obj@raw)
  expect_equal(obj_sc@raw.norm, getNormalizedData(obj))
  expect_equal(obj_sc@coex,     obj@genesCoex)
  expect_equal(obj_sc@nu,       getNu(obj))
  expect_equal(obj_sc@lambda,   getLambda(obj))
  expect_equal(obj_sc@nu,       getNu(obj))
  expect_equal(obj_sc@a,        getDispersion(obj))
  expect_equal(obj_sc@hk,       getHousekeepingGenes(obj))
  expect_equal(obj_sc@meta,     obj@metaDataset)
  expect_null(obj_sc@yes_yes)
  expect_length(obj_sc@clusters, ncol(obj_sc@raw))
  if (!all(is.na(obj_sc@clusters))) {
    expect_equal(obj_sc@clusters, getClusterizationData(obj)[["clusters"]])
    expect_equal(obj_sc@cluster_data, getClusterizationData(obj)[["coex"]])
  } else {
    expect_length(obj_sc@cluster_data, 0)
  }

  # coerce 'scCOTAN' -> 'COTAN'
  obj2 <- as(obj_sc, "COTAN")

  if (all(flagNotHousekeepingGenes(obj))) {
    #drop the hk column as it won't appear in the obj2 in this case!
    obj@metaGenes[["hkGenes"]] <- NULL
  }

  expect_equal(obj2, obj)
})

gc()
