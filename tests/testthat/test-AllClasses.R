library(Matrix)

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

#TODO: re-enable once the following file has been produced!
if (FALSE) {
test_that("'scCOTAN' converters",{
  obj <- readRDS(file.path(tm, "Completely_analysed_COTAN_obj.RDS"))
  #obj <- readRDS(file.path(tm, "temp.RDS"))

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
}

gc()
