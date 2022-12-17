tm = tempdir()
stopifnot(file.exists(tm))

test_that("Merge Cells Clusters", {

  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- proceedToCoex(obj, cores = 12, saveObj = FALSE)

  clusters <- readRDS(file.path(getwd(), "clusters1.RDS"))
  genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  list[coexDF, pValDF] <- DEAOnClusters(obj)

  obj <- addClusterizationCoex(obj, clName = "clusters", coexDF = coexDF)

  expect_equal(colnames(coexDF), unique(clusters))
  expect_equal(rownames(coexDF), getGenes(obj))

  expect_equal(colnames(pValDF), unique(clusters))
  expect_equal(rownames(pValDF), getGenes(obj))

  coexDF_exp <- readRDS(file.path(getwd(), "coex.test.cluster1.RDS"))
  pValDF_exp <- readRDS(file.path(getwd(), "pval.test.cluster1.RDS"))

  expect_equal(coexDF[genes.names.test, ], coexDF_exp)
  expect_equal(pValDF[genes.names.test, ], pValDF_exp)

  skip("Partial skip: merge_cell.clusters: p_values_clusters_merged.csv not found")

  obj <- merge_cell.clusters(obj = obj,
                             cond = "test",
                             cores = 12,
                             out_dir = tm,
                             GEO = "test",
                             sc.method = "10X")

  list[clusters_merged, clusters_merged_coex] <- getClusterizationData(obj)

  expect_true( length(unique(clusters_merged)) < length(unique(clusters)))

  #cluster_data <- readRDS(file.path(getwd(),"cluster_data_marged.RDS"))

  #expect_equal(obj@cluster_data[genes.names.test,], cluster_data)

  raw <- getRawData(obj)

  for (cl in unique(clusters_merged)) {
    cellsToDrop <- names(clusters[clusters_merged != cl])

    temp.obj <- dropGenesCells(obj, cells = cellsToDrop)

    temp.obj <- proceedToCoex(temp.obj, cores = 12, saveObj = FALSE)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_true( sum(GDI_data[["GDI"]] >= 1.5) <= 0.01 * nrow(GDI_data))
  }
})
