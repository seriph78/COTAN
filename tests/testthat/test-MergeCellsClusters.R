tm = tempdir()
stopifnot(file.exists(tm))

test_that("Merge Cells Clusters", {

  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- clean(obj)

  obj <- estimateDispersionBisection(obj, cores = 12)

  obj <- calculateCoex(obj)

  clusters <- c(readRDS(file.path(getwd(), "clusters1.RDS")), NA, NA)

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  list[coexDF, pValueDF] <- DEAOnClusters(obj)

  obj <- addClusterizationCoex(obj, clName = "clusters", coexDF = coexDF)

  expect_equal(colnames(coexDF),   unique(clusters))
  expect_equal(rownames(coexDF),   getGenes(obj)[flagNotHousekeepingGenes(obj)])

  expect_equal(colnames(pValueDF), unique(clusters))
  expect_equal(rownames(pValueDF), getGenes(obj)[flagNotHousekeepingGenes(obj)])

  coexPiece   <- c(-0.0054300268948555, -0.00654478716198637, -0.0608531124946441,
                   -0.0686569187994827, -0.00670564483533082, -0.0427999576476754,
                    0.0266846613679557, -0.0127148993053832,  -0.0953800808055625,
                   -0.0323097621779541, -0.0450472379280856,  -0.00922266136409998)
  pValuePiece <- c( 0.863898995126085,   0.320016391930115,    0.85767857630923,
                    0.870788335796612,   0.149278055414262,    0.999302205347222,
                    0.909152716370667,   0.950245644366666,    1,
                    0.730822531588433,   0.793250288753521,    0.500103002167896)

  expect_true(all.equal(coexDF  [["0"]][1:12], coexPiece  ))
  expect_true(all.equal(pValueDF[["2"]][1:12], pValuePiece))

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

    temp.obj <- clean(temp.obj)

    temp.obj <- estimateDispersionBisection(temp.obj, cores = 12)
    gc()

    temp.obj <- calculateCoex(temp.obj)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_true( sum(GDI_data[["GDI"]] >= 1.5) <= 0.01 * nrow(GDI_data))
  }
})
