tm = tempdir()
stopifnot(file.exists(tm))

test_that("Merge Cells Clusters", {

  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- proceedToCoex(obj, cores = 12, saveObj = FALSE)

  clusters <- c(readRDS(file.path(getwd(), "clusters1.RDS")), NA, NA)

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  list[coexDF, pValueDF] <- DEAOnClusters(obj)

  obj <- addClusterizationCoex(obj, clName = "clusters", coexDF = coexDF)

  expect_equal(colnames(coexDF),   unique(clusters))
  expect_equal(rownames(coexDF),   getGenes(obj)[flagNotHousekeepingGenes(obj)])

  expect_equal(colnames(pValueDF), unique(clusters))
  expect_equal(rownames(pValueDF), getGenes(obj)[flagNotHousekeepingGenes(obj)])

  coexPiece   <- c(-0.0054302803499579, -0.00654552572775876, -0.0608531866528491,
                   -0.0686568746635675, -0.00670647495847967, -0.042799959379387,
                    0.0266844833336902, -0.0127147843243859,  -0.0953799186930419,
                   -0.0323098265214701, -0.0450473968910211,  -0.00922233120011961)
  pValuePiece <- c( 0.863895661692466,   0.320043515678809,    0.857677582786368,
                    0.870787782582908,   0.149310637904463,    0.999302205563701,
                    0.909151014780023,   0.950244976264801,    1,
                    0.730821136021827,   0.793247447396414,    0.500113489568265)

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

    temp.obj <- proceedToCoex(temp.obj, cores = 12, saveObj = FALSE)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_true( sum(GDI_data[["GDI"]] >= 1.5) <= 0.01 * nrow(GDI_data))
  }
})
