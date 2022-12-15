tm = tempdir()
stopifnot(file.exists(tm))

#test_that("DEA_on_clusters_test", {
#  obj <- readRDS(file.path(tm,"temp.RDS"))
#  temp <- DEA_on_clusters(obj)
#  saveRDS(temp[[1]], file = file.path(tm,"temp.RDS") )

#  pval.cl <- readRDS(file.path(getwd(),"pval.test.cluster1.RDS"))

#  error <- sum((temp[[2]][genes.names.test,] - pval.cl)**2, na.rm = T)
#  if(error > 0.001 ){
#    warning("Error difference grater than 0.001!")
#  }

#  expect_true(error < 10^(-2))
#})


test_that("Merge Cells Clusters", {
  skip("merge_cell.clusters: p_values_clusters_merged.csv not found")

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

  obj <- addClusterizationCoex(obj, clName = "clusters", coexDF = DEA_on_clusters(obj)[[1]])

  initial.cluster.number <- length(unique(clusters))

  obj <- merge_cell.clusters(obj = obj,
                             cond = "test",
                             cores = 12,
                             out_dir = tm,
                             GEO = "test",
                             sc.method = "10X")

  list[clusters, cluster_coex] <- getClusterizationData(obj)

  expect_true( length(unique(clusters)) < initial.cluster.number)

  #cluster_data <- readRDS(file.path(getwd(),"cluster_data_marged.RDS"))

  #expect_equal(obj@cluster_data[genes.names.test,], cluster_data)

  raw <- getRawData(obj)

  for (cl in unique(clusters)) {
    cellsToDrop <- names(clusters[clusters != cl])

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
