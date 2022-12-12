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
  skip("merge_cell.clusters does not update the returned obj")

  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- clean(obj, calcExtraData = FALSE)[["objCOTAN"]]

  obj <- estimateDispersionBisection(obj, cores = 12)

  obj <- calculateCoex(obj)

  obj <- addClusterization(obj, clName = "clusters",
                           clusters = readRDS(file.path(getwd(), "clusters1.RDS")))

  obj <- DEA_on_clusters(as(obj,"scCOTAN"))[[1]]

  initial.cluster.number <- dim(obj@cluster_data)[2]
  obj <- merge_cell.clusters(obj = obj,
                             cond = "test",
                             cores = 12,
                             srat = "Seurat_obj_test_with_cotan_clusters.RDS",
                             out_dir = tm,
                             GEO = "test",
                             sc.method = "10X")

  expect_true( dim(obj@cluster_data)[2] < initial.cluster.number)

  #cluster_data <- readRDS(file.path(getwd(),"cluster_data_marged.RDS"))

  #expect_equal(obj@cluster_data[genes.names.test,], cluster_data)

  for (cl in unique(obj@clusters)) {
    cells.to_test <-  names(obj@clusters[obj@clusters == cl])

    temp.obj <- obj@raw[,colnames(obj@raw) %in% cells.to_test]

    temp.obj <- COTAN(raw = temp.obj)
    temp.obj <- initializeMetaDataset(temp.obj,
                                      GEO = "",
                                      sequencingMethod = " ",
                                      sampleCondition = "temp.clustered")

    temp.obj <- clean(temp.obj, calcExtraData = FALSE)[["objCOTAN"]]

    temp.obj <- estimateDispersionBisection(temp.obj, cores = 12)
    gc()

    temp.obj <- calculateCoex(temp.obj)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_false( dim(GDI_data[GDI_data$GDI >= 1.5,])[1]/dim(GDI_data)[1] > 0.01 )
  }
})
