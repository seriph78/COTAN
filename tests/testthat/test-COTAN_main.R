tm = tempdir()
stopifnot(file.exists(tm))

test_that("Cell Uniform Clustering", {
  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ",
                               sequencingMethod = "10X",
                               sampleCondition = "example")

  obj <- clean(obj, calcExtraData = FALSE)[["objCOTAN"]]

  obj <- estimateDispersionBisection(obj, cores = 12)

  obj <- calculateCoex(obj)

  clusters <- cell_homogeneous_clustering(obj, cond = "test", cores = 12,
                                          out_dir = paste0(tm,"/"))

  obj <- addClusterization(obj, clName = "clusters", clusters = clusters)

  saveRDS(obj, file = file.path(tm, "temp.RDS"))

  #clusters <- readRDS(file.path(getwd(),"clusters1.RDS"))

  #expect_equal(temp@clusters, clusters)
  ####################################

  # Test the low GDI (homogeneity) for each defined clusters

  ####################################
  expect_equal(getClusterizationData(obj)[["clusters"]], clusters)

  for (cl in sample(unique(clusters), size = 5)) {
    cells.to_test <-  clusters == cl

    #temp.obj <- cluster_homogeneity_check(obj = obj,cells = cells.to_test,
    #                                      out_dir = paste0(tm,"/"),
    #                                      cores = cores,
    #                                      code = 12)

    temp.obj <- COTAN(getRawData(obj)[, cells.to_test, drop = FALSE])
    temp.obj <- initializeMetaDataset(temp.obj, GEO = "",
                                      sequencingMethod = " ",
                                      sampleCondition = "temp.clustered")

    temp.obj <- clean(temp.obj, calcExtraData = FALSE)[["objCOTAN"]]

    temp.obj <- estimateDispersionBisection(temp.obj, cores = 12)
    temp.obj <- calculateCoex(temp.obj)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_false( dim(GDI_data[GDI_data$GDI >= 1.5,])[1]/dim(GDI_data)[1] > 0.01 )
  }
})

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

  obj <- readRDS(file.path(tm,"temp.RDS"))
  obj <- DEA_on_clusters(as(obj,"scCOTAN"))[[1]]

  initial.cluster.number <- dim(obj@cluster_data)[2]
  obj <- merge_cell.clusters(obj = obj,
                             cond = "test",
                             cores = 12,
                             srat = "Seurat_obj_test_with_cotan_clusters.RDS",
                             out_dir = paste0(tm,"/"),
                             GEO = "test",
                             sc.method = "10X")

  expect_true( dim(obj@cluster_data)[2] < initial.cluster.number)

  #saveRDS(obj, file = file.path(tm,"temp.RDS") )
  #cluster_data <- readRDS(file.path(getwd(),"cluster_data_marged.RDS"))

  #expect_equal(obj@cluster_data[genes.names.test,], cluster_data)

  for (cl in unique(obj@clusters)) {
    cells.to_test <-  names(obj@clusters[obj@clusters == cl])
    #temp.obj <- cluster_homogeneity_check(obj = obj,cells = cells.to_test,
    #                                     out_dir = paste0(tm,"/"),
    #                                      cores = cores,
    #                                     code = 12)

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
