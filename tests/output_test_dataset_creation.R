
# Data set features estimation done before clustering (so before dropping the cells not belonging
# to any cluster).

dataset.for.test.creation <- function(){
  tm = tempdir()
  utils::data("test.dataset.col", package = "COTAN")

  obj <- COTAN(raw = test.dataset.col)
  obj <- initializeMetaDataset(obj, GEO = " ", sequencingMethod = "10X",
                               sampleCondition = "example")
  #---------------------------------------------------

  obj <- clean(obj, calcExtraData = FALSE)[["objCOTAN"]]
  obj <- estimateDispersion(obj, cores = 10)
  obj <- calculateCoex(obj)
  obj <- as(obj, "scCOTAN")
  saveRDS(obj, file = file.path(tm,"temp.RDS") )

  cell.names.test  <- colnames(obj@raw.norm[c(1:10,990:1000,1971:1981), c(1:10,2474:2484,4990:5000)])
  genes.names.test <- rownames(obj@raw.norm[c(1:10,990:1000,1971:1981), c(1:10,2474:2484,4990:5000)])

  a.test <- obj@a[genes.names.test]
  saveRDS(a.test, "tests/testthat/a.test.RDS")

  raw.norm.test <- getNormalizedData(obj)[genes.names.test,cell.names.test]
  saveRDS(raw.norm.test, "tests/testthat/raw.norm.test.RDS")

  coex.test <- getCoex(obj, asMatrix = TRUE, genes = genes.names.test)
  saveRDS(coex.test, "tests/testthat/coex.test.RDS")

  lambda.test <- obj@lambda[genes.names.test]
  saveRDS(lambda.test, "tests/testthat/lambda.test.RDS")

  GDI.test <- get.GDI(obj)
  GDI.test <- GDI.test[genes.names.test,]
  saveRDS(GDI.test, "tests/testthat/GDI.test.RDS")

  nu.test <- obj@nu[cell.names.test]
  saveRDS(nu.test, "tests/testthat/nu.test.RDS")

  pval.test <- get.pval(obj,gene.set.col = genes.names.test,gene.set.row = genes.names.test)
  saveRDS(pval.test, "tests/testthat/pval.test.RDS")

  obj <- cell_homogeneous_clustering(cond = "test",out_dir = paste0(tm,"/"), in_dir = paste0(tm,"/"), cores = 2,
                                      dataset_type = "COTAN", dataset_name = "temp.RDS",
                                      GEO = "test",sc.method ="10X"
  )
  saveRDS(obj@clusters, "tests/testthat/clusters1.RDS")

  saveRDS(obj, file = file.path(tm,"temp.RDS") )

  temp <- DEA_on_clusters(obj)
  obj <- temp[[1]]
  pval <- temp[[2]]
  saveRDS(pval[genes.names.test,], "tests/testthat/pval.test.cluster1.RDS")
  saveRDS(obj, file = file.path(tm,"temp.RDS") )

  obj <- merge_cell.clusters(obj = obj,cond = "test",cores=10,out_dir_root = paste0(tm,"/"),
                      srat = "Seurat_obj_test_with_cotan_clusters.RDS" ,out_dir = paste0(tm,"/") ,GEO = "test",
                      sc.method = "10X",mt = TRUE, mt_prefix="^MT")

  saveRDS(obj@cluster_data[genes.names.test,], "tests/testthat/cluster_data_marged.RDS")

  saveRDS(cell.names.test,"tests/testthat/cell.names.test.RDS")
  saveRDS(genes.names.test,"tests/testthat/genes.names.test.RDS")
}
