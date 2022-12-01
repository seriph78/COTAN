tm = tempdir()
stopifnot(file.exists(tm))

#root = "tests/testthat/"
root = ""
genes.names.test <- readRDS(file.path(getwd(), "genes.names.test.RDS"))
cell.names.test <- readRDS(file.path(getwd(), "cell.names.test.RDS"))

test_that("1_initialization", {
    utils::data("test.dataset.col", package = "COTAN")

    rownames(test.dataset.col) <- test.dataset.col$V1
    test.dataset.col <- test.dataset.col[,2:ncol(test.dataset.col)]

    obj <- COTAN(raw = test.dataset.col)
    obj <- initializeMetaDataset(obj, GEO = "V",
                                 sequencingMethod = "10X",
                                 sampleCondition = "example")

    expect_s4_class(obj, "COTAN")
})

test_that("2_cleaning", {
    utils::data("test.dataset.col", package = "COTAN")

    obj <- COTAN(raw = test.dataset.col)
    obj <- initializeMetaDataset(obj, GEO = " ",
                                 sequencingMethod = "10X",
                                 sampleCondition = "example")

    obj <- clean(obj, calcExtraData = FALSE)[["objCOTAN"]]

    stopifnot(file.exists(tm))
    saveRDS(obj, file = file.path(tm,"temp.RDS"))

    raw.norm <- readRDS(file.path(getwd(), "raw.norm.test.RDS"))
    nu <- readRDS(file.path(getwd(), "nu.test.RDS"))

    expect_equal(getNormalizedData(obj)[genes.names.test,cell.names.test],
                 raw.norm)
    expect_equal(getNu(obj)[cell.names.test], nu)
})

test_that("mat_division", {
    utils::data("test.dataset.col", package = "COTAN")

    nu <- readRDS(file.path(getwd(), "nu.test.RDS"))
    raw.norm <- readRDS(file.path(getwd(), "raw.norm.test.RDS"))

    expect_equal( as.matrix(t(t(test.dataset.col[genes.names.test,cell.names.test]) * (1/nu))),
                  as.matrix(raw.norm) )
})


test_that("3_cotan_analysis_test", {
    obj <- readRDS(file.path(tm, "temp.RDS"))

    obj <- estimateDispersion(obj, cores = 12)

    saveRDS(obj, file = file.path(tm, "temp.RDS"))

    dispersion <- readRDS(file.path(getwd(), "a.test.RDS"))
    nu <- readRDS(file.path(getwd(), "nu.test.RDS"))
    lambda <- readRDS(file.path(getwd(), "lambda.test.RDS"))

    expect_equal(getDispersion(obj)[genes.names.test], dispersion)
    expect_equal(getNu(obj)[cell.names.test] , nu)
    expect_equal(getLambda(obj)[genes.names.test], lambda)
})


test_that("4_cotan_coex_test", {
    obj <- readRDS(file.path(tm, "temp.RDS"))

    obj <- calculateCoex(obj)

    coex <- getGenesCoex(obj, genes = genes.names.test)

    saveRDS(obj, file = file.path(tm, "temp.RDS"))

    coex_test <- readRDS(file.path(getwd(), "coex.test.RDS"))

    expect_equal(as.matrix(coex), coex_test)
})


test_that("PCA_test", {

  utils::data("raw.dataset", package = "COTAN")
  pca.tb = readRDS(file.path(getwd(),"pca.tb.RDS"))

  if (nrow(raw.dataset) != nrow(pca.tb) &&
      rownames(raw.dataset)[[1]] == "HK") {
    raw.dataset <- raw.dataset[2:nrow(raw.dataset),]
  }

  pca <- irlba::prcomp_irlba(raw.dataset, n=5)

  pca.raw <- pca$x
  rm(pca)

  rownames(pca.raw) <- rownames(raw.dataset)
  colnames(pca.raw) <- paste0("PC_", c(1:5))

  correlation1.value <- cor(pca.raw[,1], pca.tb[,1])
  correlation2.value <- cor(pca.raw[,2], pca.tb[,2])
  pca.tb[,1] <- correlation1.value*pca.tb[,1]
  pca.tb[,2] <- correlation2.value*pca.tb[,2]
  x1 <- pca.raw[rownames(pca.tb),1] - pca.tb[,1]
  x2 <- pca.raw[rownames(pca.tb),2] - pca.tb[,2]

  dist1 <- sqrt(sum(x1^2))
  dist2 <- sqrt(sum(x2^2))

  expect_true(dist1 < 10^(-4))
  expect_true(dist2 < 10^(-4))
})


test_that("5_calculatePValue_test", {
    obj <- readRDS(file.path(tm,"temp.RDS"))
    pval <- calculatePValue(obj, geneSubsetCol = genes.names.test)
    pval.exp  <- readRDS(file.path(getwd(), "pval.test.RDS"))

    #expect_equal(pval, pval.exp)
    error <- sqrt(mean((log(pval+10^(-10)) - log(pval.exp+10^(-10)))^2))
    error_max <- max(abs(log(pval+10^(-10)) - log(pval.exp+10^(-10))))

    if(error > 0.001 ){
        warning("Error difference grater than 0.001!")
    }

    expect_true(error < 10^(-2))
    expect_true(error_max < 10^(-1))
})


test_that("6_calculateGDI_test", {
  #utils::data("ERCC.cotan", package = "COTAN")
  object <- readRDS(file.path(tm,"temp.RDS"))

  GDI <- calculateGDI(object)[genes.names.test,]
  expect_equal(GDI, readRDS(file.path(getwd(),"GDI.test.RDS")))
})


test_that("vec2mat_rfast_test1", {
    mat <- matrix(0, nrow = 5, ncol = 5)
    mat <- Rfast::lower_tri.assign(mat,c(1:15),diag = T)
    mat <- Rfast::upper_tri.assign(mat,v = Rfast::upper_tri(Rfast::transpose(mat)))

    colnames(mat) <- paste0("row.",c(1:5))
    rownames(mat) <- paste0("row.",c(1:5))

    v <- mat2vec_rfast(mat)

    expect_equal(mat, vec2mat_rfast(v))
})

test_that("vec2mat_rfast_test2", {
    mat <- matrix(0,nrow = 10, ncol = 10)
    mat <- Rfast::lower_tri.assign(mat, c(1:55), diag = T)
    mat <- Rfast::upper_tri.assign(mat, v = Rfast::upper_tri(Rfast::transpose(mat)))

    colnames(mat) <- paste0("row.", c(1:10))
    rownames(mat) <- paste0("row.", c(1:10))

    v <- mat2vec_rfast(mat)

    genes <- c("row.1", "row.2", "row.9", "row.10")

    expect_equal(mat[, genes], vec2mat_rfast(v, genes = genes))
})

test_that("mat2vec_rfast_test", {
    vec <- c(1:15)
    names.v <- paste0("raw", c(1:5))

    names.v
    vec <- list("genes" = names.v, "values" = vec)

    m <- vec2mat_rfast(vec)

    expect_equal(vec, mat2vec_rfast(m))
})


test_that("7_cell_homogeneous_clustering", {
  #obj <- readRDS(file.path(tm,"temp.RDS"))
  temp <- cell_homogeneous_clustering(cond = "test",out_dir = paste0(tm,"/"), in_dir = paste0(tm,"/"),
                                      cores = 12,
                                      dataset_name = "temp.RDS",
                                      GEO = "test", sc.method ="10X"
                                      )
  saveRDS(temp, file = file.path(tm, "temp.RDS"))
  #clusters <- readRDS(file.path(getwd(),"clusters1.RDS"))

  #expect_equal(temp@clusters, clusters)
  ####################################

  # Test the low GDI (homogeneity) for each defined clusters

  ####################################

  for (cl in sample(unique(temp@clusters),size = 5)) {
    cells.to_test <-  names(temp@clusters[temp@clusters == cl])
    #temp.obj <- cluster_homogeneity_check(obj = obj,cells = cells.to_test,
    #                                      out_dir = paste0(tm,"/"),
    #                                      cores = cores,
    #                                      code = 12)

    temp.obj <- temp@raw[,colnames(temp@raw) %in% cells.to_test]

    temp.obj <- COTAN(raw = temp.obj)
    temp.obj <- initializeMetaDataset(temp.obj, GEO = "",
                                      sequencingMethod = " ",
                                      sampleCondition = "temp.clustered")

    temp.obj <- clean(temp.obj, calcExtraData = FALSE)[["objCOTAN"]]

    temp.obj <- estimateDispersion(temp.obj, cores = 12)
    gc()

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




test_that("8_merge_cell.clusters.test", {
  temp <- readRDS(file.path(tm,"temp.RDS"))
  temp <- DEA_on_clusters(as(temp,"scCOTAN"))[[1]]

  initial.cluster.number <- dim(temp@cluster_data)[2]
  temp <- merge_cell.clusters(obj = temp,
                              cond = "test",
                              cores=12,
                              #out_dir_root = paste0(tm,"/"),
                              srat = "Seurat_obj_test_with_cotan_clusters.RDS" ,
                              out_dir = paste0(tm,"/") ,
                              GEO = "test",
                              sc.method = "10X")#,mt = FALSE, mt_prefix="^MT")

  final.cluster.number <- dim(temp@cluster_data)[2]
  expect_true(final.cluster.number < initial.cluster.number)
  #saveRDS(temp, file = file.path(tm,"temp.RDS") )
  #cluster_data <- readRDS(file.path(getwd(),"cluster_data_marged.RDS"))

  #expect_equal(obj@cluster_data[genes.names.test,], cluster_data)

  for (cl in unique(temp@clusters)) {
    cells.to_test <-  names(temp@clusters[temp@clusters == cl])
    #temp.obj <- cluster_homogeneity_check(obj = obj,cells = cells.to_test,
    #                                     out_dir = paste0(tm,"/"),
    #                                      cores = cores,
    #                                     code = 12)

    temp.obj <- temp@raw[,colnames(temp@raw) %in% cells.to_test]

    temp.obj <- COTAN(raw = temp.obj)
    temp.obj <- initializeMetaDataset(temp.obj,
                                      GEO = "",
                                      sequencingMethod = " ",
                                      sampleCondition = "temp.clustered")

    temp.obj <- clean(temp.obj, calcExtraData = FALSE)[["objCOTAN"]]

    temp.obj <- estimateDispersion(temp.obj, cores = 12)
    gc()

    temp.obj <- calculateCoex(temp.obj)
    gc()

    GDI_data <- calculateGDI(temp.obj)

    expect_false( dim(GDI_data[GDI_data$GDI >= 1.5,])[1]/dim(GDI_data)[1] > 0.01 )

  }



})





