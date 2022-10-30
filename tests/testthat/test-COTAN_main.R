tm = tempdir()
stopifnot(file.exists(tm))

#root = "tests/testthat/"
root = ""
genes.names.test <- readRDS(file.path(getwd(),"genes.names.test.RDS"))
cell.names.test <- readRDS(file.path(getwd(),"cell.names.test.RDS"))

test_that("1.initialization", {
    
    utils::data("test.dataset.col", package = "COTAN")
    rownames(test.dataset.col) <- test.dataset.col$V1
    test.dataset.col <- test.dataset.col[,2:ncol(test.dataset.col)]
    obj.temp <- new("scCOTAN",raw = test.dataset.col)
    obj.temp <- initRaw(object = obj.temp,GEO="V" ,sc.method="10X",cond = "example")
    expect_s4_class(obj.temp,"scCOTAN")

})

test_that("2.cleaning", {
    utils::data("test.dataset.col", package = "COTAN")

    obj.temp <- new("scCOTAN",raw = test.dataset.col)
    obj.temp <- initRaw(object = obj.temp,GEO=" " ,sc.method="10X",cond = "example")
    #---------------------------------------------------

    ttm <- clean(obj.temp)
    stopifnot(file.exists(tm))
    saveRDS(ttm$object, file = file.path(tm,"temp.RDS") )
    raw.norm <- readRDS(file.path(getwd(),"raw.norm.test.RDS"))
    nu <- readRDS(file.path(getwd(),"nu.test.RDS"))
    expect_equal(as.matrix(ttm$object@raw.norm[genes.names.test,cell.names.test]),as.matrix(raw.norm))
    expect_equal(ttm$object@nu[cell.names.test],nu)

})

test_that("mat_division", {
    utils::data("test.dataset.col", package = "COTAN")
    nu <- readRDS(file.path(getwd(),"nu.test.RDS"))
    raw.norm <- readRDS(file.path(getwd(),"raw.norm.test.RDS"))

    expect_equal( as.matrix(t(t(test.dataset.col[genes.names.test,cell.names.test]) * (1/nu))),as.matrix(raw.norm) )
})


test_that("3.cotan_analysis_test", {

    obj <- readRDS(file.path(tm,"temp.RDS"))
    obj <- cotan_analysis(obj,cores = 12)

    a <- readRDS(file.path(getwd(),"a.test.RDS"))
    nu <- readRDS(file.path(getwd(),"nu.test.RDS"))
    lambda <- readRDS(file.path(getwd(),"lambda.test.RDS"))
    saveRDS(obj, file = file.path(tm,"temp.RDS") )
    expect_equal(obj@a[genes.names.test], a)
    expect_equal(obj@nu[cell.names.test] , nu)
    expect_equal(obj@lambda[genes.names.test], lambda)
})


test_that("4.cotan_coex_test", {
    obj <- readRDS(file.path(tm,"temp.RDS"))
    obj <- get.coex(obj)

    coex_test <- readRDS(file.path(getwd(),"coex.test.RDS"))
    
    coex <- extract.coex(object = obj,genes = genes.names.test)
    #coex <- vec2mat_rfast(obj@coex,genes = coex_test$genes)

    #coex <- mat2vec_rfast(coex[coex_test$genes,coex_test$genes])

    saveRDS(obj, file = file.path(tm,"temp.RDS") )

    #error <- sqrt(mean((coex$values - coex_test$values)^2))

    #if(error < 0.001 & error > 10^(-4) ){
     #   warning("Error difference grater than 0.0001!")
    #}

    expect_equal(coex, coex_test)

    #expect_true(error < 10^(-3))
})

test_that("PCA_test", {
  
  utils::data("raw.dataset", package = "COTAN")
  pca <- irlba::prcomp_irlba(raw.dataset, n=5)
  
  pca.raw <- pca$x
  rm(pca) 
  
  rownames(pca.raw) <- rownames(raw.dataset)
  colnames(pca.raw) <- paste0("PC_", c(1:5))
  pca.tb = readRDS(file.path(getwd(),"pca.tb.RDS"))
  
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



test_that("5.get_pval_test", {
    object <- readRDS(file.path(tm,"temp.RDS"))
    pval <- get.pval(object, gene.set.col = genes.names.test,gene.set.row = genes.names.test)
    pval.exp  <- readRDS(file.path(getwd(),"pval.test.RDS"))

    #expect_equal(pval, pval.exp)
    error <- sqrt(mean((log(pval+10^(-10)) - log(pval.exp+10^(-10)))^2))
    error_max <- max(abs(log(pval+10^(-10)) - log(pval.exp+10^(-10))))

    if(error > 0.001 ){
        warning("Error difference grater than 0.001!")
    }

    expect_true(error < 10^(-2))
    expect_true(error_max < 10^(-1))
})


test_that("get_GDI_test", {
  object <- readRDS(file.path(tm,"temp.RDS"))
    #utils::data("ERCC.cotan", package = "COTAN")
    GDI <- get.GDI(object)[genes.names.test,]
    expect_equal(GDI, readRDS(file.path(getwd(),"GDI.test.RDS")))

})


test_that("vec2mat_rfast_test1", {
    mat <- matrix(0,nrow = 5, ncol = 5)
    mat <- Rfast::lower_tri.assign(mat,c(1:15),diag = T)
    mat <- Rfast::upper_tri.assign(mat,v = Rfast::upper_tri(Rfast::transpose(mat)))

    colnames(mat) <- paste0("row.",c(1:5))
    rownames(mat) <- paste0("row.",c(1:5))

    v <- mat2vec_rfast(mat)

    expect_equal(mat, vec2mat_rfast(v,genes = "all"))

})

test_that("vec2mat_rfast_test2", {
    mat <- matrix(0,nrow = 10, ncol = 10)
    mat <- Rfast::lower_tri.assign(mat,c(1:55),diag = T)
    mat <- Rfast::upper_tri.assign(mat,v = Rfast::upper_tri(Rfast::transpose(mat)))

    colnames(mat) <- paste0("row.",c(1:10))
    rownames(mat) <- paste0("row.",c(1:10))

    v <- mat2vec_rfast(mat)

    genes <- c("row.1","row.2","row.9","row.10")

    expect_equal(mat[,genes], vec2mat_rfast(v,genes = genes))

})

test_that("mat2vec_rfast_test", {
    vec <- c(1:15)
    names.v <- paste0("raw",c(1:5))

    names.v
    vec <- list("genes"=names.v, "values"=vec)

    m <- vec2mat_rfast(vec)

    expect_equal(vec, mat2vec_rfast(m))

})


test_that("cell_homogeneous_clustering", {
  #obj <- readRDS(file.path(tm,"temp.RDS"))
  temp <- cell_homogeneous_clustering(cond = "test",out_dir = paste0(tm,"/"), in_dir = paste0(tm,"/"), 
                                      cores = 12, 
                                       dataset_name = "temp.RDS", 
                                      GEO = "test",sc.method ="10X"
                                      )
  saveRDS(temp, file = file.path(tm,"temp.RDS") )
  #clusters <- readRDS(file.path(getwd(),"clusters1.RDS"))
  
  #expect_equal(temp@clusters, clusters)
  ####################################
  
  # Test the low GDI (homogeneity) for each defined clusters 
  
  ####################################
  
  for (cl in sample(unique(temp@clusters),size = 5)) {
    cells.to_test <-  names(temp@clusters[temp@clusters == cl])
    #temp.obj <- cluster_homogeneity_check(obj = obj,cells = cells.to_test,
     #                                     out_dir = paste0(tm,"/"),
    #                                      cores = cores,
     #                                     code = 12)
    
    temp.obj <- temp@raw[,colnames(temp@raw) %in% cells.to_test]
    
    temp.obj <- new("scCOTAN",raw = temp.obj)
    temp.obj <- initRaw(temp.obj,GEO="" ,sc.method=" ",cond = "temp.clustered")
    
    n_cells <- length(get.cell.size(object = temp.obj))
    
    ttm <- clean(temp.obj)
    
    temp.obj <- ttm$object
    temp.obj <- hk_genes(temp.obj)
    temp.obj <- cotan_analysis(temp.obj, cores = 12)
    gc()
    temp.obj <- get.coex(temp.obj)
    gc()
    GDI_data <- get.GDI(temp.obj)
    
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




test_that("merge_cell.clusters.test", {
  temp <- readRDS(file.path(tm,"temp.RDS"))
  temp <- DEA_on_clusters(temp)
  temp <- temp[[1]]
  
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
    
    temp.obj <- new("scCOTAN",raw = temp.obj)
    temp.obj <- initRaw(temp.obj,GEO="" ,sc.method=" ",cond = "temp.clustered")
    
    n_cells <- length(get.cell.size(object = temp.obj))
    
    ttm <- clean(temp.obj)
    
    temp.obj <- ttm$object
    temp.obj <- hk_genes(temp.obj)
    temp.obj <- cotan_analysis(temp.obj, cores = 12)
    gc()
    temp.obj <- get.coex(temp.obj)
    gc()
    GDI_data <- get.GDI(temp.obj)
    
    expect_false( dim(GDI_data[GDI_data$GDI >= 1.5,])[1]/dim(GDI_data)[1] > 0.01 )
    
  }
  
  
    
})





