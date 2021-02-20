#root = "../../"
root = ""
test_that("initialization", {
    obj.temp = scCOTAN(raw = read.csv(paste0(root,"raw.csv"), header = T,row.names = 1))
    obj.temp = initRaw(obj.temp,GEO=" " ,sc.method="10X",cond = "example")

    file.py <- system.file("inst","python", "python_PCA.py", package="COTAN")
    print(file.py)

    expect_equal(obj.temp,readRDS(paste0(root,"obj.RDS")) )

})

test_that("cleaning_normal", {

    obj.temp = readRDS(paste0(root,"obj.RDS"))
    cells =as.data.frame(as.matrix(obj.temp@raw))
    #---------------------------------------------------
    cells[cells > 0] <- 1
    cells[cells <= 0] <- 0
    #cells = cells[rowSums(cells)> round((length(colnames(obj.temp@raw))/1000*3), digits = 0),]
    ttm = clean(obj.temp, cells,mean_type = "normal")
    expect_equal(ttm$object,readRDS(paste0(root,"ttmobj2.RDS")))
})

test_that("mat_division", {

    mat1 = read.csv(paste0(root,"mat1.csv"),header = T,row.names = 1)
    nu1 = read.csv(paste0(root,"nu1.csv"),header = T,row.names = 1)
    nu1 =nu1$x
    res =  read.csv(paste0(root,"res.1csv"),header = T,row.names = 1)

    expect_equal( (t(t(mat1) * (1/nu1))),as.matrix(res))
})


test_that("cotan_analysis_test", {
    obj=readRDS(paste0(root,"ttmobj2.RDS"))
    obj = cotan_analysis(obj)
    expect_equal( obj,readRDS(paste0(root,"Obj_out_cotan_an.RDS")))
})

test_that("cotan_coex_test", {
    obj=readRDS(paste0(root,"Obj_out_cotan_an.RDS"))
    obj = get.coex(obj)
    expect_equal( obj,readRDS(paste0(root,"Obj_out_cotan_coex.RDS")))
})

test_that("python_PCA_test", {
    raw = read.csv(paste0(root,"raw.csv"), header = T,row.names = 1)
    #pca.raw = python_PCA(raw)
    file.py <- system.file("inst","python", "python_PCA.py", package="COTAN")

    proc <- basiliskStart(my_env_cotan)
    on.exit(basiliskStop(proc))

    pca.raw <- basiliskRun(proc, function(arg1) {
        reticulate::source_python(file.py)
        output <- python_PCA(arg1)

        # The return value MUST be a pure R object, i.e., no reticulate
        # Python objects, no pointers to shared memory.
        output
    }, arg1=raw)

    rownames(pca.raw)=rownames(raw)
    colnames(pca.raw)=paste0("PC_", c(1:10))
    pca.raw =as.data.frame(round(pca.raw, digits = 14))

    pca.tb = (read.csv(paste0(root,"pca.mat.csv"),header = T,row.names = 1))
    expect_lt(sum((pca.raw[,1:4] - pca.tb[,1:4])**2),10**(-5))
    ### È molto strano! la PCA fatto così ha un qualche grado di randomness!
})


