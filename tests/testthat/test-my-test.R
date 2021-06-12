tm = tempdir()
stopifnot(file.exists(tm))
#root = "tests/testthat/"
root = ""
test_that("1.initialization", {
    #raw = read.csv(file.path(getwd(),"raw.csv"),header = T,row.names = 1)
    #raw = readRDS(file.path(getwd(),"raw.RDS"))
    #raw = readRDS("tests/testthat/raw.RDS")
    utils::data("ERCCraw", package = "COTAN")
    #print(data[1:5,1:5])
    rownames(data) = data$V1
    data = data[,2:ncol(data)]
    obj.temp = new("scCOTAN",raw = data)
    obj.temp = initRaw(object = obj.temp,GEO="V" ,sc.method="10X",cond = "example")


    #expect_equal(obj.temp,readRDS(file.path(getwd(),"obj.RDS")) )
    expect_s4_class(obj.temp,"scCOTAN")

})

test_that("2.cleaning", {

    #obj.temp = readRDS(file.path(getwd(),"obj.RDS"))
    utils::data("raw.dataset", package = "COTAN")

    obj.temp = new("scCOTAN",raw = raw.dataset)
    obj.temp = initRaw(object = obj.temp,GEO=" " ,sc.method="10X",cond = "example")
    #---------------------------------------------------

    ttm = clean(obj.temp)
    stopifnot(file.exists(tm))
    saveRDS(ttm$object, file = file.path(tm,"temp.RDS") )
    raw.norm = readRDS(file.path(getwd(),"raw.norm.RDS"))
    nu = readRDS(file.path(getwd(),"nu.RDS"))
    expect_equal(as.matrix(ttm$object@raw.norm),raw.norm)
    expect_equal(as.matrix(ttm$object@nu),nu)

})

test_that("mat_division", {
    utils::data("raw.dataset", package = "COTAN")
    #print(raw.dataset[1:5,1:5])
    #raw = readRDS(file.path(getwd(),"raw.RDS"))
    nu = readRDS(file.path(getwd(),"nu.RDS"))
    raw.norm = readRDS(file.path(getwd(),"raw.norm.RDS"))

    expect_equal( (t(t(raw.dataset) * (1/nu[,1]))),as.matrix(raw.norm) )
})


test_that("3.cotan_analysis_test", {

    obj=readRDS(file.path(tm,"temp.RDS"))
    obj = cotan_analysis(obj,cores = 2)

    a = readRDS(file.path(getwd(),"a.RDS"))
    nu = readRDS(file.path(getwd(),"nu.RDS"))
    lambda = readRDS(file.path(getwd(),"lambda.RDS"))
    saveRDS(obj, file = file.path(tm,"temp.RDS") )
    expect_equal(as.matrix(obj@a), a)
    expect_equal( as.matrix(obj@nu) , nu)
    expect_equal( as.matrix(obj@lambda), lambda)
})


test_that("4.cotan_coex_test", {
    obj=readRDS(file.path(tm,"temp.RDS"))

    obj = get.coex(obj)

    coex = readRDS(file.path(getwd(),"coex.RDS"))
    saveRDS(obj, file = file.path(tm,"temp.RDS") )
    expect_equal( as.matrix(obj@coex[1:500,1:200]), as.matrix(coex))
})

test_that("python_PCA_test", {
    #raw = readRDS(file.path(getwd(),"raw.RDS"))
    utils::data("raw.dataset", package = "COTAN")
    #file.py <- system.file("inst","python", "python_PCA.py", package="COTAN")

    proc <- basiliskStart(my_env_cotan)
    on.exit(basiliskStop(proc))

file.py <- system.file("python/python_PCA.py", package="COTAN",mustWork = TRUE)

    pca.raw <- basiliskRun(proc, function(arg1) {
        reticulate::source_python(file.py)
        output <- python_PCA(arg1)

        # The return value MUST be a pure R object, i.e., no reticulate
        # Python objects, no pointers to shared memory.
        output
    }, arg1=raw.dataset)

    rownames(pca.raw)=rownames(raw.dataset)
    colnames(pca.raw)=paste0("PC_", c(1:10))
    pca.raw =as.data.frame(round(pca.raw, digits = 14))
    pca.tb = readRDS(file.path(getwd(),"pca.mat.RDS"))

    expect_equal( pca.raw$PC_1, pca.tb$PC_1)
    expect_equal( pca.raw$PC_2, pca.tb$PC_2)

})


test_that("5.get_pval_test", {
    object=readRDS(file.path(tm,"temp.RDS"))
    #object = readRDS(file.path(getwd(),"Obj_out_cotan_coex_not_approx.RDS"))
    object@coex = object@coex/sqrt(object@n_cells) #Because this object was created
    #with the old method and it was divided by sqrt(cell_number)
    pval = get.pval(object, gene.set.col = rownames(object@raw)[1:100], rownames(object@raw)[1:100])
    #pval.exp = readRDS(file.path(getwd(),"pval.RDS"))
    pval.exp =readRDS(file.path(getwd(),"pval.RDS"))
    #pval.exp = as.data.frame(pval.exp)
    #rownames(pval.exp) =pval.exp$V1
    #pval.exp = pval.exp[,2:ncol(pval.exp)]
    expect_equal( pval, pval.exp)

})


test_that("get_GDI_test", {
    #object = readRDS(file.path(getwd(),"ERCC.cotan.RDS"))
    utils::data("ERCC.cotan", package = "COTAN")
    GDI = get.GDI(ERCC.cotan)
    expect_equal(GDI, readRDS(file.path(getwd(),"GDI.RDS")))

})
