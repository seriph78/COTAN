# NB NON ho copiato il test pca nell'altro file!


#tm = tempdir()
#stopifnot(file.exists(tm))
#root = "tests/testthat/"
root = ""

test_that("python_PCA_test", {
    
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




