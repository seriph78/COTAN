
#' Automatic `COTAN` shortcuts
#'
#' @description These functions take (or create) a `COTAN` object and run all
#'   the necessary steps until the genes' `COEX` matrix is calculated.
#'
#'   takes a newly created `COTAN` object (or the result of a call to
#'   [dropGenesCells()]) and applies all steps until the genes' `COEX` matrix is
#'   stored in the object
#'
#' @details `proceedToCoex()` takes a newly created `COTAN` object (or the
#'   result of a call to `dropGenesCells()`) and runs [calculateCOEX()]
#'
#' @param objCOTAN a newly created `COTAN` object
#' @param cores number of cores to be used
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output.
#'
#' @returns `proceedToCoex()` returns the updated `COTAN` object with genes'
#'   `COEX` calculated. If asked to, it will also store the object, along all
#'   relevant clean-plots, in the output directory.
#'
#' @export
#'
#' @importFrom utils write.csv
#'
#' @importFrom stats time
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom grDevices colors
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 annotate
#'
#' @examples
#' data("test.dataset")
#'
#' ## In case one needs to run more steps to clean the datatset the following
#' ## might apply
#' ##
#' ## objCOTAN <- COTAN(raw = test.dataset)
#' ## objCOTAN <- initializeMetaDataset(objCOTAN,
#' ##                                   GEO = "test",
#' ##                                   sequencingMethod = "artificial",
#' ##                                   sampleCondition = "test dataset")
#' ## objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)
#'
#' @rdname automaticCOTANObjectCreation
#'
setMethod(
  "proceedToCoex",
  "COTAN",
  function(objCOTAN, cores = 1L, saveObj = TRUE, outDir = ".") {
    startTimeAll <- Sys.time()

    logThis("Cotan analysis functions started", logLevel = 1L)

    objCOTAN <- clean(objCOTAN)

    if (isTRUE(saveObj)) {
      if (!file.exists(outDir)) {
        dir.create(file.path(outDir))
      }

      outDirCond <- file.path(outDir, cond)
      if (!file.exists(outDirCond)) {
        dir.create(outDirCond)
      }

      outDirCleaning <- file.path(outDirCond, "cleaning")
      if (!file.exists(outDirCleaning)) {
        dir.create(outDirCleaning)
      }

      plots <- cleanPlots(objCOTAN)

      sampleCondition <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])
      {
        numIter <- 1L
        pdf(file.path(outDirCleaning,
                      paste0(sampleCondition, "_",
                             numIter, "_plots_without_cleaning.pdf")))
        plot(plots[["pcaCells"]])
        plot(plots[["genes"]])
        dev.off()
      }

      {
        pdf(file.path(outDirCleaning,
                      paste0(sampleCondition,
                             "_plots_PCA_efficiency_colored.pdf")))
        plot(plots[["UDE"]])
        dev.off()
      }

      {
        pdf(file.path(outDirCleaning,
                      paste0(sampleCondition,
                             "_plots_efficiency.pdf")))
        plot(plots[["nu"]] +
               annotate(geom = "text", x = 50L, y = 0.25,
                        label = "nothing to remove ", color = "darkred"))
        dev.off()
      }

      rm(plots)
      gc()
    }

    analysisTime <- Sys.time()

    objCOTAN <- estimateDispersionBisection(objCOTAN, cores = cores)

    gc()

    genesCoexTime <- Sys.time()
    analysisTime <- difftime(genesCoexTime, analysisTime, units = "mins")

    logThis(paste0("Only analysis time ", analysisTime), logLevel = 3L)

    logThis("Cotan genes' coex estimation started", logLevel = 1L)
    objCOTAN <- calculateCoex(objCOTAN)

    gc()

    endTime <- Sys.time()

    allTime <- difftime(endTime, startTimeAll, units = "mins")
    logThis(paste0("Total time ", allTime), logLevel = 3L)

    genesCoexTime <- difftime(endTime, genesCoexTime, units = "mins")

    logThis(paste0("Only genes' coex time ", genesCoexTime), logLevel = 3L)

    if (saveObj) {
      utils::write.csv(data.frame("type"  = c("tot_time",
                                              "analysis_time",
                                              "genes_coex_time"),
                                  "times" = c(as.numeric(allTime),
                                              as.numeric(analysisTime),
                                              as.numeric(genesCoexTime)),
                                  "n.cells" = getNumCells(objCOTAN),
                                  "n.genes" = getNumGenes(objCOTAN)),
                       file = file.path(outDir, paste0(sampleCondition,
                                                       "_times.csv")))

      logThis(paste0("Saving elaborated data locally at: ",
                     file.path(outDir, paste0(sampleCondition, ".cotan.RDS"))),
              logLevel = 1L)
      saveRDS(objCOTAN, file = file.path(outDir, paste0(sampleCondition,
                                                        ".cotan.RDS")))
    }

    return(objCOTAN)
  }
)

#' @details `automaticCOTANObjectCreation()` takes a raw dataset, creates and
#'   initializes a `COTAN` objects and runs proceedToCoex()
#'
#' @param raw a matrix or dataframe with the raw counts
#' @param GEO a code reporting the GEO identification or other specific dataset
#'   code
#' @param sequencingMethod a string reporting the method used for the sequencing
#' @param sampleCondition a string reporting the specific sample condition or
#'   time point.
#' @param cores number of cores to be used
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output.
#'
#' @returns `automaticCOTANObjectCreation()` returns the new `COTAN` object with
#'   genes' `COEX` calculated. When asked, it will also store the object, along
#'   all relevant clean-plots, in the output directory.
#'
#' @export
#'
#' @examples
#'
#' ## Otherwise it is possible to run all at once.
#' objCOTAN <- automaticCOTANObjectCreation(
#'   raw = test.dataset,
#'   GEO = "code",
#'   sequencingMethod = "10X",
#'   sampleCondition = "mouse dataset",
#'   saveObj = FALSE,
#'   outDir = tempdir(),
#'   cores = 12)
#'
#' @rdname automaticCOTANObjectCreation

automaticCOTANObjectCreation <-
  function(raw, GEO, sequencingMethod, sampleCondition,
           cores = 1L, saveObj = TRUE, outDir = ".") {

  objCOTAN <- COTAN(raw = raw)
  objCOTAN <- initializeMetaDataset(objCOTAN, GEO = GEO,
                               sequencingMethod = sequencingMethod,
                               sampleCondition = sampleCondition)

  #if (isFALSE(mt)) {
  #  genes_to_rem <- getGenes(objCOTAN)[grep(mt_prefix, getGenes(objCOTAN)),])
  #  cells_to_rem <- getCells(objCOTAN)[which(getCellsSize(objCOTAN) == 0L)])
  #  dropGenesCells(objCOTAN, genes_to_rem, cells_to_rem)
  #}

  logThis(paste0("Condition ", sampleCondition), logLevel = 2L)
  logThis(paste("n cells", getNumCells(objCOTAN)), logLevel = 2L)

  return(proceedToCoex(objCOTAN, cores = cores,
                       saveObj = saveObj, outDir = outDir))
}
