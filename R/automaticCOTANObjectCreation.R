
#'
#' Automatic `COTAN` shortcuts
#'
#' @description These functions take (or create) a `COTAN` object and run all
#'   the necessary steps until the genes' `COEX` matrix is calculated.
#'
#'   takes a newly created `COTAN` object (or the result of a call to
#'   [dropGenesCells()]) and applies all steps until the genes' `COEX` matrix is
#'   stored in the object
#'
#' @name COTANObjectCreation
NULL

#'
#' @aliases proceedToCoex
#'
#' @details `proceedToCoex()` takes a newly created `COTAN` object (or the
#'   result of a call to `dropGenesCells()`) and runs [calculateCoex()]
#'
#' @param objCOTAN a newly created `COTAN` object
#' @param calcCoex a Boolean to determine whether to calculate the genes' `COEX`
#'   or stop just before at the [estimateDispersionBisection()] step
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
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
#' ## # in case the genes' `COEX` is not needed it can be skipped
#' ## # (e.g. for [cellsUniformClustering()])
#' ## objCOTAN <- proceedToCoex(objCOTAN, calcCoex = FALSE,
#' ##                           cores = 6L, optimizeForSpeed, deviceStr, saveObj = FALSE)
#'
#' @rdname COTANObjectCreation
#'
setMethod(
  "proceedToCoex",
  "COTAN",
  function(objCOTAN, calcCoex = TRUE, optimizeForSpeed = TRUE,
           deviceStr = "cuda", cores = 1L, saveObj = TRUE, outDir = ".") {
    startTimeAll <- Sys.time()

    logThis("Cotan analysis functions started", logLevel = 1L)

    objCOTAN <- clean(objCOTAN)

    if (isTRUE(saveObj)) {
      if (!file.exists(outDir)) {
        dir.create(file.path(outDir))
      }

      cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

      outDirCond <- file.path(outDir, cond)
      if (!file.exists(outDirCond)) {
        dir.create(outDirCond)
      }

      outDirCleaning <- file.path(outDirCond, "cleaning")
      if (!file.exists(outDirCleaning)) {
        dir.create(outDirCleaning)
      }

      plots <- cleanPlots(objCOTAN)

      {
        numIter <- 1L
        pdf(file.path(outDirCleaning,
                      paste0(cond, "_", numIter,
                             "_plots_without_cleaning.pdf")))
        plot(plots[["pcaCells"]])
        plot(plots[["genes"]])
        dev.off()
      }

      {
        pdf(file.path(outDirCleaning,
                      paste0(cond, "_plots_PCA_efficiency_colored.pdf")))
        plot(plots[["UDE"]])
        dev.off()
      }

      {
        pdf(file.path(outDirCleaning,
                      paste0(cond, "_plots_efficiency.pdf")))
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

    if (isTRUE(calcCoex)) {
      genesCoexTime <- Sys.time()
      analysisTime <- difftime(genesCoexTime, analysisTime, units = "secs")

      logThis(paste("Only analysis elapsed time:", analysisTime), logLevel = 3L)

      logThis("Cotan genes' coex estimation started", logLevel = 1L)
      objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE,
                                optimizeForSpeed = optimizeForSpeed,
                                deviceStr = deviceStr)

      gc()

      endTime <- Sys.time()

      genesCoexTime <- difftime(endTime, genesCoexTime, units = "secs")

      logThis(paste("Only genes' coex elapsed time:", genesCoexTime), logLevel = 3L)
    } else {
      logThis("Cotan genes' coex estimation not requested", logLevel = 2L)

      genesCoexTime <- 0.0

      gc()

      endTime <- Sys.time()
    }

    allTime <- difftime(endTime, startTimeAll, units = "secs")
    logThis(paste("Total elapsed time:", allTime), logLevel = 3L)

    if (saveObj) {
      utils::write.csv(data.frame("type"  = c("tot_time",
                                              "analysis_time",
                                              "genes_coex_time"),
                                  "times" = c(as.numeric(allTime),
                                              as.numeric(analysisTime),
                                              as.numeric(genesCoexTime)),
                                  "n.cells" = getNumCells(objCOTAN),
                                  "n.genes" = getNumGenes(objCOTAN)),
                       file = file.path(outDir, paste0(cond, "_times.csv")))

      logThis(paste("Saving elaborated data locally at:",
                    file.path(outDir, paste0(cond, ".cotan.RDS"))),
              logLevel = 1L)
      saveRDS(objCOTAN, file = file.path(outDir, paste0(cond, ".cotan.RDS")))
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
#' @param calcCoex a Boolean to determine whether to calculate the genes' `COEX`
#'   or stop just before at the [estimateDispersionBisection()] step
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
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
#'   sampleCondition = "mouse_dataset",
#'   calcCoex = TRUE,
#'   saveObj = FALSE,
#'   outDir = tempdir(),
#'   cores = 6L)
#'
#' @rdname COTANObjectCreation

automaticCOTANObjectCreation <-
  function(raw, GEO, sequencingMethod, sampleCondition,
           calcCoex = TRUE, optimizeForSpeed = TRUE, deviceStr = "cuda",
           cores = 1L, saveObj = TRUE, outDir = ".") {

  objCOTAN <- COTAN(raw = raw)
  objCOTAN <- initializeMetaDataset(objCOTAN, GEO = GEO,
                                    sequencingMethod = sequencingMethod,
                                    sampleCondition = sampleCondition)

  logThis(paste0("Condition ", sampleCondition), logLevel = 2L)
  logThis(paste("n cells", getNumCells(objCOTAN)), logLevel = 2L)

  return(proceedToCoex(objCOTAN, calcCoex = calcCoex,
                       optimizeForSpeed = optimizeForSpeed,
                       deviceStr = deviceStr, cores = cores,
                       saveObj = saveObj, outDir = outDir))
}
