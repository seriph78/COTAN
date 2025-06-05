
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
#' @param cores number of cores to use. Default is 1.
#' @param cellsCutoff `clean()` will delete from the `raw` data any gene that is
#'   expressed in less cells than threshold times the total number of cells.
#'   Default cutoff is \eqn{0.003 \; (0.3\%)}
#' @param genesCutoff `clean()` will delete from the `raw` data any cell that is
#'   expressing less genes than threshold times the total number of genes.
#'   Default cutoff is \eqn{0.002 \; (0.2\%)}
#' @param cellsThreshold any gene that is expressed in more cells than threshold
#'   times the total number of cells will be marked as **fully-expressed**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
#' @param genesThreshold any cell that is expressing more genes than threshold
#'   times the total number of genes will be marked as **fully-expressing**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
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
#' @importFrom grDevices dev.cur
#'
#' @importFrom ggplot2 annotate
#'
#' @examples
#' #
#' # In case one needs to run more steps to clean the datatset
#' # the following might apply
#' if (FALSE) {
#'   objCOTAN <- initializeMetaDataset(objCOTAN,
#'                                     GEO = "test",
#'                                     sequencingMethod = "artificial",
#'                                     sampleCondition = "test dataset")
#' #
#' # doing all the cleaning...
#' #
#' # in case the genes' `COEX` is not needed it can be skipped
#' # (e.g. when calling [cellsUniformClustering()])
#'   objCOTAN <- proceedToCoex(objCOTAN, calcCoex = FALSE,
#'                             cores = 6L, optimizeForSpeed = TRUE,
#'                             deviceStr = "cuda", saveObj = FALSE)
#' }
#'
#' @rdname COTAN_ObjectCreation
#'
setMethod(
  "proceedToCoex",
  "COTAN",
  function(objCOTAN, calcCoex = TRUE, optimizeForSpeed = TRUE,
           deviceStr = "cuda", cores = 1L,
           cellsCutoff = 0.003, genesCutoff = 0.002,
           cellsThreshold = 0.99, genesThreshold = 0.99,
           saveObj = TRUE, outDir = ".") {
    startTimeAll <- Sys.time()

    logThis("Cotan analysis functions started", logLevel = 1L)

    objCOTAN <- clean(objCOTAN, cellsCutoff, genesCutoff,
                      cellsThreshold, genesThreshold)

    if (isTRUE(saveObj)) tryCatch({
      if (!dir.exists(outDir)) {
        dir.create(file.path(outDir))
      }

      cond <- getMetadataElement(objCOTAN, datasetTags()[["cond"]])

      outDirCond <- file.path(outDir, cond)
      if (!dir.exists(outDirCond)) {
        dir.create(outDirCond)
      }

      outDirCleaning <- file.path(outDirCond, "cleaning")
      if (!dir.exists(outDirCleaning)) {
        dir.create(outDirCleaning)
      }

      plots <- cleanPlots(objCOTAN)

      if (TRUE) {
        numIter <- 1L
        pdf(file.path(outDirCleaning,
                      paste0(cond, "_", numIter,
                             "_plots_without_cleaning.pdf")))
        plot(plots[["pcaCells"]])
        plot(plots[["genes"]])
        dev.off()
      }

      if (TRUE) {
        pdf(file.path(outDirCleaning,
                      paste0(cond, "_plots_PCA_efficiency_colored.pdf")))
        plot(plots[["UDE"]])
        dev.off()
      }

      if (TRUE) {
        pdf(file.path(outDirCleaning,
                      paste0(cond, "_plots_efficiency.pdf")))
        plot(plots[["nu"]] +
               annotate(geom = "text", x = 50L, y = 0.25,
                        label = "nothing to remove ", color = "darkred"))
        dev.off()
      }

      rm(plots)
    }, error = function(err) {
      logThis(paste("While saving the clean plots", err), logLevel = 1L)
    }, finally = {
      # Check for active device
      if (dev.cur() > 1L) {
        dev.off()
      }
    })

    gc()

    analysisTime <- Sys.time()

    objCOTAN <- estimateLambdaLinear(objCOTAN)
    objCOTAN <- estimateDispersionBisection(objCOTAN, cores = cores)

    gc()

    if (isTRUE(calcCoex)) {
      genesCoexTime <- Sys.time()
      analysisTime <- difftime(genesCoexTime, analysisTime, units = "secs")

      logThis(paste("Only analysis elapsed time:", analysisTime), logLevel = 3L)

      logThis("Cotan genes' COEX estimation started", logLevel = 1L)
      objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE,
                                optimizeForSpeed = optimizeForSpeed,
                                deviceStr = deviceStr)

      gc()

      endTime <- Sys.time()

      genesCoexTime <- difftime(endTime, genesCoexTime, units = "secs")

      logThis(paste("Only genes' COEX elapsed time:", genesCoexTime),
              logLevel = 3L)
    } else {
      logThis("Cotan genes' COEX estimation not requested", logLevel = 2L)

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
#'   initializes a `COTAN` object and runs [proceedToCoex()]
#'
#' @param raw a matrix or dataframe with the raw counts
#' @param GEO a code reporting the GEO identification or other specific dataset
#'   code
#' @param sequencingMethod a string reporting the method used for the sequencing
#' @param sampleCondition a string reporting the specific sample condition or
#'   time point.
#' @param calcCoex a Boolean to determine whether to calculate the genes' `COEX`
#'   or stop just after the [estimateDispersionBisection()] step
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
#' @param cores number of cores to use. Default is 1.
#' @param cellsCutoff `clean()` will delete from the `raw` data any gene that is
#'   expressed in less cells than threshold times the total number of cells.
#'   Default cutoff is \eqn{0.003 \; (0.3\%)}
#' @param genesCutoff `clean()` will delete from the `raw` data any cell that is
#'   expressing less genes than threshold times the total number of genes.
#'   Default cutoff is \eqn{0.002 \; (0.2\%)}
#' @param cellsThreshold any gene that is expressed in more cells than threshold
#'   times the total number of cells will be marked as **fully-expressed**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
#' @param genesThreshold any cell that is expressing more genes than threshold
#'   times the total number of genes will be marked as **fully-expressing**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
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
#' @rdname COTAN_ObjectCreation

automaticCOTANObjectCreation <-
  function(raw, GEO, sequencingMethod, sampleCondition,
           calcCoex = TRUE, optimizeForSpeed = TRUE,
           deviceStr = "cuda", cores = 1L,
           cellsCutoff = 0.003, genesCutoff = 0.002,
           cellsThreshold = 0.99, genesThreshold = 0.99,
           saveObj = TRUE, outDir = ".") {

    objCOTAN <- COTAN(raw = raw)
    objCOTAN <- initializeMetaDataset(objCOTAN, GEO = GEO,
                                      sequencingMethod = sequencingMethod,
                                      sampleCondition = sampleCondition)

    logThis(paste0("Condition ", sampleCondition), logLevel = 2L)
    logThis(paste("n cells", getNumCells(objCOTAN)), logLevel = 2L)

    return(proceedToCoex(objCOTAN, calcCoex = calcCoex,
                         optimizeForSpeed = optimizeForSpeed,
                         deviceStr = deviceStr, cores = cores,
                         cellsCutoff, genesCutoff,
                         cellsThreshold, genesThreshold,
                         saveObj = saveObj, outDir = outDir))
  }
