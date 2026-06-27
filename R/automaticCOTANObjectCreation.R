
# internal implementation

.proceedToCoexImpl <- function(objCOTAN,
                                calcCoex,
                                executionOptions,
                                cleaningOptions,
                                saveObj = FALSE,
                                outDir = ".") {
  startTime <- Sys.time()

  logThis("COTAN dataset analysis: START", logLevel = 1L)

  objCOTAN <- clean(objCOTAN, cleaningOptions = cleaningOptions)

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

  startEstimTime <- Sys.time()
  cleanTime <- difftime(startEstimTime, startTime, units = "secs")
  logThis(paste("Dataset cleaning elapsed time:", cleanTime),
          logLevel = 3L)

  objCOTAN <- estimateLambdaLinear(objCOTAN)
  objCOTAN <- estimateDispersionViaSolver(
    objCOTAN,
    executionOptions = executionOptions
  )

  gc()

  startCoexTime <- Sys.time()
  analysisTime <- difftime(startCoexTime, startEstimTime, units = "secs")
  logThis(paste("Model parameter estimation elapsed time:", analysisTime),
          logLevel = 3L)

  if (isTRUE(calcCoex)) {
    logThis("COTAN genes' COEX estimation: START", logLevel = 2L)

    objCOTAN <-
      calculateCoex(
        objCOTAN,
        actOnCells = FALSE,
        executionOptions = executionOptions
      )
  } else {
    logThis("COTAN genes' COEX estimation not requested", logLevel = 2L)
  }

  gc()

  endTime <- Sys.time()
  genesCoexTime <- difftime(endTime, startCoexTime, units = "secs")
  logThis(paste("Only genes' COEX elapsed time:", genesCoexTime),
          logLevel = 3L)

  totalTime <- difftime(endTime, startTime, units = "secs")
  logThis(paste("Dataset analysis elapsed time:", totalTime), logLevel = 3L)

  logThis("COTAN dataset analysis: DONE", logLevel = 1L)

  if (saveObj) {
    utils::write.csv(data.frame("type"  = c("total_time",
                                            "clean_time",
                                            "analysis_time",
                                            "genes_coex_time"),
                                "times" = c(as.numeric(totalTime),
                                            as.numeric(cleanTime),
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


.assertDefaultCleaningArgs <- function(cellsCutoff,
                                       genesCutoff,
                                       cellsThreshold,
                                       genesThreshold) {
  assert_that(
    identical(cellsCutoff, 0.003),
    identical(genesCutoff, 0.002),
    identical(cellsThreshold, 0.99),
    identical(genesThreshold, 0.99),
    msg = paste("Do not mix `cleaningOptions` with",
                "legacy cleaning arguments",
                "(`cellsCutoff`, `genesCutoff`,",
                "`cellsThreshold`, `genesThreshold`)."))
}

.assertDefaultExecutionArgs <- function(cores,
                                        optimizeForSpeed,
                                        deviceStr) {
  assert_that(
    identical(optimizeForSpeed, TRUE),
    identical(deviceStr, "cuda"),
    identical(as.integer(cores), 1L),
    msg = paste("Do not mix `executionOptions` with",
                "legacy execution arguments",
                "(`cores`, `optimizeForSpeed`, `deviceStr`)."))
}

# ------ proceedToCoex -------

#'
#' @aliases proceedToCoex
#'
#' @details `proceedToCoex()` takes a newly created `COTAN` object (or the
#'   result of a call to `dropGenesCells()`) and runs [calculateCoex()]
#'
#' @param objCOTAN a newly created `COTAN` object
#' @param calcCoex a Boolean to determine whether to calculate the genes' `COEX`
#'   or stop just before at the [estimateDispersionViaSolver()] step
#' @param optimizeForSpeed `r lifecycle::badge("deprecated")` Legacy execution
#'   scalar. When `TRUE`, `COTAN` tries to use accelerated `torch`-based matrix
#'   calculations when available; otherwise it falls back to the slower legacy
#'   code. Use
#'   `executionOptions = ExecutionOptions(optimizeForSpeed = ...)` instead.
#' @param deviceStr `r lifecycle::badge("deprecated")` Legacy execution scalar.
#'   Requested `torch` device string, for example `"cpu"`, `"cuda"`, or
#'   `"cuda:0"`. Use `executionOptions = ExecutionOptions(deviceStr = ...)`
#'   instead.
#' @param cores `r lifecycle::badge("deprecated")` Legacy execution scalar.
#'   Requested number of CPU cores. The effective value is bounded by the
#'   available cores. Use `executionOptions = ExecutionOptions(cores = ...)`
#'   instead.
#' @param cellsCutoff `r lifecycle::badge("deprecated")` Legacy cleaning scalar.
#'   `clean()` deletes from the `raw` data any gene expressed in fewer cells than
#'   this fraction times the total number of cells. Default cutoff is
#'   \eqn{0.003 \; (0.3\%)}. Use
#'   `cleaningOptions = CleaningOptions(cellsCutoff = ...)` instead.
#' @param genesCutoff `r lifecycle::badge("deprecated")` Legacy cleaning scalar.
#'   `clean()` deletes from the `raw` data any cell expressing fewer genes than
#'   this fraction times the total number of genes. Default cutoff is
#'   \eqn{0.002 \; (0.2\%)}. Use
#'   `cleaningOptions = CleaningOptions(genesCutoff = ...)` instead.
#' @param cellsThreshold `r lifecycle::badge("deprecated")` Legacy cleaning
#'   scalar. Any gene expressed in more cells than this fraction times the total
#'   number of cells is marked as **fully-expressed**. Default threshold is
#'   \eqn{0.99 \; (99.0\%)}. Use
#'   `cleaningOptions = CleaningOptions(cellsThreshold = ...)` instead.
#' @param genesThreshold `r lifecycle::badge("deprecated")` Legacy cleaning
#'   scalar. Any cell expressing more genes than this fraction times the total
#'   number of genes is marked as **fully-expressing**. Default threshold is
#'   \eqn{0.99 \; (99.0\%)}. Use
#'   `cleaningOptions = CleaningOptions(genesThreshold = ...)` instead.
#' @param cleaningOptions A `CleaningOptions` object bundling cleaning cutoffs
#'   and fully-expressed / fully-expressing thresholds. This is the preferred
#'   interface for new code.
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output.
#' @param executionOptions An `ExecutionOptions` object bundling the execution
#'   controls. This is the preferred interface for new code.
#'
#' @section Lifecycle:
#' Legacy scalar cleaning and execution arguments are soft-deprecated as of
#' COTAN 2.13.3. Use `cleaningOptions = CleaningOptions(...)` and
#' `executionOptions = ExecutionOptions(...)` in new code.
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
#' options(parallelly.fork.enable = TRUE)
#'
#' #
#' # In case one needs to run more steps to clean the datatset
#' # the following might apply
#' if (FALSE) {
#'   objCOTAN <- initializeMetaDataset(objCOTAN,
#'                                     GEO = "test",
#'                                     sequencingMethod = "artificial",
#'                                     sampleCondition = "test dataset")
#' #
#' # doing all the cleaning and analysis...
#' #
#'   exec <- ExecutionOptions(cores = 6L,
#'                            optimizeForSpeed = TRUE,
#'                            deviceStr = "cuda",
#'                            chunkSize = 1024L)
#'   cleanOpt <- CleaningOptions()
#'
#' # in case the genes' `COEX` is not needed it can be skipped
#' # (e.g. when calling [cellsUniformClustering()])
#' #
#'   objCOTAN <- proceedToCoex(
#'     objCOTAN,
#'     calcCoex = FALSE,
#'     executionOptions = exec,
#'     cleaningOptions = cleanOpt,
#'     saveObj = FALSE
#'   )
#' }
#'
#' @rdname COTAN_ObjectCreation
#'

setMethod(
  "proceedToCoex",
  signature(
    objCOTAN = "COTAN",
    executionOptions = "missing",
    cleaningOptions = "missing"),
  function(objCOTAN,
           calcCoex = TRUE,
           optimizeForSpeed = TRUE,
           deviceStr = "cuda",
           cores = 1L,
           cellsCutoff = 0.003,
           genesCutoff = 0.002,
           cellsThreshold = 0.99,
           genesThreshold = 0.99,
           cleaningOptions = NULL,
           saveObj = FALSE,
           outDir = ".",
           executionOptions = NULL) {
    callArgs <- names(as.list(match.call(expand.dots = FALSE))[-1L])

    .warnDeprecatedPackArgs(
      callArgs = callArgs,
      packClass = "CleaningOptions",
      replacementArg = "cleaningOptions",
      details = paste(
        "Use `cleaningOptions = CleaningOptions(...)` to configure",
        "cleaning thresholds."
      )
    )

    cleaningOptions <- resolveCleaningOptions(
      cellsCutoff = cellsCutoff,
      genesCutoff = genesCutoff,
      cellsThreshold = cellsThreshold,
      genesThreshold = genesThreshold,
      cleaningOptions = cleaningOptions
    )

    .warnDeprecatedPackArgs(
      callArgs = callArgs,
      packClass = "ExecutionOptions",
      replacementArg = "executionOptions",
      details = paste(
        "Use `executionOptions = ExecutionOptions(...)` to configure",
        "execution/backend parameters."
      )
    )

    executionOptions <- legacyExecutionOptions(
      cores = cores,
      optimizeForSpeed = optimizeForSpeed,
      deviceStr = deviceStr
    )

    objCOTAN <- .proceedToCoexImpl(
      objCOTAN = objCOTAN,
      calcCoex = calcCoex,
      executionOptions = executionOptions,
      cleaningOptions = cleaningOptions,
      saveObj = saveObj,
      outDir = outDir
    )

    return(objCOTAN)
  }
)



#' @details Alternative interface using an `ExecutionOptions` object.
#'
#' @rdname COTAN_ObjectCreation
#' @aliases proceedToCoex,COTAN,ExecutionOptions,missing-method
setMethod(
  "proceedToCoex",
  signature(
    objCOTAN = "COTAN",
    executionOptions = "ExecutionOptions",
    cleaningOptions = "missing"),
  function(objCOTAN,
           calcCoex = TRUE,
           optimizeForSpeed = TRUE,
           deviceStr = "cuda",
           cores = 1L,
           cellsCutoff = 0.003,
           genesCutoff = 0.002,
           cellsThreshold = 0.99,
           genesThreshold = 0.99,
           cleaningOptions = NULL,
           saveObj = FALSE,
           outDir = ".",
           executionOptions) {
    callArgs <- names(as.list(match.call(expand.dots = FALSE))[-1L])

    .warnDeprecatedPackArgs(
      callArgs = callArgs,
      packClass = "CleaningOptions",
      replacementArg = "cleaningOptions",
      details = paste(
        "Use `cleaningOptions = CleaningOptions(...)` to configure",
        "cleaning thresholds."
      )
    )

    cleaningOptions <- resolveCleaningOptions(
      cellsCutoff = cellsCutoff,
      genesCutoff = genesCutoff,
      cellsThreshold = cellsThreshold,
      genesThreshold = genesThreshold,
      cleaningOptions = cleaningOptions
    )

    .assertDefaultExecutionArgs(
      cores = cores,
      optimizeForSpeed = optimizeForSpeed,
      deviceStr = deviceStr
    )

    objCOTAN <- .proceedToCoexImpl(
      objCOTAN = objCOTAN,
      calcCoex = calcCoex,
      executionOptions = executionOptions,
      cleaningOptions = cleaningOptions,
      saveObj = saveObj,
      outDir = outDir
    )

    return(objCOTAN)
  }
)


#' @details Alternative interface using a `CleaningOptions` object.
#'
#' @rdname COTAN_ObjectCreation
#' @aliases proceedToCoex,COTAN,missing,CleaningOptions-method
setMethod(
  "proceedToCoex",
  signature(
    objCOTAN = "COTAN",
    executionOptions = "missing",
    cleaningOptions = "CleaningOptions"),
  function(objCOTAN,
           calcCoex = TRUE,
           optimizeForSpeed = TRUE,
           deviceStr = "cuda",
           cores = 1L,
           cellsCutoff = 0.003,
           genesCutoff = 0.002,
           cellsThreshold = 0.99,
           genesThreshold = 0.99,
           cleaningOptions,
           saveObj = FALSE,
           outDir = ".",
           executionOptions = NULL) {
    callArgs <- names(as.list(match.call(expand.dots = FALSE))[-1L])

    .assertDefaultCleaningArgs(
      cellsCutoff = cellsCutoff,
      genesCutoff = genesCutoff,
      cellsThreshold = cellsThreshold,
      genesThreshold = genesThreshold
    )

    .warnDeprecatedPackArgs(
      callArgs = callArgs,
      packClass = "ExecutionOptions",
      replacementArg = "executionOptions",
      details = paste(
        "Use `executionOptions = ExecutionOptions(...)` to configure",
        "execution/backend parameters."
      )
    )

    executionOptions <- legacyExecutionOptions(
      cores = cores,
      optimizeForSpeed = optimizeForSpeed,
      deviceStr = deviceStr
    )

    objCOTAN <- .proceedToCoexImpl(
      objCOTAN = objCOTAN,
      calcCoex = calcCoex,
      executionOptions = executionOptions,
      cleaningOptions = cleaningOptions,
      saveObj = saveObj,
      outDir = outDir
    )

    return(objCOTAN)
  }
)


#' @details Alternative interface using `ExecutionOptions` and `CleaningOptions`
#'   objects.
#'
#' @rdname COTAN_ObjectCreation
#' @aliases proceedToCoex,COTAN,ExecutionOptions,CleaningOptions-method
setMethod(
  "proceedToCoex",
  signature(
    objCOTAN = "COTAN",
    executionOptions = "ExecutionOptions",
    cleaningOptions = "CleaningOptions"),
  function(objCOTAN,
           calcCoex = TRUE,
           optimizeForSpeed = TRUE,
           deviceStr = "cuda",
           cores = 1L,
           cellsCutoff = 0.003,
           genesCutoff = 0.002,
           cellsThreshold = 0.99,
           genesThreshold = 0.99,
           cleaningOptions,
           saveObj = FALSE,
           outDir = ".",
           executionOptions) {
    .assertDefaultCleaningArgs(
      cellsCutoff = cellsCutoff,
      genesCutoff = genesCutoff,
      cellsThreshold = cellsThreshold,
      genesThreshold = genesThreshold
    )

    .assertDefaultExecutionArgs(
      cores = cores,
      optimizeForSpeed = optimizeForSpeed,
      deviceStr = deviceStr
    )

    objCOTAN <- .proceedToCoexImpl(
      objCOTAN = objCOTAN,
      calcCoex = calcCoex,
      executionOptions = executionOptions,
      cleaningOptions = cleaningOptions,
      saveObj = saveObj,
      outDir = outDir
    )

    return(objCOTAN)
  }
)


#' @details `automaticCOTANObjectCreation()` creates a `COTAN` object,
#'   initializes its dataset metadata, and then calls [proceedToCoex()] on it.
#'
#' @param raw a matrix or dataframe with the raw counts
#' @param GEO a code reporting the GEO identification or other specific dataset
#'   code
#' @param sequencingMethod a string reporting the method used for the sequencing
#' @param sampleCondition a string reporting the specific sample condition or
#'   time point
#' @param ... Additional arguments forwarded to [proceedToCoex()]
#'
#' @section Forwarded arguments:
#'   Additional arguments are forwarded to [proceedToCoex()]. Typical forwarded
#'   arguments include `calcCoex`, `cleaningOptions`, `executionOptions`,
#'   `saveObj` and `outDir`.
#'
#'   Legacy scalar cleaning and execution arguments accepted by [proceedToCoex()]
#'   are still forwarded for compatibility, but are soft-deprecated as of COTAN
#'   2.13.3.
#'
#' @returns `automaticCOTANObjectCreation()` returns a new `COTAN` object after
#'   initialization and analysis via [proceedToCoex()].
#'
#' @export
#'
#' @examples
#'
#' ## Otherwise it is possible to run all at once.
#' exec <- ExecutionOptions(cores = 6L,
#'                          optimizeForSpeed = TRUE,
#'                          deviceStr = "cuda",
#'                          chunkSize = 1024L)
#' cleanOpt <- CleaningOptions()
#'
#' objCOTAN <- automaticCOTANObjectCreation(
#'   raw = test.dataset,
#'   GEO = "code",
#'   sequencingMethod = "10X",
#'   sampleCondition = "mouse_dataset",
#'   calcCoex = TRUE,
#'   executionOptions = exec,
#'   cleaningOptions = cleanOpt,
#'   saveObj = FALSE
#' )
#'
#' @rdname COTAN_ObjectCreation

automaticCOTANObjectCreation <- function(raw,
                                         GEO,
                                         sequencingMethod,
                                         sampleCondition,
                                         ...) {
    objCOTAN <- COTAN(raw = raw)
    objCOTAN <- initializeMetaDataset(objCOTAN, GEO = GEO,
                                      sequencingMethod = sequencingMethod,
                                      sampleCondition = sampleCondition)

    logThis(paste0("Condition ", sampleCondition), logLevel = 2L)
    logThis(paste("n cells", getNumCells(objCOTAN)), logLevel = 2L)

    dots <- list(...)

    checkDotsAgainstFunction(
      dots = dots,
      fun = "proceedToCoex",
      allowRemaining = FALSE,
      allowUnnamed = FALSE
    )

    objCOTAN <- do.call(proceedToCoex, c(list(objCOTAN = objCOTAN), dots))

    return(objCOTAN)
  }

