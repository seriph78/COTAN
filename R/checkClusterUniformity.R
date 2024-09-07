
#' @details `checkClusterUniformity()` takes a `COTAN` object and a cells'
#'   *cluster* and checks whether the latter is **uniform** by `GDI`. The
#'   function runs `COTAN` to check whether the `GDI` is lower than the given
#'   `GDIThreshold` (1.43) for all but at the most `ratioAboveThreshold`
#'   (\eqn{1\%}) genes. If the `GDI` results to be too high for too many genes,
#'   the *cluster* is deemed
#'   **non-uniform**.
#'
#' @param objCOTAN a `COTAN` object
#' @param clusterName the tag of the *cluster*
#' @param cells the cells belonging to the *cluster*
#' @param checker the object that defines the method and the threshold to
#'   discriminate whether a *cluster* is *uniform transcript*. See
#'   [UniformTranscriptCheckers] for more details
#' @param cores number of cores to use. Default is 1.
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file(s)
#' @param outDir an existing directory for the analysis output. The effective
#'   output will be paced in a sub-folder.
#'
#' @returns `checkClusterUniformity` returns a checker object of the same type
#'   as the input one, that contains both threshold and results of the check:
#'   see [UniformTranscriptCheckers] for more details
#'
#' @importFrom utils head
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @importFrom withr local_options
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname UniformClusters
#'

checkClusterUniformity <- function(
    objCOTAN,
    clusterName,
    cells,
    checker,
    cores = 1L,
    optimizeForSpeed = TRUE,
    deviceStr = "cuda",
    saveObj = TRUE,
    outDir = ".") {
  # handle legacy usage
  cellsToDrop <- getCells(objCOTAN)[!getCells(objCOTAN) %in% cells]

  objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)

  objCOTAN <- proceedToCoex(objCOTAN, cores = cores,
                            optimizeForSpeed = optimizeForSpeed,
                            deviceStr = deviceStr, saveObj = FALSE)
  gc()

  checker@clusterSize <- getNumCells(objCOTAN)

  logThis(paste0("Checking uniformity for the cluster '", clusterName,
                 "' with ", checker@clusterSize, " cells"), logLevel = 2L)

  GDIData <- calculateGDI(objCOTAN)
  objCOTAN <- storeGDI(objCOTAN, GDIData)

  # Plots
  if (isTRUE(saveObj) && !dir.exists(outDir)) {
    saveObj <- FALSE
    warning(paste("Asked to save check results,",
                  "but given output folder does not exist"))
  }

  if (isTRUE(saveObj)) tryCatch({
      # this will be restored to previous value on function exit
      local_options(list(ggrepel.max.overlaps = Inf))

      pdf(file.path(outDir, paste0("cluster_", clusterName, "_plots.pdf")))

      c(..., nuPlot, zoomedNuPlot) %<-%
        cleanPlots(objCOTAN, includePCA = FALSE)

      genesToLabel <-
        head(rownames(GDIData[order(GDIData[["GDI"]],
                                    decreasing = TRUE), ]), n = 10L)
      gdiPlot <- GDIPlot(objCOTAN, GDIIn = GDIData,
                         GDIThreshold = getCheckerThreshold(checker),
                         genes = list("top 10 GDI genes" = genesToLabel))

      plot(nuPlot)
      plot(zoomedNuPlot)
      plot(gdiPlot)

      rm(nuPlot, zoomedNuPlot, gdiPlot)
      dev.off()
    }, error = function(err) {
      logThis(paste("While saving cluster plots", err), logLevel = 0L)
      dev.off()
    }
  )

  checker <- checkObjIsUniform(currentC = checker, previousC = NULL,
                               objCOTAN = objCOTAN)
  rm(objCOTAN)
  gc()

  logThis(paste0(
    "Cluster ", clusterName, ", with size ", checker@clusterSize, ", is ",
    ifelse(checker@isUniform, "", "not "), "uniform"), logLevel = 1)

  if (TRUE) {
    dumpDF <- checkersToDF(checker)
    logThis(paste0(colnames(dumpDF), " = ", dumpDF, collapse = ", "),
            logLevel = 3L)
    rm(dumpDF)
  }

  if (isTRUE(saveObj)) tryCatch({
      outFile <- file.path(outDir,
                           paste0(ifelse(checker@isUniform, "", "non-"),
                                  "uniform_cluster_", clusterName, ".csv"))
      write.csv(cells, file = outFile)
    },
    error = function(err) {
      logThis(paste("While saving current clusterization", err),
              logLevel = 0L)
    }
  )

  return(checker)
}
