
setMethod(
  "checkObjIsUniform",
  signature(objCOTAN = "COTAN",
            currentChecker = "SimpleGDIUniformityCheck",
            usedChecker = "SimpleGDIUniformityCheck"),
  function(objCOTAN, currentChecker, usedChecker) {
    invisible(validObject(currentChecker))
    if(validObject(usedChecker)) {

    } else {
      GDIData <- calculateGDI(objCOTAN)

      GDIThreshold
      ratioAboveThreshold
    }


  }
)


#'
#' @details `isClusterUniform()` takes in the current thresholds and used them
#'   to check whether the calculated cluster parameters are sufficient to
#'   determine whether the cluster is **uniform** and in the positive scenario
#'   the corresponding answer
#'
#' @param GDIThreshold the threshold level that discriminates uniform
#'   *clusters*. It defaults to \eqn{1.43}
#' @param ratioAboveThreshold the fraction of genes allowed to be above the
#'   `GDIThreshold`. It defaults to \eqn{1\%}
#' @param ratioQuantile the `GDI` quantile corresponding to the `usedRatioAbove`
#' @param fractionAbove the fraction of genes above the `usedGDIThreshold`
#' @param usedGDIThreshold the threshold level actually used to calculate fourth
#'   argument
#' @param usedRatioAbove the fraction of genes actually used to calculate the
#'   third argument
#'
#' @returns a single `Boolean` value when it is possible to decide the answer
#'   with the given information and `NA` otherwise
#'
#' @importFrom assertthat assert_that
#'
#' @rdname UniformClusters
#'

isClusterUniform <- function(GDIThreshold, ratioAboveThreshold,
                             ratioQuantile, fractionAbove,
                             usedGDIThreshold, usedRatioAbove) {
  assert_that(!is.na(GDIThreshold), !is.na(ratioAboveThreshold),
              !is.na(usedGDIThreshold), !is.na(usedRatioAbove),
              GDIThreshold >= 0.0, ratioAboveThreshold >= 0.0,
              ratioAboveThreshold <= 1.0, msg = "wrong thresholds passed in")

  if (!is.na(fractionAbove) && GDIThreshold == usedGDIThreshold) {
    return(fractionAbove <= ratioAboveThreshold)
  } else if (!is.na(ratioQuantile) && ratioAboveThreshold == usedRatioAbove) {
    return(ratioQuantile < GDIThreshold)
  } else {
    return(NA)
  }
}



#'
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
#' @param GDIThreshold the threshold level that discriminates uniform
#'   *clusters*. It defaults to \eqn{1.43}
#' @param ratioAboveThreshold the fraction of genes allowed to be above the
#'   `GDIThreshold`. It defaults to \eqn{1\%}
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
#' @returns `checkClusterUniformity` returns a list with:
#'   * `"isUniform"`: a flag indicating whether the *cluster* is **uniform**
#'   * `"fractionAbove"`: the percentage of genes with `GDI` above the threshold
#'   * `"ratioQuantile"`: the quantile associated to the high quantile
#'   associated to given ratio
#'   * `"size"`: the number of cells in the cluster
#'   * `"GDIThreshold"` the used `GDI` threshold
#'   * `"ratioAboveThreshold"` the used fraction of genes above threshold
#'     allowed in **uniform** *clusters*
#'
#' @importFrom utils head
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @importFrom withr local_options
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
    checker = NULL,
    GDIThreshold = NULL,
    cores = 1L,
    optimizeForSpeed = TRUE,
    deviceStr = "cuda",
    saveObj = TRUE,
    outDir = ".") {
  cellsToDrop <- getCells(objCOTAN)[!getCells(objCOTAN) %in% cells]

  objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)

  objCOTAN <- proceedToCoex(objCOTAN, cores = cores,
                            optimizeForSpeed = optimizeForSpeed,
                            deviceStr = deviceStr, saveObj = FALSE)
  gc()

  clusterSize <- getNumCells(objCOTAN)

  logThis(paste0("Checking uniformity for the cluster '", clusterName,
                 "' with ", clusterSize, " cells"), logLevel = 2L)

  GDIData <- calculateGDI(objCOTAN)

  # Plots
  if (isTRUE(saveObj)) tryCatch({
      # this will be restored to previous value on function exit
      local_options(list(ggrepel.max.overlaps = Inf))

      pdf(file.path(outDir, paste0("cluster_", clusterName, "_plots.pdf")))

      c(..., nuPlot, zoomedNuPlot) %<-%
        cleanPlots(objCOTAN, includePCA = FALSE)

      genesToLabel <-
        head(rownames(GDIData[order(GDIData[["GDI"]],
                                    decreasing = TRUE), ]), n = 10L)
      gdiPlot <- GDIPlot(objCOTAN, GDIIn = GDIData, GDIThreshold = GDIThreshold,
                         genes = list("top 10 GDI genes" = genesToLabel))

      plot(nuPlot)
      plot(zoomedNuPlot)
      plot(gdiPlot)

      rm(nuPlot, zoomedNuPlot, gdiPlot)
      dev.off()
    },
    error = function(err) {
      logThis(paste("While saving cluster plots", err),
              logLevel = 0L)
    }
  )

  rm(objCOTAN)
  gc()

  # A cluster is deemed uniform if the number of genes
  # with [GDI > GDIThreshold] is at the most ratioAboveThreshold
  gdi <- getColumnFromDF(GDIData, "GDI")

  quantAboveThr <- quantile(gdi, probs = 1.0 - ratioAboveThreshold)
  percAboveThr <- sum(gdi >= GDIThreshold) / length(gdi)

  clusterIsUniform <- isClusterUniform(GDIThreshold, ratioAboveThreshold,
                                       quantAboveThr, percAboveThr,
                                       GDIThreshold, ratioAboveThreshold)

  logThis(paste0(
    "Cluster ", clusterName, ", with size ", clusterSize, ", is ",
    (if (clusterIsUniform) {""} else {"not "}), "uniform\n",
    round(percAboveThr * 100.0, digits = 2L), "% of the genes is above ",
    "the given GDI threshold ", GDIThreshold, "\n",
    "highest ", round(ratioAboveThreshold * 100.0, digits = 2L),
    "% GDI quantile is at ", round(quantAboveThr, digits = 4L)), logLevel = 3L)

  if (isTRUE(saveObj)) tryCatch({
      pre <- ""
      if (!clusterIsUniform) {
        pre <- "non-"
      }
      outFile <- file.path(outDir,
                           paste0(pre, "uniform_cluster_", clusterName, ".csv"))
      write.csv(cells, file = outFile)
    },
    error = function(err) {
      logThis(paste("While saving current clusterization", err),
              logLevel = 0L)
    }
  )

  return(list("isUniform" = clusterIsUniform,
              "fractionAbove" = percAboveThr,
              "ratioQuantile" = quantAboveThr[[1L]],
              "size" = clusterSize,
              "GDIThreshold" = GDIThreshold,
              "ratioAboveThreshold" = ratioAboveThreshold))
}
