
#'
#' @details `checkClusterUniformity()` takes a `COTAN` object and a cells'
#'   *cluster* and checks whether the latter is **uniform** by `GDI`. The
#'   function runs `COTAN` to check whether the `GDI` is lower than the given
#'   `GDIThreshold` for the \eqn{99\%} of the genes. If the `GDI` results to be
#'   too high for too many genes, the *cluster* is deemed **non-uniform**.
#'
#' @param objCOTAN a `COTAN` object
#' @param clusterName the tag of the *cluster*
#' @param cells the cells belonging to the *cluster*
#' @param GDIThreshold the threshold level that discriminates uniform
#'   *clusters*. It defaults to \eqn{1.43}
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
#'   * `"firstPercentile"`: the quantile associated to the highest percentile
#'   * `"size"`: the number of cells in the cluster
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

checkClusterUniformity <- function(objCOTAN, clusterName, cells,
                                   GDIThreshold = 1.43, cores = 1L,
                                   optimizeForSpeed = TRUE, deviceStr = "cuda",
                                   saveObj = TRUE, outDir = ".") {

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
  # with [GDI > GDIThreshold] is not more than 1%
  gdi <- GDIData[["GDI"]]

  quantAboveThr <- quantile(gdi, probs = 0.99)
  percAboveThr <- sum(gdi >= GDIThreshold) / length(gdi)

  clusterIsUniform <- percAboveThr <= 0.01

  logThis(paste0("Cluster ", clusterName, ", with size ", clusterSize, ", is ",
                 (if (clusterIsUniform) {""} else {"not "}), "uniform\n",
                 round(percAboveThr * 100.0, digits = 2L),
                 "% of the genes is above the given GDI threshold ",
                 GDIThreshold, "\n", "GDI 99% quantile is at ",
                 round(quantAboveThr, digits = 4L)), logLevel = 3L)

  if (isTRUE(saveObj)) tryCatch({
      pre <- ""
      if(!clusterIsUniform) {
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
              "firstPercentile" = quantAboveThr[[1L]],
              "size" = clusterSize))
}
