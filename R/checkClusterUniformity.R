
#' @details `checkClusterUniformity()` takes a `COTAN` object and a cells'
#'   *cluster* and checks whether the latter is **uniform** by `GDI`. The
#'   function runs `COTAN` to check whether the `GDI` is lower than the given
#'   `GDIThreshold` for the \eqn{99\%} of the genes. If the `GDI` results to be
#'   too high for too many genes, the *cluster* is deemed **non-uniform**.
#'
#' @param objCOTAN a `COTAN` object
#' @param cluster the tag of the *cluster*
#' @param cells the cells belonging to the *cluster*
#' @param GDIThreshold the threshold level that discriminates uniform
#'   *clusters*. It defaults to \eqn{1.4}
#' @param cores number of cores used
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file(s)
#' @param outDir an existing directory for the analysis output. The effective
#'   output will be paced in a sub-folder.
#'
#' @returns `checkClusterUniformity` returns a list with:
#'   * `"isUniform"`: a flag indicating whether the *cluster* is **uniform**
#'   * `"fractionAbove"`: the percentage of genes with `GDI` above the threshold
#'   * `"1stPercentile"`: the quantile associated to the highest percentile
#'
#' @importFrom utils head
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @importFrom withr local_options
#'
#' @export
#'
#' @rdname UniformClusters
#'

checkClusterUniformity <- function(objCOTAN, cluster, cells,
                                   GDIThreshold = 1.4, cores = 1L,
                                   saveObj = TRUE, outDir = ".") {

  cellsToDrop <- getCells(objCOTAN)[!getCells(objCOTAN) %in% cells]

  objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)

  objCOTAN <- proceedToCoex(objCOTAN, cores = cores, saveObj = FALSE)
  gc()

  logThis(paste0("Checking uniformity for the cluster '", cluster,
                 "' with ", getNumCells(objCOTAN), " cells"), logLevel = 2L)

  GDIData <- calculateGDI(objCOTAN)

  # Plots
  if (saveObj) {
    # this will be restored to previous value on function exit
    local_options(list(ggrepel.max.overlaps = Inf))

    pdf(file.path(outDir, paste0("cluster_", cluster, "_plots.pdf")))

    plots <- cleanPlots(objCOTAN)

    plot(plots[["pcaCells"]])
    plot(plots[["genes"]])
    plot(plots[["UDE"]])
    plot(plots[["nu"]])
    rm(plots)

    genesToLabel <- head(rownames(GDIData[order(GDIData[["GDI"]],
                                                decreasing = TRUE), ]), n = 10L)
    plot(GDIPlot(objCOTAN, GDIIn = GDIData, GDIThreshold = GDIThreshold,
                 genes = list("top 10 GDI genes" = genesToLabel)))

    dev.off()
  }

  rm(objCOTAN)
  gc()

  # A cluster is deemed uniform if the number of genes
  # with [GDI > GDIThreshold] is not more than 1%
  gdi <- GDIData[["GDI"]]

  quantAboveThr <- quantile(gdi, probs = 0.99)
  percAboveThr <- sum(gdi >= GDIThreshold) / length(gdi)

  clusterIsUniform <- percAboveThr <= 0.01

  logThis(paste0("Cluster ", cluster, " is ",
                 (if(clusterIsUniform) {""} else {"not"}), " uniform\n",
                 round(percAboveThr * 100.0, digits = 2L),
                 "% of the genes is above the given GDI threshold ",
                 GDIThreshold, "\n", "GDI 99% quantile is at ",
                 round(quantAboveThr, digits = 4L), logLevel = 3L)

  if (!clusterIsUniform && saveObj) {
    outFile <- file.path(outDir,
                         paste0("non-uniform_cluster_", cluster, ".csv"))
    write.csv(cells, file = outFile)
  }

  return(list("isUniform" = clusterIsUniform,
              "fractionAbove" = percAboveThr,
              "1stPercentile" = quantAboveThr))
}
