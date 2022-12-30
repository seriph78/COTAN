#' checkClusterUniformity
#'
#' @description This function takes a `COTAN` object and a cells' cluster and
#'   checks whether the latter is **uniform** by GDI.
#'
#' @details It runs `COTAN` to check whether the GDI is lower than \eqn{1.5} for
#'   the \eqn{99\%} of the genes. If the GDI results to be too high, the cluster
#'   is deemed **non-uniform**
#'
#' @param objCOTAN a `COTAN` object
#' @param cluster the tag of the cluster
#' @param cells the cells in the cluster
#' @param cores number of cores used
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file
#' @param outDir an existing directory for the analysis output.
#'
#' @return an array of cells that need to be re-clustered or nothing
#'
#' @importFrom utils head
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- automaticCOTANObjectCreation(raw = raw.dataset,
#'                                          GEO = "S",
#'                                          sequencingMethod = "10X",
#'                                          sampleCondition = "Test",
#'                                          cores = 12,
#'                                          saveObj = FALSE,
#'                                          outDir = tempdir())
#' clusters <- cellsUniformClustering(objCOTAN, cores = 12,
#'                                    saveObj = FALSE,
#'                                    outDir = tempdir())
#' isUniform <- checkClusterUniformity(objCOTAN, cluster = clusters[1],
#'                                     cells = getCells(objCOTAN)[clusters %in% clusters[1]],
#'                                     cores = 12, saveObj = FALSE, outDir = tempdir())
#'
#' @rdname checkClusterUniformity
#'

checkClusterUniformity <- function(objCOTAN, cluster, cells,
                                   cores = 1, saveObj = TRUE, outDir = ".") {

  cellsToDrop <- getCells(objCOTAN)[!getCells(objCOTAN) %in% cells]

  objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)

  objCOTAN <- proceedToCoex(objCOTAN, cores = cores, saveObj = FALSE)
  gc()

  logThis(paste0("Checking uniformity for the cluster '", cluster,
                 "' with ", getNumCells(objCOTAN), " cells"), logLevel = 2)

  GDIData <- calculateGDI(objCOTAN)

  # Plots
  if (saveObj) {
    prevOptVale <- options(ggrepel.max.overlaps = Inf)

    pdf(file.path(outDir, paste0("cluster_", cluster, "_plots.pdf")))

    plots <- cleanPlots(objCOTAN)

    plot(plots[["pcaCells"]])
    plot(plots[["genes"]])
    plot(plots[["UDE"]])
    plot(plots[["nu"]])
    rm(plots)

    genesToLabel = head(rownames(GDIData[order(GDIData[["GDI"]], decreasing = T), ]), n = 10)
    plot(GDIPlot(objCOTAN, GDI.df = GDIData, genes = list("top 10 GDI genes" = genesToLabel)))

    dev.off()
    options(ggrepel.max.overlaps = prevOptVale)
  }

  rm(objCOTAN)
  gc()

  # A cluster is deemed uniform if the number of genes
  # with [GDI > 1.5] is not more than 1%
  clusterIsUniform <- (nrow(GDIData[GDIData[["GDI"]] >= 1.5, ]) <=
                         0.01 * nrow(GDIData))

  if (!clusterIsUniform && saveObj) {
    outFile <- file.path(outDir, paste0("non-uniform_cluster_", cluster, ".csv"))
    write.csv(cells, file = outFile)
  }

  return(invisible(clusterIsUniform))
}


