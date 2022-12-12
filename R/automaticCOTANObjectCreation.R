#' automaticCOTANObjectCreation
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
#' @returns It return the COTAN object. It will also store it directly in the
#'   output directory
#'
#' @export
#'
#' @import grDevices
#' @import gsubfn
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom utils write.csv
#' @importFrom methods new
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' data("raw.dataset")
#' obj <- automaticCOTANObjectCreation(
#'   raw = raw.dataset,
#'   GEO = "test_GEO",
#'   sequencingMethod = "test_method",
#'   sampleCondition = "test",
#'   saveObj = TRUE,
#'   outDir = tempdir(),
#'   cores = 12)
#'
#' @rdname automaticCOTANObjectCreation

automaticCOTANObjectCreation <-
  function(raw, GEO, sequencingMethod, sampleCondition,
           cores = 1, saveObj = TRUE, outDir = ".") {
    start_time_all <- Sys.time()

    obj <- COTAN(raw = raw)
    obj <- initializeMetaDataset(obj, GEO = GEO,
                                 sequencingMethod = sequencingMethod,
                                 sampleCondition = sampleCondition)

    #if (isFALSE(mt)) {
    #  genes_to_rem <- getGenes(obj)[grep(mt_prefix, getGenes(obj)),])
    #  cells_to_rem <- getCells(obj)[which(getCellsSize(obj) == 0)])
    #  dropGenesCells(obj, genes_to_rem, cells_to_rem)
    #}

    logThis(paste0("Condition ", sampleCondition), logLevel = 2)
    logThis(paste("n cells", getNumCells(obj)), logLevel = 2)


    #--------------------------------------
    logThis("Cotan analysis functions started", logLevel = 1)

    list[obj, data] <- clean(obj, calcExtraData = saveObj)

    if (isTRUE(saveObj)) {
      plots <- cleanPlots(obj, data[["pcaCells"]], data[["D"]])

      if (!file.exists(outDir)) {
        dir.create(file.path(outDir))
      }

      if (!file.exists(paste0(outDir, "cleaning"))) {
        dir.create(file.path(outDir, "cleaning"))
      }
      {
        numIter <- 1
        pdf(file.path(outDir, "cleaning",
                      paste0(sampleCondition, "_", numIter, "_plots_without_cleaning.pdf")))
        plot(plots[["pcaCells"]])
        plot(plots[["genes"]])
        dev.off()
      }

      {
        pdf(file.path(outDir, "cleaning",
                      paste0(sampleCondition, "_plots_PCA_efficiency_colored.pdf")))
        plot(plots[["UDE"]])
        dev.off()
      }

      {
        pdf(file.path(outDir, "cleaning",
                      paste0(sampleCondition,"_plots_efficiency.pdf")))
        plot(plots[["nu"]] +
             annotate(geom = "text", x = 50, y = 0.25,
                      label = "nothing to remove ", color = "darkred"))
        dev.off()
      }

      rm(plots)
      gc()
    }

    analysis_time <- Sys.time()

    obj <- estimateDispersionBisection(obj, cores = cores)

    genes_coex_time <- Sys.time()
    analysis_time <- difftime(genes_coex_time, analysis_time, units = "mins")

    logThis(paste0("Only analysis time ", analysis_time), logLevel = 3)

    logThis("Cotan genes' coex estimation started", logLevel = 1)
    obj <- calculateCoex(obj)

    end_time <- Sys.time()

    all.time <- difftime(end_time, start_time_all, units = "mins")
    logThis(paste0("Total time ", all.time), logLevel = 3)

    genes_coex_time <- difftime(end_time, genes_coex_time, units = "mins")

    logThis(paste0("Only genes' coex time ", genes_coex_time), logLevel = 3)

    if (saveObj) {
      utils::write.csv(data.frame("type" = c("tot_time",
                                             "analysis_time",
                                             "genes_coex_time"),
                                  "times"= c(as.numeric(all.time),
                                             as.numeric(analysis_time),
                                             as.numeric(genes_coex_time) ),
                                  "n.cells"=getNumCells(obj),
                                  "n.genes"=getNumGenes(obj) ),
                       file = file.path(outDir, paste0(sampleCondition, "_times.csv")))

      logThis(paste0("Saving elaborated data locally at ",
                     outDir, sampleCondition, ".cotan.RDS"),
              logLevel = 1)
      saveRDS(obj, file = file.path(outDir, paste0(sampleCondition, ".cotan.RDS")))
    }

    return(obj)
  }

