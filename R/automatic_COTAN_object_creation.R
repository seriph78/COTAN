#' automatic.COTAN.object.creation
#'
#' @param df dataframe with the row counts
#' @param out_dir directory for the output
#' @param save.obj the created object is not automatically writen on disk (default NO).
#' If "Yes" the object is written in out_dir
#' @param GEO GEO or other code that identify the dataset
#' @param sc.method Type of single cell RNA-seq method used
#' @param cond A string that will identify the sample or condition. It will be part of the
#' final file name.
#' @param cores number of cores to be used
#' @return It return the COTAN object. It will also store it directly in the output directory
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
#'
#' data("raw.dataset")
#' obj <- automatic.COTAN.object.creation(
#'   df = raw.dataset,
#'   out_dir = tempdir(),
#'   GEO = "test_GEO",
#'   sc.method = "test_method",
#'   cond = "test"
#' )
#'
#' @rdname automatic.COTAN.object.creation

automatic.COTAN.object.creation <-
  function(df, out_dir, save.obj = "NO",
           GEO, sc.method, cond,
           #mt = FALSE, mt_prefix="^mt",
           cores = 1) {
    start_time_all <- Sys.time()

    obj <- COTAN(raw = df)
    obj <- initializeMetaDataset(obj, GEO = GEO,
                                 sequencingMethod = sc.method,
                                 sampleCondition = cond)

    #if (isFALSE(mt)) {
    #  genes_to_rem <- getGenes(obj)[grep(mt_prefix, getGenes(obj)),])
    #  cells_to_rem <- getCells(obj)[which(getCellsSize(obj) == 0)])
    #  dropGenesCells(obj, genes_to_rem, cells_to_rem)
    #}

    print(paste0("Condition ", cond))
    print(paste("n cells", getNumCells(obj)))


    #--------------------------------------
    if (!file.exists(out_dir)) {
      dir.create(file.path(out_dir))
    }

    if (!file.exists(paste0(out_dir, "cleaning"))) {
      dir.create(file.path(out_dir, "cleaning"))
    }

    list[obj, data] <- clean(obj)
    plots <- cleanPlots(obj, data[["pcaCells"]], data[["D"]])

    means <- PC1 <- PC2 <- nu <- NULL

    {
      numIter <- 1
      pdf(file.path(out_dir, "cleaning",
                    paste0(cond, "_", numIter, "_plots_without_cleaning.pdf")))
      plot(plots[["pcaCells"]])
      plot(plots[["genes"]])
      dev.off()
    }

    {
      pdf(file.path(out_dir, "cleaning",
                    paste0(cond,"_plots_PCA_efficiency_colored.pdf")))
      plot(plots[["UDE"]])
      dev.off()
    }

    {
      pdf(file.path(out_dir, "cleaning",
                    paste0(cond,"_plots_efficiency.pdf")))
      plot(plots[["nu"]] +
           annotate(geom = "text", x = 50, y = 0.25,
                    label = "nothing to remove ", color = "darkred"))
      dev.off()
    }

    rm(plots)
    gc()

    print("Cotan analysis function started")
    analysis_time <- Sys.time()

    obj <- estimateDispersion(obj, cores = cores)

    genes_coex_time <- Sys.time()
    analysis_time <- difftime(genes_coex_time, analysis_time, units = "mins")

    print(paste0("Only analysis time ", analysis_time))

    print("Cotan genes' coex estimation started")
    obj <- calculateCoex(obj)

    end_time <- Sys.time()

    all.time <- difftime(end_time, start_time_all, units = "mins")
    print(paste0("Total time ", all.time))

    genes_coex_time <- difftime(end_time, genes_coex_time, units = "mins")

    print(paste0("Only genes' coex time ", genes_coex_time))

    utils::write.csv(data.frame("type" = c("tot_time",
                                           "analysis_time",
                                           "genes_coex_time"),
                                "times"= c(as.numeric(all.time),
                                           as.numeric(analysis_time),
                                           as.numeric(genes_coex_time) ),
                                "n.cells"=getNumCells(obj),
                                "n.genes"=getNumGenes(obj) ),
                     file = file.path(out_dir, paste0(cond, "_times.csv")))

              obj <- as(obj, "scCOTAN")

              if (save.obj == "yes" | save.obj == "Yes" | save.obj == "YES") {
                print(paste0("Saving elaborated data locally at ",
                             out_dir, cond, ".cotan.RDS"))
                saveRDS(obj,file = file.path(out_dir,paste0(cond, ".cotan.RDS")))
              }
              return(obj)

  }

