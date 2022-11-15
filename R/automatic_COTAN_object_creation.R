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
#' @rdname automatic.COTAN.object.creation
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
setGeneric(
  "automatic.COTAN.object.creation",
  function(df, out_dir,save.obj = "NO",
           GEO, sc.method, cond, 
           #mt = FALSE, mt_prefix="^mt", 
           cores = 1) {
    standardGeneric("automatic.COTAN.object.creation")
  }
)
#' @rdname automatic.COTAN.object.creation
setMethod(
  "automatic.COTAN.object.creation",
  "data.frame",
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

              ttm <- clean(obj)

              means <- PC1 <- PC2 <- nu <- NULL
              
              obj <- ttm$object

              n_it <- 1

              pdf(file.path(out_dir, "cleaning",
                            paste0(cond, "_", n_it, "_plots_without_cleaning.pdf")))
              plot(ttm$pca.cell.2)
              plot.fig = ggplot(ttm$D, aes(x = n, y = means)) +
                         geom_point() +
                         geom_text_repel(
                           data = subset(ttm$D, n > (max(ttm$D$n)- 15)),
                           aes(n, means, label = rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
                           nudge_y = 0.05, nudge_x = 0.05, direction = "x",
                           angle = 90, vjust = 0, segment.size = 0.2) +
                         ggtitle(label = "B cell group genes mean expression",
                                 subtitle = " - B group NOT removed -") +
                         plotTheme("genes")
              plot(plot.fig)

              dev.off()

              nu_est <- round(obj@nu, digits = 7)

              plot_nu <- ggplot(ttm$pca_cells,
                                aes(x = PC1,y = PC2, colour = log(nu_est))) +
                         geom_point(size = 1, alpha=  0.8) +
                         scale_color_gradient2(low = "#E64B35B2", mid = "#4DBBD5B2", high = "#3C5488B2",
                                               midpoint = log(mean(nu_est)), name = "ln(nu)") +
                         ggtitle("Cells PCA coloured by cells efficiency") +
                         plotTheme("UDE")

              pdf(file.path(out_dir, "cleaning",
                            paste0(cond,"_plots_PCA_efficiency_colored.pdf")))
              plot(plot_nu)
              dev.off()

              nu_df <- data.frame("nu"= sort(obj@nu), "n"=seq_along(obj@nu))

              pdf(file.path(out_dir, "cleaning",
                            paste0(cond,"_plots_efficiency.pdf")))
              plot(ggplot(nu_df, aes(x = n, y = nu)) +
                   geom_point(colour = "#8491B4B2", size = 1) +
                   annotate(geom = "text", x = 50, y = 0.25,
                            label = "nothing to remove ", color = "darkred") +
                   plotTheme("common") )
              dev.off()

              analysis_time <- Sys.time()

              print("Cotan analysis function started")
              obj <- cotan_analysis(obj,cores = cores)

              coex_time <- Sys.time()
              analysis_time <- difftime(Sys.time(), analysis_time, units = "mins")

              print(paste0("Only analysis time ", analysis_time))

              print("Cotan coex estimation started")
              obj <- get.coex(obj)

              end_time <- Sys.time()

              all.time <-  difftime(end_time, start_time_all, units = "mins")
              print(paste0("Total time ", all.time))

              coex_time <- difftime(end_time, coex_time, units = "mins")

              print(paste0("Only coex time ",coex_time))

              utils::write.csv(data.frame("type" = c("tot_time",
                                                     "analysis_time",
                                                     "coex_time"),
                                          "times"= c(as.numeric(all.time),
                                                     as.numeric(analysis_time),
                                                     as.numeric(coex_time) ),
                                          "n.cells"=getNumCells(obj),
                                          "n.genes"=getNumGenes(obj) ),
                               file = file.path(out_dir, paste0(cond, "_times.csv")))

              
              if (save.obj == "yes" | save.obj == "Yes" | save.obj == "YES") {
                print(paste0("Saving elaborated data locally at ",
                             out_dir, cond, ".cotan.RDS"))
                saveRDS(obj,file = file.path(out_dir,paste0(cond, ".cotan.RDS")))
              }
              return(obj)

          }

)
