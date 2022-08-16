#' automatic.COTAN.object.creation
#'
#' @param df dataframe with the row counts
#' @param out_dir directory for the output
#' @param GEO GEO or other code that identify the dataset
#' @param sc.method Type of single cell RNA-seq method used
#' @param cond A string that will identify the sample or condition. It will be part of the
#' final file name.
#' @param cores number of cores to be used
#' @param mt A boolean (default F). If T mitochondrial  genes will be kept in the analysis,
#' otherwise they will be removed.
#' @param mt_prefix is the prefix that identify the mitochondrial  genes (default is the mouse
#' prefix: "^mt")
#' @return It return the COTAN object. It will also store it directly in the output directory
#' @export
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom utils write.csv
#' @importFrom methods new
#' @import grDevices
#' @import ggplot2
#' @import ggrepel
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
setGeneric("automatic.COTAN.object.creation", function(df, out_dir, GEO, sc.method,
                                                       cond, mt = FALSE, mt_prefix = "^mt", cores = 1) {
  standardGeneric("automatic.COTAN.object.creation")
})
#' @rdname automatic.COTAN.object.creation
setMethod(
  "automatic.COTAN.object.creation", "data.frame",
  function(df, out_dir, GEO, sc.method, cond, mt = FALSE, mt_prefix = "^mt", cores = 1) {
    start_time_all <- Sys.time()

    means <- PC1 <- PC2 <- nu <- NULL

    mycolours <- c("A" = "#8491B4B2", "B" = "#E64B35FF")
    my_theme <- theme(
      axis.text.x = element_text(
        size = 14, angle = 0, hjust = .5,
        vjust = .5,
        face = "plain", colour = "#3C5488FF"
      ),
      axis.text.y = element_text(
        size = 14, angle = 0, hjust = 0,
        vjust = .5,
        face = "plain", colour = "#3C5488FF"
      ),
      axis.title.x = element_text(
        size = 14, angle = 0, hjust = .5,
        vjust = 0,
        face = "plain", colour = "#3C5488FF"
      ),
      axis.title.y = element_text(
        size = 14, angle = 90, hjust = .5,
        vjust = .5,
        face = "plain", colour = "#3C5488FF"
      )
    )
    obj <- methods::new("scCOTAN", raw = df)
    obj <- initRaw(obj, GEO = GEO, sc.method = sc.method, cond = cond)
    if (mt == FALSE) {
      genes_to_rem <- rownames(obj@raw[grep(mt_prefix, rownames(obj@raw)), ])
      obj@raw <- obj@raw[!rownames(obj@raw) %in% genes_to_rem, ]
      cells_to_rem <- colnames(obj@raw[which(colSums(obj@raw) == 0)])
      obj@raw <- obj@raw[, !colnames(obj@raw) %in% cells_to_rem]
    }
    t <- cond

    print(paste("Condition ", t, sep = ""))
    #--------------------------------------
    n_cells <- length(colnames(obj@raw))
    print(paste("n cells", n_cells, sep = " "))

    n_it <- 1

    if (!file.exists(out_dir)) {
      dir.create(file.path(out_dir))
    }

    if (!file.exists(paste(out_dir, "cleaning", sep = ""))) {
      dir.create(file.path(out_dir, "cleaning"))
    }

    ttm <- clean(obj)

    obj <- ttm$object

    pdf(file.path(out_dir, "cleaning", paste(t, "_", n_it, "_plots_without_cleaning.pdf",
      sep = ""
    )))
    plot(ttm$pca.cell.2)
    plot.fig <- ggplot(ttm$D, aes(x = n, y = means)) +
      geom_point() +
      geom_text_repel(
        data = subset(ttm$D, n > (max(ttm$D$n) - 15)),
        aes(n, means, label = rownames(ttm$D[ttm$D$n > (max(ttm$D$n) - 15), ])),
        nudge_y = 0.05,
        nudge_x = 0.05,
        direction = "x",
        angle = 90,
        vjust = 0,
        segment.size = 0.2
      ) +
      ggtitle(
        label = "B cell group genes mean expression", subtitle =
          " - B group NOT removed -"
      ) +
      my_theme +
      theme(
        plot.title = element_text(
          color = "#3C5488FF",
          size = 20, face = "italic",
          vjust = -10, hjust = 0.02
        ) # ,
        # plot_subtitle = element_text(color = "darkred",vjust = - 15,hjust = 0.01 )
      )
    plot(plot.fig)

    dev.off()

    nu_est <- round(obj@nu, digits = 7)

    plot_nu <- ggplot(ttm$pca_cells, aes(x = PC1, y = PC2, colour = log(nu_est)))

    plot_nu <- plot_nu + geom_point(size = 1, alpha = 0.8) +
      scale_color_gradient2(
        low = "#E64B35B2", mid = "#4DBBD5B2", high = "#3C5488B2",
        midpoint = log(mean(nu_est)), name = "ln(nu)"
      ) +
      ggtitle("Cells PCA coloured by cells efficiency") +
      my_theme + theme(
        plot.title = element_text(color = "#3C5488FF", size = 20),
        legend.title = element_text(
          color = "#3C5488FF", size = 14,
          face = "italic"
        ),
        legend.text = element_text(color = "#3C5488FF", size = 11),
        legend.key.width = unit(2, "mm"),
        legend.position = "right"
      )

    pdf(file.path(out_dir, "cleaning", paste(t, "_plots_PCA_efficiency_colored.pdf", sep = "")))
    plot(plot_nu)
    dev.off()

    nu_df <- data.frame("nu" = sort(obj@nu), "n" = seq_along(obj@nu))

    pdf(file.path(out_dir, "cleaning", paste(t, "_plots_efficiency.pdf", sep = "")))
    plot(ggplot(nu_df, aes(x = n, y = nu)) +
      geom_point(colour = "#8491B4B2", size = 1) +
      my_theme +
      annotate(geom = "text", x = 50, y = 0.25, label = "nothing to remove ", color = "darkred"))
    dev.off()

    analysis_time <- Sys.time()

    print("Cotan analysis function started")
    obj <- cotan_analysis(obj, cores = cores)

    coex_time <- Sys.time()
    analysis_time <- difftime(Sys.time(), analysis_time, units = "mins")

    print(paste0("Only analysis time ", analysis_time))

    print("Cotan coex estimation started")
    obj <- get.coex(obj)

    end_time <- Sys.time()

    all.time <- difftime(end_time, start_time_all, units = "mins")
    print(paste0("Total time ", all.time))

    coex_time <- difftime(end_time, coex_time, units = "mins")

    print(paste0("Only coex time ", coex_time))

    utils::write.csv(data.frame(
      "type" = c("tot_time", "analysis_time", "coex_time"),
      "times" =
        c(
          as.numeric(all.time),
          as.numeric(analysis_time),
          as.numeric(coex_time)
        ),
      "n.cells" = n_cells, "n.genes" = dim(obj@raw)[1]
    ),
    file = file.path(out_dir, paste(t, "_times.csv", sep = ""))
    )

    print(paste0("Saving elaborated data locally at ", out_dir, t, ".cotan.RDS"))
    saveRDS(obj, file = file.path(out_dir, paste(t, ".cotan.RDS", sep = "")))

    return(obj)
  }
)
