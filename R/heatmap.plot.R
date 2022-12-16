## Main functions


#' plot_heatmap
#'
#' This is the function that create the heatmap of one or more COTAN object.
#'
#' @param p_val.tr p-value threshold. Default is 0.05
#' @param df_genes this is a list of gene array. The first array will define genes in the columns.
#' @param sets This is a numeric array indicating from which fields of the previous list
#' will be considered
#' @param conditions An array of prefixes indicating the different files.
#' @param dir The directory in which are all COTAN files (corresponding to the previous prefixes)
#'
#' @return a ggplot object
#'
#' @importFrom Matrix forceSymmetric
#' @importFrom tidyr pivot_longer
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_fill_gradient2
#'
#' @import scales
#'
#' @export
#'
#' @examples
#' data("ERCCraw")
#' rownames(ERCCraw) = ERCCraw$V1
#' ERCCraw = ERCCraw[,2:ncol(ERCCraw)]
#' objCOTAN <- COTAN(raw = ERCCraw)
#' data_dir <- tempdir()
#' saveRDS(objCOTAN, file = file.path(data_dir, "ERCC.cotan.RDS"))
#' # some genes
#' primary.markers <- c("ERCC-00154", "ERCC-00156", "ERCC-00164")
#' # a example of named list of different gene set
#' gene.sets.list <- list(
#'   "primary.markers" = primary.markers,
#'   "2.R" = c("ERCC-00170", "ERCC-00158"),
#'   "3.S" = c("ERCC-00160", "ERCC-00162")
#' )
#' plot_heatmap(
#'   p_v = 0.05, df_genes = gene.sets.list,
#'   sets = c(2, 3), conditions = c("ERCC"), dir = paste0(data_dir, "/")
#' )
#'
#' @rdname plot_heatmap
setGeneric("plot_heatmap", function(p_val.tr = 0.05, df_genes, sets, conditions, dir) {
  standardGeneric("plot_heatmap")
})
setMethod(
  "plot_heatmap", "ANY",
  function(p_val.tr = 0.05, df_genes, sets, conditions, dir) {
    time <- g2 <- NULL
    print("plot heatmap")
    gr <- df_genes[[1]]
    ge <- unique(array(sort(unlist(df_genes[sets]))))
    df.to.print <- data.frame()
    for (ET in conditions) {
      print(paste0("Loading condition", ET))
      obj <- readRDS(paste0(dir, ET, ".cotan.RDS"))
      obj <- as(obj, "scCOTAN")

      if (any(gr %in% rownames(obj@coex)) == FALSE) {
        print(paste0("primary markers all absent in ", ET))
        stop()
      }
      p_val <- calculatePValue(obj, geneSubsetCol = gr, geneSubsetRow = ge)
      p_val <- as.data.frame(p_val)

      # this to add some eventually effective housekeeping genes
      if (any(ge %in% obj@hk)) {
        genes.to.add <- ge[ge %in% obj@hk]
        temp.hk.rows <- as.data.frame(matrix(
                          ncol = ncol(p_val),
                          nrow = length(genes.to.add)
                        ))
        rownames(temp.hk.rows) <- genes.to.add
        colnames(temp.hk.rows) <- colnames(p_val)
        temp.hk.rows <- 1
        p_val <- rbind(p_val, temp.hk.rows)
      }

      if (any(gr %in% obj@hk)) {
        genes.to.add <- gr[gr %in% obj@hk]
        temp.hk.cols <- as.data.frame(matrix(
                          ncol = length(genes.to.add),
                          nrow = nrow(p_val)
                        ))
        colnames(temp.hk.cols) <- genes.to.add
        rownames(temp.hk.cols) <- rownames(p_val)
        temp.hk.cols <- 1
        p_val <- cbind(p_val, temp.hk.cols)
      }

      p_val$g2 <- as.vector(rownames(p_val))
      df.temp.pval <- pivot_longer(p_val, cols = seq_along(colnames(p_val)) - 1, names_to = "g1", values_to = "p_val")

      coex <- getGenesCoex(obj, genes = gr)
      diag(coex) <- 0
      coex <- coex[rownames(coex) %in% ge, ]
      # this to add some eventually effective housekeeping genes
      if (any(ge %in% obj@hk)) {
        temp.hk.rows <- 0
        coex <- rbind(coex, temp.hk.rows)
      }

      if (any(gr %in% obj@hk)) {
        temp.hk.cols <- 0
        coex <- cbind(coex, temp.hk.cols)
      }
      #---------------------------------------------------------
      coex <- as.data.frame(coex)
      coex$g2 <- as.vector(rownames(coex))
      df.temp.coex <- pivot_longer(coex, cols = seq_along(colnames(p_val)) - 1, names_to = "g1", values_to = "coex")
      df.temp <- merge(df.temp.coex, df.temp.pval)
      df.temp$time <- ET
      df.temp$type <- NA
      df.temp$absent <- NA
      df.temp2 <- data.frame()
      for (type in names(df_genes)[sets]) {
        for (g1 in gr) {
          tt <- df.temp[df.temp$g2 %in% df_genes[[type]] & df.temp$g1 == g1, ]
          # control if the subset is smaller than the number of wanted genes
          if (dim(tt)[1] < length(df_genes[[type]])) {
            n.row <- length(df_genes[[type]]) - dim(tt)[1]
            t.rows <- as.data.frame(matrix(nrow = n.row, ncol = 7))
            colnames(t.rows) <- colnames(tt)
            t.rows[, "g1"] <- g1
            t.rows[, "time"] <- ET
            t.rows[, "absent"] <- "yes"
            t.rows[, "p_val"] <- 1
            t.rows[, "g2"] <- df_genes[[type]][!df_genes[[type]] %in% tt$g2]
            tt <- rbind(tt, t.rows)
          }
          tt$type <- type
          df.temp2 <- rbind(df.temp2, tt)
        }
        print(type)
      }
      df.temp <- df.temp2
      df.temp$t_hk <- ifelse((df.temp$g2 %in% obj@hk) | (df.temp$g1 %in% obj@hk), "hk", "n")
      df.temp[df.temp$p_val > p_val.tr, ]$coex <- 0
      df.to.print <- rbind(df.to.print, df.temp)
    }
    print(paste("min coex:", min(df.to.print$coex, na.rm = TRUE),
                "max coex",  max(df.to.print$coex, na.rm = TRUE)))

    heatmap <- ggplot(data = subset(df.to.print, type %in% names(df_genes)[sets]),
                      aes(time, factor(g2, levels = rev(levels(factor(g2)))))) +
               geom_tile(aes(fill = coex), colour = "black", show.legend = TRUE) +
               facet_grid(type ~ g1, scales = "free", space = "free") +
               scale_fill_gradient2(low = "#E64B35FF", mid = "gray93", high = "#3C5488FF",
                                    midpoint = 0,na.value = "grey80", space = "Lab",
                                    guide = "colourbar", aesthetics = "fill", oob = scales::squish) +
               plotTheme("heatmap", textSize = 9)

    return(heatmap)
  }
)


#' plot_general.heatmap
#'
#' This function is used to plot an heatmap made using only some genes, as markers,
#' and collecting all other genes correlated
#' with these markers with a p-value smaller than the set threshold. Than all relations are plotted.
#' Primary markers will be plotted as groups of rows. Markers list will be plotted as columns.
#' @param prim.markers A set of genes plotted as rows.
#' @param markers.list A set of genes plotted as columns.
#' @param dir The directory where the COTAN object is stored.
#' @param condition The prefix for the COTAN object file.
#' @param p_value The p-value threshold
#' @param symmetric A boolean: default F. If T the union of prim.markers and marker.list
#' is sets as both rows and column genes
#'
#' @return A ggplot2 object
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @importFrom  stats quantile
#' @importFrom circlize colorRamp2
#' @importFrom Matrix forceSymmetric
#' @importFrom rlang is_empty
#' @export
#' @rdname plot_general.heatmap
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' data_dir <- tempdir()
#' saveRDS(objCOTAN, file = file.path(data_dir, "ERCC.cotan.RDS"))
#' # some genes
#' primary.markers <- c("ERCC-00154", "ERCC-00156", "ERCC-00164")
#' # a example of named list of different gene set
#' gene.sets.list <- list(
#'   "primary.markers" = primary.markers,
#'   "2.R" = c("ERCC-00170", "ERCC-00158"),
#'   "3.S" = c("ERCC-00160", "ERCC-00162")
#' )
#' plot_general.heatmap(
#'   prim.markers = primary.markers, p_value = 0.05, markers.list = gene.sets.list,
#'   condition = "ERCC", dir = paste0(data_dir, "/")
#' )
setGeneric("plot_general.heatmap", function(prim.markers = c("Satb2", "Bcl11b", "Cux1", "Fezf2", "Tbr1"),
                                            markers.list = c(),
                                            dir, condition,
                                            p_value = 0.001,
                                            symmetric = TRUE) {
  standardGeneric("plot_general.heatmap")
})
#' @rdname plot_general.heatmap
setMethod(
  "plot_general.heatmap", "ANY",
  function(prim.markers = c("Satb2", "Bcl11b", "Cux1", "Fezf2", "Tbr1"), markers.list = c(), dir,
           condition, p_value = 0.001, symmetric = TRUE) {
    print("ploting a general heatmap")
    ET <- NULL

    if (symmetric == TRUE) {
      markers.list <- as.list(c(unlist(prim.markers), unlist(markers.list)))
    }

    if (is.null(markers.list)) {
      markers.list <- as.list(prim.markers)
    } else {
      markers.list <- as.list(markers.list)
    }

    obj <- readRDS(paste0(dir, condition, ".cotan.RDS"))

    no_genes <- unique(c(unlist(markers.list), prim.markers))[!unique(c(
      unlist(markers.list),
      prim.markers
    ))
    %in% rownames(obj@coex)]

    if (!rlang::is_empty(no_genes)) {
      print(paste0(no_genes, " not present!"))
    }

    pval <- calculatePValue(object = obj)
    diag(pval) <- 1
    pval <- pval[, unique(c(unlist(markers.list), prim.markers))]

    pval.red <- apply(pval, 1, FUN = min)
    genes.row <- names(pval.red[pval.red < p_value])

    genes.row <- unique(c(unique(c(unlist(markers.list), prim.markers)), genes.row))
    pval <- as.data.frame(pval)

    coex <- getGenesCoex(obj)
    diag(coex) <- 0
    if (symmetric == TRUE) {
      coex <- coex[rownames(coex) %in% genes.row, colnames(coex) %in% genes.row]
    } else {
      coex <- coex[rownames(coex) %in% genes.row, ]
    }

    list.rows <- c()
    for (m in unlist(markers.list)) {
      genes <- rownames(pval[pval[, m] < p_value, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0, ])]
      list.rows[[m]] <- genes
    }

    list.cols <- c()
    for (m in prim.markers) {
      genes <- rownames(pval[pval[, m] < p_value, ])
      genes <- genes[genes %in% rownames(coex[coex[, m] > 0, ])]
      list.cols[[m]] <- genes
    }

    cl.genes.rows <- c()
    for (ll in names(list.rows)) {
      tmp <- data.frame("genes" = list.rows[[ll]], "cl" = rep(ll, length(list.rows[[ll]])))
      cl.genes.rows <- rbind(cl.genes.rows, tmp)
    }

    cl.genes.rows <- cl.genes.rows[cl.genes.rows$genes %in% rownames(coex), ]

    reorder_idx_row <- match(cl.genes.rows$gene, rownames(coex))


    if (symmetric == TRUE) {
      cl.genes.cols <- data.frame()
      for (ll in names(list.rows)) {
        tmp <- data.frame("genes" = list.rows[[ll]], "cl" = rep(ll, length(list.rows[[ll]])))
        cl.genes.cols <- rbind(cl.genes.cols, tmp)
      }
    } else {
      cl.genes.cols <- data.frame()
      for (ll in names(list.cols)) {
        tmp <- data.frame("genes" = list.cols[[ll]], "cl" = rep(ll, length(list.cols[[ll]])))
        cl.genes.cols <- rbind(cl.genes.cols, tmp)
      }
    }
    cl.genes.cols <- cl.genes.cols[cl.genes.cols$genes %in% colnames(coex), ]

    reorder_idx_col <- match(cl.genes.cols$gene, colnames(coex))


    to.plot <- coex[reorder_idx_row, reorder_idx_col]

    col_fun <- circlize::colorRamp2(
      c(
        round(stats::quantile(as.matrix(to.plot), probs = 0.001),
          digits = 3
        ), 0,
        round(stats::quantile(as.matrix(to.plot), probs = 0.999),
          digits = 3
        )
      ),
      c("#E64B35FF", "gray93", "#3C5488FF")
    )

    # The next line is to set the columns and raws order
    # need to be implemented
    part1 <- ComplexHeatmap::Heatmap(as.matrix(to.plot),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_split = cl.genes.rows$cl,
      column_split = cl.genes.cols$cl,
      col = col_fun,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title_gp = grid::gpar(
        fill = "#8491B44C", font = 3,
        col = "#3C5488FF"
      ),
      row_title_gp = grid::gpar(fill = "#8491B44C", font = 3, col = "#3C5488FF")
    )
    lgd <- ComplexHeatmap::Legend(
      col_fun = col_fun, title = "coex", grid_width =
        unit(0.3, "cm"),
      direction = "horizontal", title_position = "topcenter",
      title_gp = grid::gpar(fontsize = 10, fontface = "bold", col = "#3C5488FF"),
      labels_gp = grid::gpar(col = "#3C5488FF", font = 3)
    )
    ComplexHeatmap::draw(part1,
      show_heatmap_legend = FALSE,
      annotation_legend_list = lgd, annotation_legend_side = "bottom"
    )
  }
)


