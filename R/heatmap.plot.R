
#' heatmapPlot
#'
#' @description This function creates the `heatmap` of one or more `COTAN`
#'   objects.
#'
#' @param genesLists A `list` of genes' `array`s. The first `array` defines the
#'   genes in the columns.
#' @param sets A numeric array indicating which fields in the
#'   previous `list` should be used.
#' @param conditions An array of prefixes indicating the different files.
#' @param pValueThreshold The p-value threshold. Default is 0.05.
#' @param dir The directory in which are all `COTAN` files (corresponding to the
#'   previous prefixes)
#'
#' @returns a `ggplot` object
#'
#' @importFrom Matrix forceSymmetric
#'
#' @importFrom tidyr pivot_longer
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_fill_gradient2
#'
#' @importFrom scales squish
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
#' heatmapPlot(
#'   pValueThreshold = 0.05, genesLists = gene.sets.list,
#'   sets = c(2, 3), conditions = c("ERCC"), dir = paste0(data_dir, "/")
#' )
#'
#' @rdname heatmapPlot
#'
heatmapPlot <- function(genesLists, sets, conditions,
                        pValueThreshold = 0.05, dir) {
  #time <- g2 <- NULL
  logThis("heatmap plot: START", logLevel = 2)

  colGenes <- genesLists[[1]]
  allGenes <- unique(sort(unlist(genesLists[sets])))

  df.to.print <- data.frame()
  for (cond in conditions) {
    logThis(paste0("Loading condition '", cond, "'"), logLevel = 3)
    obj <- readRDS(file.path(dir, paste0(cond, ".cotan.RDS")))

    stopifnot("primary markers all absent" <- any(colGenes %in% getGenes(obj)))

    pValue <- calculatePValue(obj, geneSubsetCol = colGenes, geneSubsetRow = allGenes)
    pValue <- as.data.frame(pValue)

    pValue$g2 <- as.vector(rownames(pValue))
    df.temp.pval <- pivot_longer(pValue, cols = seq_along(colnames(pValue)) - 1,
                                 names_to = "g1", values_to = "pValue")

    #---------------------------------------------------------
    coex <- getGenesCoex(obj)
    diag(coex) <- 0
    coex <- coex[getGenes(obj) %in% allGenes, getGenes(obj) %in% colGenes]
    coex <- as.data.frame(coex)

    coex$g2 <- as.vector(rownames(coex))
    df.temp.coex <- pivot_longer(coex, cols = seq_along(colnames(pValue)) - 1,
                                 names_to = "g1", values_to = "coex")

    df.temp <- merge(df.temp.coex, df.temp.pval)
    df.temp$time <- cond
    df.temp$type <- NA
    df.temp$absent <- NA
    df.temp2 <- data.frame()
    for (type in names(genesLists)[sets]) {
      for (g1 in colGenes) {
        tt <- df.temp[df.temp$g2 %in% genesLists[[type]] & df.temp$g1 == g1, ]
        # control if the subset is smaller than the number of wanted genes
        if (dim(tt)[1] < length(genesLists[[type]])) {
          n.row <- length(genesLists[[type]]) - dim(tt)[1]
          t.rows <- as.data.frame(matrix(nrow = n.row, ncol = 7))
          colnames(t.rows) <- colnames(tt)
          t.rows[, "g1"] <- g1
          t.rows[, "time"] <- cond
          t.rows[, "absent"] <- "yes"
          t.rows[, "pValue"] <- 1
          t.rows[, "g2"] <- genesLists[[type]][!genesLists[[type]] %in% tt$g2]
          tt <- rbind(tt, t.rows)
        }
        tt$type <- type
        df.temp2 <- rbind(df.temp2, tt)
      }
      print(type)
    }
    df.temp <- df.temp2
    df.temp$t_hk <- ifelse((df.temp$g2 %in% obj@hk) | (df.temp$g1 %in% obj@hk), "hk", "n")
    df.temp[df.temp$pValue > pValueThreshold, ]$coex <- 0
    df.to.print <- rbind(df.to.print, df.temp)
  }
  print(paste("min coex:", min(df.to.print$coex, na.rm = TRUE),
              "max coex",  max(df.to.print$coex, na.rm = TRUE)))

  heatmap <- ggplot(data = subset(df.to.print, type %in% names(genesLists)[sets]),
                    aes(time, factor(g2, levels = rev(levels(factor(g2)))))) +
             geom_tile(aes(fill = coex), colour = "black", show.legend = TRUE) +
             facet_grid(type ~ g1, scales = "free", space = "free") +
             scale_fill_gradient2(low = "#E64B35FF", mid = "gray93", high = "#3C5488FF",
                                  midpoint = 0,na.value = "grey80", space = "Lab",
                                  guide = "colourbar", aesthetics = "fill", oob = squish) +
             plotTheme("heatmap", textSize = 9)

  return(heatmap)
}



#' plot_general.heatmap
#'
#' This function is used to plot an heatmap made using only some genes, as markers,
#' and collecting all other genes correlated
#' with these markers with a p-value smaller than the set threshold. Than all relations are plotted.
#' Primary markers will be plotted as groups of rows. Markers list will be plotted as columns.
#' @param prim.markers A set of genes plotted as rows.
#' @param markers.list A set of genes plotted as columns.
#' @param dir The directory where the `COTAN` object is stored.
#' @param condition The prefix for the `COTAN` object file.
#' @param p_value The p-value threshold
#' @param symmetric A boolean: default F. If T the union of prim.markers and marker.list
#' is sets as both rows and column genes
#'
#' @return A ggplot2 object
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap Legend
#' @importFrom ComplexHeatmap draw
#' @importFrom grid gpar
#' @importFrom stats quantile
#' @importFrom circlize colorRamp2
#' @importFrom Matrix forceSymmetric
#' @importFrom rlang is_empty
#'
#' @rdname plot_general.heatmap
#'
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
setGeneric("plot_general.heatmap", function(prim.markers,
                                            markers.list = c(),
                                            dir, condition,
                                            p_value = 0.001,
                                            symmetric = TRUE) {
  standardGeneric("plot_general.heatmap")
})
#' @rdname plot_general.heatmap
setMethod(
  "plot_general.heatmap", "ANY",
  function(prim.markers, markers.list = c(), dir,
           condition, p_value = 0.001, symmetric = TRUE) {
    print("ploting a general heatmap")

    if (symmetric == TRUE) {
      markers.list <- as.list(c(unlist(prim.markers), unlist(markers.list)))
    }

    if (is.null(markers.list)) {
      markers.list <- as.list(prim.markers)
    } else {
      markers.list <- as.list(markers.list)
    }

    obj <- readRDS(file.path(dir, paste0(condition, ".cotan.RDS")))
    obj <- as(obj, "scCOTAN")

    no_genes <- unique(c(unlist(markers.list), prim.markers))[!unique(c(
      unlist(markers.list),
      prim.markers
    ))
    %in% rownames(obj@coex)]

    if (!is_empty(no_genes)) {
      print(paste0(no_genes, " not present!"))
    }

    pval <- calculatePValue(obj)
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

    col_fun <- colorRamp2(
      c(
        round(quantile(as.matrix(to.plot), probs = 0.001),
          digits = 3
        ), 0,
        round(quantile(as.matrix(to.plot), probs = 0.999),
          digits = 3
        )
      ),
      c("#E64B35FF", "gray93", "#3C5488FF")
    )

    # The next line is to set the columns and raws order
    # need to be implemented
    part1 <- Heatmap(as.matrix(to.plot),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_split = cl.genes.rows$cl,
      column_split = cl.genes.cols$cl,
      col = col_fun,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title_gp = gpar(
        fill = "#8491B44C", font = 3,
        col = "#3C5488FF"
      ),
      row_title_gp = gpar(fill = "#8491B44C", font = 3, col = "#3C5488FF")
    )
    lgd <- Legend(
      col_fun = col_fun, title = "coex", grid_width =
        unit(0.3, "cm"),
      direction = "horizontal", title_position = "topcenter",
      title_gp = gpar(fontsize = 10, fontface = "bold", col = "#3C5488FF"),
      labels_gp = gpar(col = "#3C5488FF", font = 3)
    )
    draw(part1,
      show_heatmap_legend = FALSE,
      annotation_legend_list = lgd, annotation_legend_side = "bottom"
    )
  }
)


