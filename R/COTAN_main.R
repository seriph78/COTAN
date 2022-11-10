## Main functions


#' initRaw
#'
#' It starts to fill some fields of cotan object.
#' @param object the dataframe containing the raw data: it should be a data.frame, not a matrix.
#' @param GEO a code reporting the GEO identification or other specific dataset
#' code
#' @param sc.method a string reporting the method used for the sequencing
#' @param cond a string reporting the specific sample condition or time point
#'
#' @return A COTAN object to store all information
#' @import methods
#' @export
#' @rdname initRaw
#' @examples
#'
#' data("raw.dataset")
#' obj <- new("scCOTAN", raw = raw.dataset)
#' obj <- initRaw(obj, GEO = "code", sc.method = "10X", cond = "mouse dataset")
#'
setGeneric("initRaw", function(object, GEO, sc.method = "10X", cond) {
  standardGeneric("initRaw")
})
#' @rdname initRaw
setMethod("initRaw","scCOTAN",
            function(object,GEO,sc.method,cond) {
                print("Initializing S4 object")
                if( (!attr(object@raw, "class")[1] == "dgCMatrix") |
                   is.null(attr(object@raw, "class")) ){
                    object@raw  <- methods::as(as.matrix(object@raw),
                                             "sparseMatrix")
                }
                if(! all((object@raw - round(object@raw)) == 0)){
                    print("WARNING! Input data contains not integer numbers!")

                }
                object@meta[1,seq_len(2)] = c("GEO:",GEO)
                object@meta[2,seq_len(2)] = c("scRNAseq method:",sc.method)
                object@meta[3,1] = "starting n. of cells:"
                object@meta[3,2] = ncol(object@raw)
                object@meta[4,seq_len(2)] = c("Condition sample:",cond)

              object@clusters = rep(NA,ncol(object@raw))
              names(object@clusters)=getCells(object)
                return(object)
          }
)


#' clean
#'
#' Main function that can be used to check and clean the dataset. It also
#' produce
#' (using the function fun_linear) and store
#' the estimators for nu and lambda. It also fill the raw.norm (raw / nu) and
#' n_cell
#' (the initial number of cells in the dataset)
#' @param object COTAN object
#' @return a list of objects containing: "cl1" is the first cell cluster, "cl2"
#' is the
#' second cell cluster,
#' "pca.cell.2" is a ggplot2 cell pca plot, "object" is the COTAN object with
#' saved the
#' estimated lambda and mu, "mu_estimator", "D"
#' "pca_cells" pca numeric data.
#' @export
#'
#' @import ggplot2
#' dplyr
#'
#' @importFrom tibble rownames_to_column
#' @importFrom stats hclust
#' @importFrom stats cutree
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix t
#' @importFrom utils head
#' @importFrom Matrix colMeans
#' @rdname clean
#' @examples
#' data("ERCC.cotan")
#' ttm <- clean(ERCC.cotan)
#'
setGeneric("clean", function(object) standardGeneric("clean"))
#' @rdname clean
setMethod("clean","scCOTAN",

              function(object) {

                 mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")

                 my_theme  <- theme(axis.text.x = element_text(size = 14,
                                                             angle = 0, hjust = .5,
                                                             vjust = .5,
                                                              face = "plain",
                                                             colour ="#3C5488FF" ),
                                   axis.text.y = element_text( size = 14, angle = 0, hjust = 0,
                                                               vjust = .5,
                                                               face = "plain",
                                                               colour ="#3C5488FF"),
                                   axis.title.x = element_text( size = 14,
                                                                angle = 0, hjust = .5,
                                                                vjust = 0,
                                                                face = "plain",
                                                                colour ="#3C5488FF"),
                                   axis.title.y = element_text( size = 14, angle = 90, hjust = .5,
                                                                vjust = .5,
                                                                face = "plain",
                                                                colour ="#3C5488FF"))
                 print("Starting")
                 pc1 <- PC2 <- PC1 <- NULL


                cells  <-  get.zero_one.cells(object)
                
                object@raw  <-  object@raw[rownames(cells), colnames(cells)]
                print("Cells/genes selection done")
                rm(cells)
                gc()

                list1  <- fun_linear(object)
                print("Fun linear DONE")


                dist_cells  <- list1$dist_cells
                pca_cells  <-  list1$pca_cells
                
                #t_to_clust  <-  as.matrix(list1$t_to_clust)
                #mu_estimator  <-  list1$mu_estimator

                object  <-  list1$object
                to_clust <-  list1$to_clust
                rm(list1)
                gc()



                raw_norm <- Matrix::t(Matrix::t(object@raw) *
                                  (1/(as.vector(object@nu))))
                object@raw.norm <- raw_norm
                rm(raw_norm)
                gc()
                
                print("starting hclust")

                hc_cells <- hclust(dist_cells, method = "complete")
                gc()
                
                groups <- cutree(hc_cells, k=2)

                if(length(groups[groups == 1]) < length(groups[groups == 2])  ){
                  groups[groups == 1] <- "B"
                  groups[groups == 2] <- "A"

                }else{
                  groups[groups == 1] <- "A"
                  groups[groups == 2] <- "B"
                }

                cl2 <- names(which(cutree(hc_cells, k = 2) == 2))
                cl1 <- names(which(cutree(hc_cells, k = 2) == 1))

                if (length(cl2) > length(cl1) ) {
                  cl2  <-  names(which(cutree(hc_cells, k = 2) == 1))
                  cl1  <-  names(which(cutree(hc_cells, k = 2) == 2))
                }
                
                t_to_clust <- Matrix::t(to_clust)
                t_to_clust <- round(t_to_clust,digits = 4)
                t_to_clust <- as.data.frame(as.matrix(t_to_clust))

                if (identical(rownames(t_to_clust),names(groups))) {
                  t_to_clust <- cbind(t_to_clust,groups)
                  
                }else{
                  stop("Error in the cell names")
                }
                
              # ---- next: to check which genes are specific for the B group of cells

                to_clust <- as.matrix(round(to_clust,digits = 4))
                
                B <- as.data.frame(to_clust[,colnames(to_clust) %in% cl2])
                rm(to_clust)
                gc()
                
                colnames(B)<-cl2
                B <- rownames_to_column(B)
                if (dim(B)[2]>2) {
                  B <- B[order(rowMeans(B[,2:length(colnames(B))]),
                              decreasing = TRUE), ]
                }else{
                  B <- B[order(B[,2],decreasing = TRUE), ]
                }

                #print(utils::head(B, 15))

                C <-  arrange(B,rowMeans(B[2:length(colnames(B))]))
                rownames(C) <-  C$rowname
                D <-  data.frame("means" = rowMeans(C[2:length(colnames(C))]),
                               "n" = NA )
                D <-  D[D$means>0,]
                D$n <-  seq_along(D$means)

                #check if the pca plot is clean enought and from the printed genes,
                #if the smalest group of cells are caratterised by particular genes

                pca_cells  <-  cbind(pca_cells,"groups"=t_to_clust$groups)

                pca.cell.1  <-  ggplot(subset(pca_cells,groups == "A" ),
                                    aes(x=PC1, y=PC2,colour =groups)) +
                  geom_point(alpha = 0.5, size=3)

                pca.cell.2  <- pca.cell.1 + geom_point(data = subset(pca_cells,
                            groups != "A" ),aes(x=PC1, y=PC2,colour =groups),
                                                   alpha = 0.8, size=3)+
                  scale_color_manual("groups", values = mycolours)  +
                  my_theme + theme(legend.title = element_blank(),
                                   legend.text = element_text( size = 12,
                                                            color = "#3C5488FF",
                                                            face ="italic" ),
                                   legend.position="bottom")

                #genes plot 
                pl <- ggplot(D, aes(x=n,y=means)) + geom_point() +
                  geom_text_repel(data=subset(D, n > (max(D$n) - 15) ), aes(n,means,label=rownames(D[D$n > (max(D$n)- 15),])),
                                  nudge_y      = 0.05,
                                  nudge_x      = 0.05,
                                  direction    = "x",
                                  angle        = 90,
                                  vjust        = 0,
                                  segment.size = 0.2)+
                  ggtitle(label = "B cell group genes mean expression")+my_theme +
                  theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 1,hjust = 0.02 ),
                        plot.subtitle = element_text(color = "darkred",vjust = - 15,hjust = 0.01 ))
                
              
                # UDE/nu plot
                nu_est = round(getNu(object), digits = 7)
                
                plot.nu <-ggplot(pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))
                
                plot.nu = plot.nu + geom_point(size = 1,alpha= 0.8)+
                  scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                                        midpoint = log(mean(nu_est)),name = "ln(nu)")+
                  ggtitle("Cells PCA coloured by cells efficiency") +
                  my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                                    legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                                    legend.text = element_text(color = "#3C5488FF", size = 11),
                                    legend.key.width = unit(2, "mm"),
                                    legend.position="right")
                

                output  <- list("cl1"=cl1,
                                "cl2"=cl2,
                                "pca.cell.2"=pca.cell.2,
                                "object"=object,
                                #"mu_estimator"=mu_estimator,
                                "D"=D,
                                "pca_cells"=pca_cells,
                                "genes.plot"=pl,
                                "UDE.plot"= plot.nu)
                return(output)
          }
)


#' cotan_analysis
#'
#' This is the main function that estimates the a vector to store all the
#' negative binomial
#' dispersion factors. It need to be run after \code{\link{clean}}
#' @param object A COTAN object
#' @param cores number of cores to use. Default is 11.
#' @importFrom parallel mclapply
#' @return a COTAN object
#' @export
#' @rdname cotan_analysis
#' @examples
#' data("ERCC.cotan")
#' ERCC.cotan <- cotan_analysis(ERCC.cotan)
#'
setGeneric("cotan_analysis", function(object, cores = 1) {
  standardGeneric("cotan_analysis")
})
#' @rdname cotan_analysis
setMethod("cotan_analysis","scCOTAN",
          function(object, cores= 1) {

              print("cotan analysis")
              if (Sys.info()['sysname'] == "Windows") {
                  print("On windows the numebr of cores used will be 1! Multicore is not supported.")
                  cores= 1
              }

              cells  <-  object@raw
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              #cells[cells <= 0] <- 0

              # exclude the effectively ubiquitous genes and saved in a separate file
              mu_estimator <- estimateMu(object)

              object <- as(housekeepingGenes(object), "scCOTAN")
              
              hk <- object@hk

              mu_estimator <- mu_estimator[!rownames(mu_estimator) %in% hk,]
              cells <- cells[!rownames(cells) %in% hk, ]

              mu_estimator <- as.matrix(mu_estimator)
              print("start a minimization")
              gc()
              tot <- list()

              p_begin <- 1
              n_genes = nrow(mu_estimator)
              while(p_begin <= n_genes) {
                  p_end <- p_begin + 200
                  if(p_end >= n_genes){
                      print("Final trance!")
                      p_end = n_genes
                  }
                  tot1 <- parallel::mclapply(
                                  rownames(mu_estimator)[p_begin:p_end],
                                  dispersionBisection, zeroOneMatrix = cells,
                                  muEstimator = mu_estimator, mc.cores = cores)
                  tot <- append(tot, tot1)
                  p_begin <- p_end+1
                  if((p_begin %% 10)==0){
                      print(paste("Next gene:", rownames(mu_estimator)[p_begin],
                                  "number", p_begin,
                                  sep = " "))
                  }
              }
              gc()

              tot2 <- tot[[1]]
              for (tt in 2:length(tot)) {
                  tot2<- rbind(tot2,tot[[tt]])
              }

              object@a <- tot2$dispersion
              names(object@a) <- rownames(tot2)
              print(paste("a min:",  min(tot2$dispersion),
                          "| a max", max(tot2$dispersion),
                          "| negative a %:", sum(tot2$dispersion <0)/nrow(tot2)*100,
                          sep=" "))
              
              gc()

              return(object)
          }
)






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
#' @importFrom Matrix forceSymmetric
#' @export
#'
#' @import ggplot2
#' @import tidyr
#' @import scales
#' @rdname plot_heatmap
#' @examples
#' data("ERCC.cotan")
#' data_dir <- tempdir()
#' saveRDS(ERCC.cotan, file = file.path(data_dir, "ERCC.cotan.RDS"))
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
setGeneric("plot_heatmap", function(p_val.tr = 0.05, df_genes, sets, conditions, dir) {
  standardGeneric("plot_heatmap")
})
#' @rdname plot_heatmap
setMethod(
  "plot_heatmap", "ANY",
  function(p_val.tr = 0.05, df_genes, sets, conditions, dir) {
    time <- g2 <- NULL
    print("plot heatmap")
    gr <- df_genes[[1]]
    ge <- unique(array(sort(unlist(df_genes[sets]))))
    df.to.print <- data.frame()
    for (ET in conditions) {
      print(paste("Loading condition", ET, sep = " "))
      obj <- readRDS(paste(dir, ET, ".cotan.RDS", sep = ""))
      if (is(class(obj@coex)[1], "dtCMatrix")) {
        print("COTAN object in the old format! Converting...")
        obj <- get.coex(obj)
        print(paste("Saving as new file as ", dir, ET, "new.cotan.RDS", sep = ""))
        saveRDS(obj, paste(dir, ET, "new.cotan.RDS", sep = ""))
      }
      if (any(gr %in% obj@coex$genes) == FALSE) {
        paste0("primary markers all absent in ", ET)
        stop()
      }
      p_val <- get.pval(obj, gene.set.col = gr, gene.set.row = ge)
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

      coex <- vec2mat_rfast(obj@coex, genes = gr)
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
    print(paste("min coex:", min(df.to.print$coex, na.rm = TRUE), "max coex", max(df.to.print$coex, na.rm = TRUE), sep = " "))
    heatmap <- ggplot(data = subset(df.to.print, type %in% names(df_genes)[sets]), aes(time, factor(g2, levels = rev(levels(factor(g2)))))) +
      geom_tile(aes(fill = coex), colour = "black", show.legend = TRUE) +
      facet_grid(type ~ g1, scales = "free", space = "free") +
      scale_fill_gradient2(
        low = "#E64B35FF", mid = "gray93", high = "#3C5488FF",
        midpoint = 0,
        na.value = "grey80", space = "Lab", guide = "colourbar",
        aesthetics = "fill", oob = scales::squish
      ) +
      theme(
        axis.title.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect(fill = "#8491B44C"),
        strip.text.y = element_text(size = 9, colour = "#3C5488FF"),
        strip.text.x = element_text(size = 9, angle = 90, colour = "#3C5488FF"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(
          size = 9, angle = 0, hjust = 0, vjust = .5,
          face = "plain", colour = "#3C5488FF"
        ),
        legend.text = element_text(color = "#3C5488FF", face = "italic"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.height = unit(2, "mm")
      )
    heatmap
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
#' data("ERCC.cotan")
#' data_dir <- tempdir()
#' saveRDS(ERCC.cotan, file = file.path(data_dir, "ERCC.cotan.RDS"))
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

    obj <- readRDS(paste(dir, condition, ".cotan.RDS", sep = ""))

    if (is(class(obj@coex)[1], "dtCMatrix")) {
      print("COTAN object in the old format! Converting...")
      obj <- get.coex(obj)
      print(paste("Saving as new file as ", dir, ET, "new.cotan.RDS", sep = ""))
      saveRDS(obj, paste(dir, ET, "new.cotan.RDS", sep = ""))
    }

    no_genes <- unique(c(unlist(markers.list), prim.markers))[!unique(c(
      unlist(markers.list),
      prim.markers
    ))
    %in% obj@coex$genes]

    if (!rlang::is_empty(no_genes)) {
      print(paste(no_genes, "not present!", sep = " "))
    }

    pval <- get.pval(object = obj)
    diag(pval) <- 1
    pval <- pval[, unique(c(unlist(markers.list), prim.markers))]

    pval.red <- apply(pval, 1, FUN = min)
    genes.row <- names(pval.red[pval.red < p_value])

    genes.row <- unique(c(unique(c(unlist(markers.list), prim.markers)), genes.row))
    pval <- as.data.frame(pval)

    coex <- vec2mat_rfast(obj@coex)
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


#' plot_GDI
#'
#' This function directly evaluate and plot the GDI for a sample.
#'
#' @param object A COTAN object
#' @param cond A string corresponding to the condition/sample (it is used only for the title)
#' @param type Type of statistic to be used. Default is "S":
#' Pearson's chi-squared test statistics. "G" is G-test statistics
#'
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#' @importFrom  stats quantile
#' @importFrom Matrix forceSymmetric
#' @rdname plot_GDI
#' @examples
#' data("ERCC.cotan")
#' plot_GDI(ERCC.cotan, cond = "ERCC")
setGeneric("plot_GDI", function(object, cond, type = "S") standardGeneric("plot_GDI"))
#' @rdname plot_GDI
setMethod(
  "plot_GDI", "scCOTAN",
  function(object, cond, type = "S") {
    ET <- sum.raw.norm <- NULL

    if (is(class(object@coex)[1], "dtCMatrix")) {
      print("COTAN object in the old format! Converting...")
      object <- get.coex(object)
      print(paste("Saving as new file as ", dir, ET, "new.cotan.RDS", sep = ""))
      saveRDS(object, paste(dir, ET, "new.cotan.RDS", sep = ""))
    }

    print("GDI plot ")
    if (type == "S") {
      GDI <- get.GDI(object, type = "S")
    } else if (type == "G") {
      print("Using G")
      GDI <- get.GDI(object, type = "G")
    }

    si <- 12
    plot <- ggplot(GDI, aes(x = sum.raw.norm, y = GDI)) +
      geom_point(size = 2, alpha = 0.5, color = "#8491B4B2") +
      geom_hline(yintercept = 1.5, linetype = "dotted", color = "darkred", size = 1) +
      geom_hline(
        yintercept = stats::quantile(GDI$GDI)[4], linetype = "dashed",
        color = "darkblue"
      ) +
      geom_hline(
        yintercept = stats::quantile(GDI$GDI)[3], linetype = "dashed",
        color = "darkblue"
      ) +
      xlab("log normalized reads sum") +
      ylab("global p val index (GDI)") +
      ggtitle(paste("GDI ", cond, sep = " ")) +
      theme(
        axis.text.x = element_text(
          size = si, angle = 0, hjust = .5, vjust = .5,
          face = "plain", colour = "#3C5488FF"
        ),
        axis.text.y = element_text(
          size = si, angle = 0, hjust = 0, vjust = .5,
          face = "plain", colour = "#3C5488FF"
        ),
        axis.title.x = element_text(
          size = si, angle = 0, hjust = .5, vjust = 0,
          face = "plain", colour = "#3C5488FF"
        ),
        axis.title.y = element_text(
          size = si, angle = 90, hjust = .5, vjust = .5,
          face = "plain", colour = "#3C5488FF"
        ),
        legend.title = element_blank(),
        plot.title = element_text(color = "#3C5488FF", size = 14, face = "bold.italic"),
        legend.text = element_text(color = "#3C5488FF", face = "italic"),
        legend.position = "bottom"
      )

    return(plot)
  }
)


#' get.observed.ct
#'
#' This function is used to get the observed contingency table for a given pair of genes.
#'
#' @param object The cotan object
#' @param g1 A gene
#' @param g2 The other gene
#'
#' @return A contingency table as dataframe
#' @export
#' @rdname get.observed.ct
#' @examples
#' data("ERCC.cotan")
#' g1 <- getGenes(ERCC.cotan)[sample(getNumGenes(ERCC.cotan), 1)]
#' g2 <- getGenes(ERCC.cotan)[sample(getNumGenes(ERCC.cotan), 1)]
#' get.observed.ct(object = ERCC.cotan, g1 = g1, g2 = g2)
setGeneric("get.observed.ct", function(object, g1, g2) standardGeneric("get.observed.ct"))
#' @rdname get.observed.ct
setMethod(
  "get.observed.ct", "scCOTAN",
  function(object, g1, g2) {
    yes_yes <- observedContingencyYY(object)
    si_si <- yes_yes[g1, g2]

    si_any <- max(yes_yes[g1, ])

    any_si <- max(yes_yes[g2, ])
    si_no <- si_any - si_si
    no_si <- any_si - si_si
    no_no <- getNumCells(object) - (si_si + si_no + no_si)

    ct <- as.data.frame(matrix(ncol = 2, nrow = 2))
    colnames(ct) <- c(paste(g1, "yes", sep = "."), paste(g1, "no", sep = "."))
    rownames(ct) <- c(paste(g2, "yes", sep = "."), paste(g2, "no", sep = "."))
    ct[1, 1] <- si_si
    ct[1, 2] <- no_si
    ct[2, 2] <- no_no
    ct[2, 1] <- si_no

    return(ct)
  }
)


#' get.expected.ct
#'
#' This function is used to get the expected contingency table for a given pair of genes.
#' @param object The cotan object
#' @param g1 A gene
#' @param g2 The other gene
#'
#' @return A contingency table as dataframe
#' @export
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @rdname get.expected.ct
#' @examples
#' data("ERCC.cotan")
#' g1 <- getGenes(ERCC.cotan)[sample(getNumGenes(ERCC.cotan), 1)]
#' g2 <- getGenes(ERCC.cotan)[sample(getNumGenes(ERCC.cotan), 1)]
#' while (g1 %in% getHousekeepingGenes(ERCC.cotan)) {
#'   g1 <- getGenes(ERCC.cotan)[sample(getNumGenes(ERCC.cotan), 1)]
#' }
#'
#' while (g2 %in% getHousekeepingGenes(ERCC.cotan)) {
#'   g2 <- getGenes(ERCC.cotan)[sample(getNumGenes(ERCC.cotan), 1)]
#' }
#' get.expected.ct(object = ERCC.cotan, g1 = g1, g2 = g2)
setGeneric("get.expected.ct", function(object, g1, g2) standardGeneric("get.expected.ct"))
#' @rdname get.expected.ct
setMethod(
  "get.expected.ct", "scCOTAN",
  function(object, g1, g2) {
    stopifnot("a gene is constitutive!" = !(g1 %in% object@hk | g2 %in% object@hk))

    mu_estimator <- object@lambda[c(g1, g2)] %*% t(object@nu)

    cells <- as.matrix(object@raw)[c(g1, g2), ]
    #---------------------------------------------------

    cells[cells > 0] <- 1
    cells[cells <= 0] <- 0

    M <- funProbZero(object@a[c(g1, g2)], mu_estimator)
    rownames(M) <- c(g1, g2)
    N <- 1 - M

    n_zero_esti <- rowSums(M)
    n_zero_obs <- rowSums(cells[!rownames(cells) %in% object@hk, ] == 0)

    dist_zeros <- sqrt(sum((n_zero_esti - n_zero_obs)^2))
    stopifnot("Errore: some Na in matrix M " = !any(is.na(M)))

    gc()
    estimator_no_no <- M %*% t(M)
    estimator_no_si <- M %*% t(N)
    estimator_si_no <- t(estimator_no_si)
    estimator_si_si <- N %*% t(N)


    ct <- as.data.frame(matrix(ncol = 2, nrow = 2))
    colnames(ct) <- c(paste(g1, "yes", sep = "."), paste(g1, "no", sep = "."))
    rownames(ct) <- c(paste(g2, "yes", sep = "."), paste(g2, "no", sep = "."))
    ct[1, 1] <- estimator_si_si[g1, g2]
    ct[1, 2] <- estimator_no_si[g1, g2]
    ct[2, 2] <- estimator_no_no[g1, g2]
    ct[2, 1] <- estimator_si_no[g1, g2]

    return(ct)
  }
)


#' drop.genes.cells
#'
#' This function remove an array of genes and/or cells from the original object raw matrix.
#'
#' @param object a COTAN object
#' @param genes an array of gene names
#' @param cells an array of cell names
#'
#' @return the original object but with the raw matrix without the indicated cells and/or genes.
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' genes.to.rem <- getGenes(ERCC.cotan)[grep('^MT', getGenes(ERCC.cotan))]
#' cells.to.rem <- getCells(ERCC.cotan)[which(getCellsSize(ERCC.cotan) == 0)]
#' ERCC.cotan <- drop.genes.cells(ERCC.cotan, genes.to.rem, cells.to.rem)
setGeneric("drop.genes.cells", function(object, genes = c(), cells = c()) standardGeneric("drop.genes.cells"))
#' @rdname drop.genes.cells
setMethod(
  "drop.genes.cells", "scCOTAN",
  function(object, genes, cells) {
    object@raw <- object@raw[!rownames(object@raw) %in% genes, ]
    object@raw <- object@raw[, !colnames(object@raw) %in% cells]
    return(object)
  }
)

#' add.row.to.meta
#'
#' This function is used to add a line of information to the information data frame (metadata).
#'
#' @param object a COTAN object
#' @param text.line an array containing the information
#'
#' @return the updated COTAN object
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' text <- c("Test", "This is a test")
#' ERCC.cotan <- add.row.to.meta(ERCC.cotan, text)
#' getMetadataDataset(ERCC.cotan)
setGeneric("add.row.to.meta", function(object, text.line) standardGeneric("add.row.to.meta"))
#' @rdname add.row.to.meta
setMethod(
  "add.row.to.meta", "scCOTAN",
  function(object, text.line) {
    object@meta[(nrow(object@meta) + 1), seq_len(length(text.line))] <- text.line
    return(object)
  }
)
