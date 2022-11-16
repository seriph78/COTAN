#' cluster_homogeneity
#'
#' # Function that takes a seurat object and the whole raw data (after cleaning).
#' It runs cotan on each cluster and check if the GDI is lower than 1.5 for the 99% of genes.
#' If it is too high,
#' the cluster is not uniform and so it save the cells in an array to cluster again.
#' @param data.seurat a Seurat object with clusters already estimated
#' @param data.raw raw data matrix
#' @param out_dir path to the output directory
#' @param cond string defining the condition
#' @param cores number of cores used
#'
#' @import Seurat
#' @import ggrepel
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @return an array of cells that need to be re-clustered
#' @export
#'
#' @examples
setGeneric("cluster_homogeneity", function(data.seurat,data.raw,out_dir, cond, cores=1) standardGeneric("cluster_homogeneity"))
#' @rdname cluster_homogeneity
setMethod("cluster_homogeneity","scCOTAN",
          function(data.seurat,data.raw,out_dir, cond, cores){

            options(ggrepel.max.overlaps = Inf)

            # This detect the number of clusters
            tot.clusters = max(as.numeric(as.vector(unique(data.seurat@meta.data$seurat_clusters[!data.seurat@meta.data$seurat_clusters == "singleton" ]))))
            print(paste("Number of clusters:", tot.clusters))
            to_rec = NA
            for (cl in c(0:tot.clusters)) { #0

              print(paste("Working on cluster", cl, "out of", tot.clusters))
              raw = as.data.frame(data.raw[,rownames(data.seurat@meta.data[data.seurat@meta.data$seurat_clusters == cl,])])
              t = paste0(cond, "_cl.", cl)

              obj = COTAN(raw = raw)
              obj = initializeMetaDataset(
                      obj, GEO="", sequencingMethod = " ",
                      sampleCondition = paste("temp.", cond,
                                              "clustered by Seurat"))

              #obj = readRDS(paste0(out_dir, t, ".cotan.RDS"))

              cells_to_rem = getCells(obj)[which(getCellsSize(obj) == 0)]
              obj = dropGenesCells(obj, cells = cells_to_rem)

              print(paste0("Condition ", t))
              #--------------------------------------
              n_cells = getNumCells(obj)
              print(paste("n cells", n_cells))

              obj <- as(obj, "scCOTAN")

              if (n_cells <= 10) {
                print("Cell cluster too small!")
                to_rec = c(to_rec,colnames(raw))
              }else{

                numIter = 1

                list[obj, data] <- clean(obj)
                plots <- cleanPlots(obj, data[["pcaCells"]], data[["D"]])

                #---------- run this when B cells are to be removed
                pdf(paste0(out_dir, "/", t, "_", numIter, "_plots.pdf"))

                plot(plots[["pcaCells"]])
                plot(plots[["genes"]])
                plot(plots[["UDE"]])
                plot(plots[["nu"]])

                rm(plots)
                gc()

                #dev.off()

                obj <- as(obj, "scCOTAN")
                obj = cotan_analysis(obj, cores = cores)
                # saving the structure
                #saveRDS(obj,file = paste0(out_dir, t, ".cotan.RDS"))

                # COEX evaluation and storing
                gc()
                obj = get.coex(obj)

                gc()

                GDI_data_wt1 = get.GDI(obj)

                #obj = readRDS(paste0(out_dir, cond, "_cl.", cl, ".cotan.RDS"))
                #GDI_data_wt1 = get.GDI(obj)

                #Test if the number of genes with GDI > 1.5 is more than 1%
                if (dim(GDI_data_wt1[GDI_data_wt1$GDI >= 1.5,])[1]/dim(GDI_data_wt1)[1] > 0.01) {
                  print(paste("Cluster", cl, "too high GDI! Recluster!"))
                  cells_to_cluster = colnames(obj@raw)
                  write.csv(cells_to_cluster,
                            file = paste0(out_dir, "to_recluster_",
                                          cond, "_cl.", cl, ".csv"))
                  to_rec = c(to_rec,cells_to_cluster)
                }


                genes.to.label = rownames(GDI_data_wt1[order(GDI_data_wt1$GDI,decreasing = T),][1:20,])

                #pdf(paste0(out_dir, cond, "_cl.", cl, ".GDI_plots.pdf"), onefile=TRUE)
                #my.plots[[(cl+1)]] <-
                plot(plot_GDI(obj, genes = list("top 20 GDI genes"=genes.to.label)))

                dev.off()
                #graphics.off()
                rm(obj)
                gc()
              }

            }


            cells.to.add = rownames(data.seurat@meta.data[data.seurat@meta.data$seurat_clusters == "singleton",])
            if (length(cells.to.add)>0) {
              to_rec = c(to_rec,cells.to.add)
            }
            return(to_rec)
          }
)

