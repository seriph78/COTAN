#' merge_cell.clusters
#'
#' This function takes in input a COTAN object with an already homogeneous clusters
#' and, through the cosine distance and the hclust (with ward.D2 method), checks if merging
#' two leaf clusters will form a still homogeneous cluster (this is done iteratively).
#' All structures are saved on disk.
#'
#' @param obj COTAN object
#' @param cond sample condition name
#' @param cores number cores used
#' @param out_dir path to a directory for output path to the directory in which there is the Seurat object
#' @param GEO GEO or other data set code
#' @param sc.method scRNAseq method
#' @param markers a list of marker genes. Default NULL
#'
#' @return a COTAN object
#'
#' @importFrom Matrix t
#' @importFrom Matrix rowSums
#'
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats as.dendrogram
#'
#' @importFrom dendextend get_nodes_attr
#'
#' @importFrom stringr str_remove
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
#' @examples
#'
#' @rdname merge_cell.clusters
setGeneric("merge_cell.clusters", function(obj, cond,
                                           cores = 1,  out_dir ,GEO,
                                           sc.method, #mt = FALSE, mt_prefix="^mt",
                                           markers = NULL) standardGeneric("merge_cell.clusters"))
#' @rdname merge_cell.clusters
setMethod("merge_cell.clusters","COTAN",
         function(obj, cond, cores, out_dir, GEO,
                  sc.method, #mt, #mt_prefix,
                  markers){

           #######srat <- readRDS(paste0(out_dir,srat))

           #Check if there are NA left clusters

           #if (any(is.na(obj@clusters))) {
          #   print("Problems! Some NA left!")
          #   break
          # }


           ###############
            cl1_not_mergeable <- c()
            cl2_not_mergeable <- c()
            cl1_not_mergeable_old <- 1
            round <- 0
            while (!identical(cl1_not_mergeable,cl1_not_mergeable_old)) {
              round <- round + 1
              print(paste0("Start merging smallest clusters: round ", round))
              cl1_not_mergeable_old <- cl1_not_mergeable

              #merge small cluster based on distances
              D_sim <- cosineDissimilarity(as.matrix(getClustersCoex(obj)[[length(getClustersCoex(obj))]]))

              tree <- stats::hclust(D_sim,method = "ward.D2")
              #plot(tree)
              dend <- stats::as.dendrogram(tree)


              ############### This part check if any little two pair of leaf clusters could be merged
              # based on the tree

              id <- NA
              for (i in c(1:length(dendextend::get_nodes_attr(dend,"members")))) {
                if(dendextend::get_nodes_attr(dend,"members")[i] == 2){
                  id <- c(id,i+1,i+2)
                }
              }
              id <- id[2:length(id)]
              print(paste0("Created leafs id form marging: ",
                           paste(dendextend::get_nodes_attr(dend,"label",id = id ),collapse=" ")))

              if (!file.exists(file.path(out_dir, "cond"))) {
                dir.create(file.path(out_dir, "cond"))
              }

              dir <- file.path(out_dir, "cond", "leafs_merge")
              if (!file.exists(dir)) {
                dir.create(dir)
              }

              meta.cells <- getClusterizationData(obj)[[1]]

              p = 1
              while (p < length(id)) {
                p1 <- p
                p2 <- p+1
                p <- p2+1

                cl.1 <- stringr::str_remove(dendextend::get_nodes_attr(dend,"label",id = id[p1]),pattern = "cl.")
                cl.2 <- stringr::str_remove(dendextend::get_nodes_attr(dend,"label",id = id[p2]),pattern = "cl.")


                cond.merge <- paste0("merge_cl_", cl.1, "_", cl.2)
                print(cond.merge)
                if(cl.1 %in% cl1_not_mergeable){
                  print(paste0("Clusters ", cl.1, " ", cl.2, " already analyzed and not mergeable: skip."))
                  next
                }

                #mat <- srat@assays$RNA@counts[,colnames(srat@assays$RNA@counts) %in%
                 #                               names(obj@clusters[obj@clusters %in% c(cl.1,cl.2)])]
                                                #rownames(srat@meta.data[srat@meta.data$cotan %in% c(cl.1,cl.2),])]



                cells.to.merge <- names(meta.cells[meta.cells %in% c(cl.1,cl.2)])
                merged.obj <- automaticCOTANObjectCreation(raw = getRawData(obj)[,cells.to.merge],
                                                           GEO = GEO,
                                                           sequencingMethod = sc.method,
                                                           sampleCondition = cond.merge,
                                                           cores = cores,
                                                           saveObj = FALSE,
                                                           outDir = out_dir)

                GDI.data = calculateGDI(merged.obj)

                #Test if the number of genes with GDI > 1.5 is more than 1%
                if (dim(GDI.data[GDI.data$GDI >= 1.5,])[1]/dim(GDI.data)[1] > 0.01) {
                  print(paste("Clusters", cl.1, "and", cl.2, "too high GDI!"))
                  cl1_not_mergeable <- c(cl1_not_mergeable,cl.1)
                  cl2_not_mergeable <- c(cl2_not_mergeable,cl.2)
                  #cells_to_cluster = colnames(merged.obj@raw)
                  #write.csv(cells_to_cluster,
                  #          file = paste0(out_dir, "to_recluster_",
                  #                        cond, "_cl.", cl, ".csv"))
                  #to_rec = c(to_rec,cells_to_cluster)
                }else{
                  print(paste("Clusters", cl.1, "and", cl.2, "can be merged."))
                #  write.csv(colnames(mat),paste0(dir, "merged_clusters_",
                #                                 cl.1, "_", cl.2, "cell_ids.csv"))
                  min.cl <- min(as.numeric(c(cl.1,cl.2)))
                  max.cl <- max(as.numeric(c(cl.1,cl.2)))
                  #srat@meta.data[srat@meta.data$cotan == max.cl,]$cotan <- min.cl
                  #obj@clusters[obj@clusters == max.cl] <- min.cl
                  meta.cells[meta.cells == max.cl] <- min.cl
                }


                GDI.data$color = "normal"
                genes.to.label = GDI.data[order(GDI.data$GDI,decreasing = T),][1:20,]

                #genes.to.label = rbind(genes.to.label,GDI.data[more_genes,])
                genes.to.label$color = "dif"
                #genes.to.label[more_genes,]$color = "mk"
                GDI.data[rownames(genes.to.label),]$color = "dif"
                #GDI.data[more_genes,]$color = "mk"

                # from here it is just a plot example
                mycolours <- c("dif" = "#3C5488B2","normal"="#F39B7FE5","hk"="#7E6148B2","mk"="#E64B35B2")

                GDI_plot <- ggplot(subset(GDI.data, !rownames(GDI.data) %in% unique(rownames(genes.to.label))),
                                       aes(x = sum.raw.norm, y = GDI)) +
                                geom_point(alpha = 0.4, color = "#8491B4B2", size = 2) +
                                geom_point(data = subset(GDI.data, rownames(GDI.data) %in% rownames(genes.to.label)),
                                                         aes(x = sum.raw.norm, y = GDI, color = color), alpha = 1, size = 2) +
                                # geom_hline(yintercept = quantile(GDI$GDI)[4], linetype = "dashed", color = "darkblue") +
                                # geom_hline(yintercept = quantile(GDI$GDI)[3], linetype = "dashed", color = "darkblue") +
                                geom_hline(yintercept = 1.5, linetype="dotted", color = "#3C5488B2", linewidth = 0.5) +
                                scale_color_manual("color", values = mycolours)  +
                                scale_fill_manual("color", values = mycolours)  +
                                xlab("log normalized counts") +
                                ylab("GDI") +
                                geom_label(data = genes.to.label,
                                                 aes(x = sum.raw.norm, y = GDI,
                                                     label = rownames(genes.to.label), fill = color),
                                                 label.size = NA, alpha = 0.5,
                                                 #direction = "both",
                                           na.rm = TRUE,
                                           #seed = 1234
                                           ) +
                                ggtitle(paste(cond.merge,"- cell number: ", getNumCells(merged.obj))) +
                                plotTheme("GDI", textSize = 12)

                pdf(file.path(dir, paste0(cond.merge, ".GDI_plots.pdf")), onefile=TRUE)
                plot(GDI_plot)

                graphics.off()

                rm(merged.obj)
                gc()
              }

              #saveRDS(srat,paste0(out_dir, "Seurat_obj_",
              #                    cond, "_with_cotan_clusters_merged.RDS"))
              #srat <- readRDS(paste0(out_dir_root, "Seurat_obj_",
              #                       cond, "_with_cotan_clusters_merged.RDS"))
              gc()


              # Update the homogenus cluster array in cotan object

              #obj.clusters.new <- set_names(srat@meta.data$cotan, rownames(srat@meta.data))
              #if(all(getCells(obj) %in% names(obj.clusters.new))){
              #  obj@clusters <- obj.clusters.new
              #}else{
              #  print("Problem! not all(getCells(obj) %in% names(obj.clusters.new)")
              #  break
              #}

              # New DEA on clusters
              #obj <- addClusterization(obj,clName = "first_merged",clusters = meta.cells)

              list[coexDF, pvalDF] <- DEAOnClusters(obj, clusterization = meta.cells)

              obj <- addClusterization(obj, clName = paste0("merged_", round),
                                       clusters = meta.cells, coexDF = coexDF)

              #clusters.names = unique(obj@clusters)[!is.na(unique(obj@clusters))]
              #list.clusters = list(names(obj@clusters[obj@clusters %in% clusters.names[1]]))
              #names(list.clusters)=clusters.names[1]
              #for (c in c(2:length(clusters.names))) {
              #  tmp = list(names(obj@clusters[obj@clusters %in% clusters.names[c]]))
              #  names(tmp)= clusters.names[c]
              #  list.clusters = c(list.clusters,tmp)
              #}

              #rm(srat)
              gc()



              #srat <- readRDS(paste0(out_dir, "Seurat_obj_", cond,
               #                      "_with_cotan_clusters_merged.RDS"))
              #obj = obj_list[[1]]

              #p_value = obj_list[[2]]

              #rm(obj_list)
              #gc()

              write.csv(pvalDF, file =
                          file.path(out_dir, cond, "p_values_clusters_merged.csv"))

              #write.csv(obj@cluster_data,
              #          file = paste0(out_dir, cond, "/coex_clusters_merged.csv"))
              #if (!is.null(markers)) {
               # write.csv(p_value[unlist(markers),],
              #            file = paste0(out_dir, cond, "/p_values_clusters_merged_markers.csv"))
              #  write.csv(obj@cluster_data[unlist(markers),],
              #            file = paste0(out_dir, cond, "/coex_clusters_merged_markers.csv"))
              #}

            }
            #saveRDS(obj, file = paste0(out_dir, cond, "_merged_cotan.RDS"))

            return(obj)

          }
)
