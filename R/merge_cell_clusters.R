#' merge_cell.clusters
#'
#'This function takes in input a COTAN object with an already homogeneous clusters
#'and, through the cosine distance and the hclust (with ward.D2 method), checks it merging
#'two leaf clusters it is formed a still homogeneous cluster (this is done iteratively).
#'All structures are saved on disk. 
#'
#' @param obj COTAN object
#' @param cond sample condition name
#' @param cores number cores used
#' @param srat the Seurat object for the same data set (the name should be "cond")
#' @param out_dir path to a directory for output path to the directory in which there is the Seurat object
#' @param GEO GEO or other data set code
#' @param sc.method scRNAseq method
#' @param markers a list of marker genes. Default NULL
#' 
#' @importFrom Matrix t
#' @importFrom Matrix rowSums
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats as.dendrogram
#' @importFrom dendextend get_nodes_attr
#' @importFrom stringr str_remove
#' @import ggplot2
#'
#' @return a COTAN object 
#' @export
#'
#' @examples
setGeneric("merge_cell.clusters", function(obj,cond,cores=1,srat,out_dir ,GEO,
                                           sc.method,#mt = FALSE, mt_prefix="^mt", 
                                           markers = NULL) standardGeneric("merge_cell.clusters"))
#' @rdname merge_cell.clusters
setMethod("merge_cell.clusters","scCOTAN",
         function(obj,cond,cores,srat,out_dir ,GEO,
                   sc.method,#mt, #mt_prefix, 
                  markers){
            
           srat <- readRDS(paste0(out_dir,srat))
           
           #Check if there are NA left clusters
           
           if (any(is.na(obj@clusters))) {
             print("Problems! Some NA left!")
             break
           }
           
           
           ###############
            cl1_not_mergiable <- c()
            cl2_not_mergiable <- c()
            cl1_not_mergiable_old <- 1
            round <- 0
            while (!identical(cl1_not_mergiable,cl1_not_mergiable_old)) {
              round <- round + 1
              print(paste0("Start merging smallest clusters: round ", round))
              cl1_not_mergiable_old <- cl1_not_mergiable
              
              #merge small cluster based on distances
              cluster_data <- obj@cluster_data
              
              ######## This is the best: cosine dissimilarity
              Matrix <- as.matrix(Matrix::t(cluster_data))
              sim <- Matrix / sqrt(Matrix::rowSums(Matrix * Matrix))
              sim <- sim %*% Matrix::t(sim)
              D_sim <- stats::as.dist(1 - sim)
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
              
              dir <- paste0(out_dir,"cond/","leafs_merge/")
              if(!file.exists(paste0(out_dir,"cond"))){
                dir.create(paste0(out_dir,"cond"))
              }
              
              
              if(!file.exists(dir)){
                dir.create(dir)
              }
              
              
              p = 1
              while (p < length(id)) {
                p1 <- p
                p2 <- p+1
                p <- p2+1
                
                cl.1 <- stringr::str_remove(dendextend::get_nodes_attr(dend,"label",id = id[p1]),pattern = "cl.")
                cl.2 <- stringr::str_remove(dendextend::get_nodes_attr(dend,"label",id = id[p2]),pattern = "cl.")
                
                
                cond.merge <- paste0("merge_cl_",cl.1,"_",cl.2)
                print(cond.merge)
                if(cl.1 %in% cl1_not_mergiable){
                  print(paste0("Clusters ", cl.1," ", cl.2,"already analyzed and not mergiable: skip."))
                  next
                }
                
                mat <- srat@assays$RNA@counts[,colnames(srat@assays$RNA@counts) %in% 
                                                rownames(srat@meta.data[srat@meta.data$cotan %in% c(cl.1,cl.2),])]
                
                merged.obj <- automatic.COTAN.object.creation(df = as.data.frame(mat),out_dir = out_dir ,GEO = GEO
                                                              ,sc.method=sc.method,cond = cond.merge,cores = cores,
                                                              #mt = mt,mt_prefix = mt_prefix  
                )
                GDI_data_wt1 = get.GDI(merged.obj)
                
                #Test if the number of genes with GDI > 1.5 is more than 1%
                if (dim(GDI_data_wt1[GDI_data_wt1$GDI >= 1.5,])[1]/dim(GDI_data_wt1)[1] > 0.01) {
                  print(paste("Clusters ",cl.1, " and ",cl.2," too high GDI!",sep = " "))
                  cl1_not_mergiable <- c(cl1_not_mergiable,cl.1)
                  cl2_not_mergiable <- c(cl2_not_mergiable,cl.2)
                  #cells_to_cluster = colnames(merged.obj@raw)
                  #write.csv(cells_to_cluster, file = paste(out_dir,"to_recluster_",cond,"_cl.",cl,".csv",sep = ""))
                  #to_rec = c(to_rec,cells_to_cluster)
                }else{
                  print(paste("Clusters ",cl.1, " and ",cl.2," can be merged.",sep = " "))
                  write.csv(colnames(mat),paste0(dir,"merged_clusters_",cl.1,"_",cl.2,"cell_ids.csv"))
                  min.cl <- min(as.numeric(c(cl.1,cl.2)))
                  max.cl <- max(as.numeric(c(cl.1,cl.2)))
                  srat@meta.data[srat@meta.data$cotan == max.cl,]$cotan <- min.cl
                }
                
                
                GDI_data_wt1$color = "normal"
                genes.to.label = GDI_data_wt1[order(GDI_data_wt1$GDI,decreasing = T),][1:20,]
                
                #genes.to.label = rbind(genes.to.label,GDI_data_wt1[more_genes,])
                genes.to.label$color = "dif"
                #genes.to.label[more_genes,]$color = "mk"
                GDI_data_wt1[rownames(genes.to.label),]$color = "dif"
                #GDI_data_wt1[more_genes,]$color = "mk"
                
                # from here it is just a plot example
                mycolours <- c("dif" = "#3C5488B2","normal"="#F39B7FE5","hk"="#7E6148B2","mk"="#E64B35B2")
                
                f1 = ggplot( subset(GDI_data_wt1,!rownames(GDI_data_wt1) %in% unique(rownames(genes.to.label))),  aes(x=sum.raw.norm, y=GDI)) +  geom_point(alpha = 0.4, color = "#8491B4B2", size=2)
                
                si=12
                GDI_plot_wt1 = f1 + geom_point(data = subset(GDI_data_wt1,rownames(GDI_data_wt1) %in% c(rownames(genes.to.label),"Lhx1os","5330434G04Rik")),aes(x=sum.raw.norm, y=GDI, color=color),alpha = 1, size=2)+
                  #geom_hline(yintercept=quantile(GDI$GDI)[4], linetype="dashed", color = "darkblue") +
                  #geom_hline(yintercept=quantile(GDI$GDI)[3], linetype="dashed", color = "darkblue") +
                  geom_hline(yintercept=1.5, linetype="dotted", color = "#3C5488B2", size= 0.5) +
                  scale_color_manual("color", values = mycolours)  +
                  scale_fill_manual("color", values = mycolours)  +
                  xlab("log normalized counts")+ylab("GDI")+
                  geom_label_repel(data =genes.to.label , aes(x=sum.raw.norm, y=GDI, label = rownames(genes.to.label),
                                                              fill=color),
                                   label.size = NA,
                                   alpha = 0.5,
                                   direction = "both",
                                   na.rm=TRUE,
                                   seed = 1234) +
                  theme(axis.text.x = element_text(size = si, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                        axis.text.y = element_text( size = si, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                        axis.title.x = element_text( size = si, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                        axis.title.y = element_text( size = si, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
                        legend.title = element_blank(),
                        legend.text = element_text(color = "#3C5488FF",face ="italic" ),
                        legend.position = "none") + ggtitle(paste(cond.merge,getNumCells(obj),sep = " "))
                
                
                pdf(paste(dir,cond.merge, ".GDI_plots.pdf", sep = ""), onefile=TRUE)
                plot(GDI_plot_wt1)
                graphics.off()
                
                rm(merged.obj)
                gc()
              }
              
              saveRDS(srat,paste(out_dir,"Seurat_obj_",cond,"_with_cotan_clusters_merged.RDS",sep = ""))
              #srat <- readRDS(paste(out_dir_root,"Seurat_obj_",cond,"_with_cotan_clusters_merged.RDS",sep = ""))
              gc()
              
              
              # Update the homogenus cluster array in cotan object
              obj.clusters.new <- srat@meta.data$cotan
              names(obj.clusters.new) <- rownames(srat@meta.data)
              if(all(colnames(obj@raw) %in% names(obj.clusters.new))){
                obj@clusters <- obj.clusters.new
              }else{
                print("Problem! not all all(colnames(obj@raw) %in% names(obj.clusters.new)")
                break
              }
              
              
              # New DEA on clusters
              clusters.names = unique(obj@clusters)[!is.na(unique(obj@clusters))]
              list.clusters = list(names(obj@clusters[obj@clusters %in% clusters.names[1]]))
              names(list.clusters)=clusters.names[1]
              for (c in c(2:length(clusters.names))) {
                tmp = list(names(obj@clusters[obj@clusters %in% clusters.names[c]]))
                names(tmp)= clusters.names[c]
                list.clusters = c(list.clusters,tmp)
              }
              
              rm(srat)
              gc()
              
              
              obj_list = DEA_on_clusters(obj,list.clusters)
              gc()
              
              srat <- readRDS(paste(out_dir,"Seurat_obj_",cond,"_with_cotan_clusters_merged.RDS",sep = ""))
              obj = obj_list[[1]]
              
              p_value = obj_list[[2]]
              
              rm(obj_list)
              gc()
              
              write.csv(p_value,file = paste(out_dir,cond,"/p_values_clusters_merged.csv", sep = ""))
              write.csv(obj@cluster_data,file = paste(out_dir,cond,"/coex_clusters_merged.csv", sep = ""))
              if (!is.null(markers)) {
                write.csv(p_value[unlist(markers),],file = paste(out_dir,cond,"/p_values_clusters_merged_markers.csv", sep = ""))
                write.csv(obj@cluster_data[unlist(markers),],file = paste(out_dir,cond,"/coex_clusters_merged_markers.csv", sep = ""))
                
              }
              
            }
            saveRDS(obj,file = paste(out_dir,cond,"_merged_cotan.RDS", sep = ""))
            
            return(obj)
            
          }
)
