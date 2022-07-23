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
#' @import ggplot2
#' @return an array of cells that need to be re-clustered
#' @export
#'
#' @examples
setGeneric("cluster_homogeneity", function(data.seurat,data.raw,out_dir, cond, cores=1) standardGeneric("cluster_homogeneity"))
#' @rdname cluster_homogeneity
setMethod("cluster_homogeneity","Seurat",
          function(data.seurat,data.raw,out_dir, cond, cores){
            options(ggrepel.max.overlaps = Inf)
            # This detect the number of clusters
            mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
            my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                             axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                             axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                             axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))

            tot.clusters = max(as.numeric(as.vector(unique(data.seurat@meta.data$seurat_clusters[!data.seurat@meta.data$seurat_clusters == "singleton" ]))))
            print(paste("Number of clusters:",tot.clusters,sep=" "))
            to_rec = NA
            for (cl in c(0:tot.clusters)) { #0

              print(paste("Working on cluster",cl,"out of",tot.clusters,sep=" " ))
              raw = as.data.frame(data.raw[,rownames(data.seurat@meta.data[data.seurat@meta.data$seurat_clusters == cl,])])
              t = paste(cond,"_cl.",cl,sep = "")


              obj = new("scCOTAN",raw = raw)
              obj = initRaw(obj,GEO="" ,sc.method=" ",cond = paste("temp.",cond, "clustered by Seurat", sep=" "))

              #obj = readRDS(paste(out_dir,t,".cotan.RDS", sep = ""))
              genes_to_rem = get.genes(obj)[grep('^mt', get.genes(obj))]
              cells_to_rem = names(get.cell.size(obj)[which(get.cell.size(obj) == 0)])
              obj = drop.genes.cells(obj,genes_to_rem,cells_to_rem )

              print(paste("Condition ",t,sep = ""))
              #--------------------------------------
              n_cells = length(get.cell.size(object = obj))
              print(paste("n cells", n_cells, sep = " "))
              if (n_cells <= 10) {
                print("Cell cluster too small!")
                to_rec = c(to_rec,colnames(raw))
              }else{

                n_it = 1

                #obj@raw = obj@raw[rownames(cells), colnames(cells)]

                ttm = clean(obj)
                #ttm = clean.sqrt(obj, cells)

                obj = ttm$object

                ttm$pca.cell.2

                #---------- run this when B cells are to be removed
                pdf(paste(out_dir,"cleaning/",t,"_",n_it,"_plots.pdf", sep = ""))

                plot(ttm$pca.cell.2)

                plot(
                  ggplot(ttm$D, aes(x=n,y=means)) + geom_point() +
                    geom_text_repel(data=subset(ttm$D, n > (max(ttm$D$n)- 15) ), aes(n,means,label=rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
                                    nudge_y      = 0.05,
                                    nudge_x      = 0.05,
                                    direction    = "x",
                                    angle        = 90,
                                    vjust        = 0,
                                    segment.size = 0.2)+
                    ggtitle("B cell group genes mean expression")+my_theme +
                    theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 5,hjust = 0.02 ))
                )
                dev.off()


                nu_est = round(get.nu(object = obj), digits = 7)

                p<-ggplot(ttm$pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))


                pdf(paste(out_dir,"cleaning/",t,"_plots_PCA_efficiency_colored.pdf", sep = ""))
                # or tiff("plot.tiff")
                plot(
                  p+geom_point(size = 1,alpha= 0.8)+
                    scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                                          midpoint = log(mean(nu_est)),name = "ln(nu)")+
                    ggtitle("Cells PCA coloured by cells efficiency") +
                    my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                                      legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                                      legend.text = element_text(color = "#3C5488FF", size = 11),
                                      legend.key.width = unit(2, "mm"),
                                      legend.position="right")
                )
                dev.off()

                nu_df = data.frame("nu"= sort(obj@nu), "n"=c(1:length(obj@nu)))

                pdf(paste(out_dir,"cleaning/",t,"_plots_efficiency.pdf", sep = ""))
                plot(ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1)+my_theme )
                dev.off()

                obj = hk_genes(obj)

                obj = cotan_analysis(obj, cores = cores)
                # saving the structure
                #saveRDS(obj,file = paste(out_dir,t,".cotan.RDS", sep = ""))

                # COEX evaluation and storing
                gc()
                obj = get.coex(obj)

                gc()
                # saving the structure
                #saveRDS(obj,file = paste(out_dir,t,".cotan.RDS", sep = ""))

                #}

                #for (cl in c(0:tot.clusters)) {#tot.clusters
                #obj = readRDS(paste(out_dir,cond,"_cl.", cl,".cotan.RDS",sep = ""))
                GDI_data_wt1 = get.GDI(obj)

                #obj = readRDS(paste(out_dir,cond,"_cl.", cl,".cotan.RDS",sep = ""))
                #GDI_data_wt1 = get.GDI(obj)

                #Test if the number of genes with GDI > 1.5 is more than 1%
                if (dim(GDI_data_wt1[GDI_data_wt1$GDI >= 1.5,])[1]/dim(GDI_data_wt1)[1] > 0.01) {
                  print(paste("Cluster",cl,"too high GDI!Recluster!",sep = " "))
                  cells_to_cluster = colnames(obj@raw)
                  write.csv(cells_to_cluster, file = paste(out_dir,"to_recluster_",cond,"_cl.",cl,".csv",sep = ""))
                  to_rec = c(to_rec,cells_to_cluster)
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
                        legend.position = "none")  +ggtitle(paste("Cluster n.",cl,"N.cells",obj@n_cells,sep = " "))


                pdf(paste(out_dir,cond,"_cl.",cl, ".GDI_plots.pdf", sep = ""), onefile=TRUE)
                #my.plots[[(cl+1)]] <-
                plot(GDI_plot_wt1)
                graphics.off()
                #graphics.off()
              }

            }

            cells.to.add = rownames(data.seurat@meta.data[data.seurat@meta.data$seurat_clusters == "singleton",])
            if (length(cells.to.add)>0) {
              to_rec = c(to_rec,cells.to.add)
            }
            return(to_rec)
          }
)

