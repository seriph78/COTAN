#' cluster_homogeneity_check
#'
#' # Function that takes a COTAN object and a set of cells to check if that 
#' cluster is or not homogeneous.
#' It runs COTAN check if the GDI is lower than 1.5 for the 99% of genes. 
#' If it is too high,
#' the cluster is not uniform and so it return the cells in an array to cluster 
#' again.
#' @param obj a COTAN object 
#' @param cells an array of cell codes
#' @param out_dir path to the output directory (to print a pdf with the plots
#' or store the cell codes if the cluster need to be split)
#' @param code string defining the condition/cluster number
#' @param cores number of cores used
#' 
#' @import ggrepel
#' @import ggplot2
#' @return an array of cells that need to be re-clustered or nothing
#' @export
#'
#' @examples
setGeneric("cluster_homogeneity_check", function(obj,cells,out_dir, cores=1,code) standardGeneric("cluster_homogeneity_check"))
#' @rdname cluster_homogeneity_check
setMethod("cluster_homogeneity_check","scCOTAN",
          function(obj,cells,out_dir, cores,code){
            
            options(ggrepel.max.overlaps = Inf)
            
              
            raw <- obj@raw[,colnames(obj@raw) %in% cells]
  
            obj <- new("scCOTAN",raw = raw)
            obj <- initRaw(obj,GEO="" ,sc.method=" ",cond = "temp.clustered")
            
            cells_to_rem <- getCells(obj)[which(getCellsSize(obj) == 0)]
            obj <- drop.genes.cells(obj, cells = cells_to_rem )
            
            #--------------------------------------
            print(paste("n cells", getNumCells(obj), sep = " "))
            
          
            ttm <- clean(obj)
              
            obj <- ttm$object
              
            ttm$pca.cell.2
              
            #---------- run this when B cells are to be removed
            
            obj <- cotan_analysis(obj, cores = cores)
            # COEX evaluation and storing
            gc()
            obj <- get.coex(obj)
            gc()
            
            GDI_data_wt1 <- get.GDI(obj)
            
            # Plots
            genes.to.label = rownames(GDI_data_wt1[order(GDI_data_wt1$GDI,decreasing = T),][1:10,])
            pdf(paste(out_dir,"/", code ,"_plots.pdf", sep = ""))
            
            plot(ttm$pca.cell.2)
            plot(ttm$genes.plot)
            plot(ttm$UDE.plot)
            nu_df = data.frame("nu"= sort(obj@nu), "n"=c(1:length(obj@nu)))
            plot(ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1) )
            plot(plot_GDI(obj,GDI.df = GDI_data_wt1, genes = list("top 20 GDI genes"=genes.to.label)))
            
            dev.off()
            #graphics.off()
            rm(obj)
            rm(ttm)
            gc()
            
            #Test if the number of genes with GDI > 1.5 is more than 1%
            if (dim(GDI_data_wt1[GDI_data_wt1$GDI >= 1.5,])[1]/dim(GDI_data_wt1)[1] > 0.01) {
              print(paste("Cluster",code,"too high GDI!Recluster!",sep = " "))
              write.csv(cells, file = paste(out_dir,"to_recluster_cl_",code,".csv",sep = ""))
              return(cells)
            }else{
              print(paste("Cluster",code,"homogeneous.",sep = " "))
            }
            
          }
)

