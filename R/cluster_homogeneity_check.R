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
setMethod("cluster_homogeneity_check","COTAN",
          function(obj,cells,out_dir, cores,code){

            options(ggrepel.max.overlaps = Inf)


            raw <- obj@raw[,colnames(obj@raw) %in% cells]

            obj <- COTAN(raw = raw)
            obj <- initializeMetaDataset(obj, GEO = "", sequencingMethod = " ",
                                         sampleCondition = "temp.clustered")

            cells_to_rem <- getCells(obj)[which(getCellsSize(obj) == 0)]
            obj <- dropGenesCells(obj, cells = cells_to_rem)

            #--------------------------------------
            print(paste("n cells", getNumCells(obj)))

            list[obj, data] <- clean(obj)
            plots <- cleanPlots(obj, data[["pcaCells"]], data[["D"]])

            #---------- run this when B cells are to be removed

            obj <- estimateDispersion(obj, cores = cores)
            gc()

            # COEX evaluation and storing
            obj <- calculateCoex(obj)
            gc()

            GDI_data_wt1 <- calculateGDI(obj)

            # Plots
            genes.to.label = rownames(GDI_data_wt1[order(GDI_data_wt1$GDI,decreasing = T),][1:10,])
            pdf(paste0(out_dir, "/", code, "_plots.pdf"))

            plot(plots[["pcaCells"]])
            plot(plots[["genes"]])
            plot(plots[["UDE"]])
            plot(plots[["nu"]])
            obj <- as(obj, "scCOTAN")
            plot(plot_GDI(obj,GDI.df = GDI_data_wt1, genes = list("top 20 GDI genes"=genes.to.label)))

            dev.off()
            #graphics.off()
            rm(obj)
            rm(plots)
            gc()

            #Test if the number of genes with GDI > 1.5 is more than 1%
            if (dim(GDI_data_wt1[GDI_data_wt1$GDI >= 1.5,])[1]/dim(GDI_data_wt1)[1] > 0.01) {
              print(paste("Cluster", code, "too high GDI! Recluster!"))
              write.csv(cells, file = paste0(out_dir, "to_recluster_cl_", code, ".csv"))
              return(cells)
            }else{
              print(paste("Cluster", code, "homogeneous."))
            }

          }
)

