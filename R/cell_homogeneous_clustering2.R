#' cell_homogeneous_clustering
#'
#'This function imports a Seurat object from in_dir with a file name specified by data set_name.
#'From this, and using COTAN functions, it splits the data set until all clusters are homogeneous by GDI
#'(calling the cluster_homogeneity function).
#' @param cond a string specifying the experiment condition
#' @param out_dir an existing directory for the analysis output. In this directory is saved also the final Seurat object.
#' @param in_dir directory in which there is the Seurat object
#' @param cores number of cores (NB for windows system no more that 1 can be used)
#' @param dataset_name Seurat object file name
#' @param GEO GEO or other dataset official code
#' @param sc.method single cell method used fot the experiment 
#' @import Seurat
#' @importFrom stringr str_detect
#' @return the scCOTAN object and it saves, in the out_dir, the Seurat elaborated data file
#' @export 
#'
#' @examples
setGeneric("cell_homogeneous_clustering", function(cond,out_dir,in_dir,cores=1, dataset_name,
                                                   GEO, sc.method )
  standardGeneric("cell_homogeneous_clustering"))
#' @rdname cell_homogeneous_clustering
setMethod("cell_homogeneous_clustering","character",
 function(cond,out_dir,in_dir,cores, dataset_name){
  
  out_dir_root <- paste0(out_dir,"/",cond,"/")
  if(!file.exists(out_dir_root)){
    dir.create(file.path(out_dir_root))
  }
  #Step 1
  obj <- readRDS(paste0(in_dir,dataset_name))
  srat <- CreateSeuratObject(counts = as.data.frame(obj@raw), project = cond, 
                             min.cells = 3, min.features = 4)
  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(srat)
  srat <- ScaleData(srat, features = all.genes)
  srat <- RunPCA(srat, features = VariableFeatures(object = srat))

  out_dir <- out_dir_root
  
  srat <- FindNeighbors(srat, dims = 1:25)
  srat <- FindClusters(srat, resolution = 0.5,algorithm = 2)
  
  srat <- RunUMAP(srat,umap.method = "uwot",metric = "cosine", 
                        dims = 1:min(c(50,(nrow(srat@meta.data)-1))))
  print(paste0("PDF UMAP! ", out_dir, "pdf_umap.pdf"))
  pdf(paste0(out_dir,"pdf_umap.pdf"))
  plot(DimPlot(srat, reduction = "umap",label = TRUE))
  dev.off()
  gc()
  
  # Step 2
  df.cell.clustering <- as.data.frame(matrix(nrow = ncol(obj@raw),ncol = 1))
  rownames(df.cell.clustering) <- colnames(obj@raw)
  colnames(df.cell.clustering) <- "cl_1"
  
  to_recluster_new <- NA
  for (cl in unique(srat@meta.data$seurat_clusters)) {
  print(cl)
    cells <- rownames(srat@meta.data[srat@meta.data$seurat_clusters == cl,])
    to_recluster_new <- c(to_recluster_new,cluster_homogeneity_check(obj = obj,
                                                                     cells = cells,
                                                                     out_dir = out_dir,
                                                                     cores = cores,
                                                                     code = cl) 
    )
  }
    
  homog.clusters <- to_recluster_new[stringr::str_detect(to_recluster_new,pattern = "Cluster")]
  print(paste(homog.clusters,collapse = " "))
  to_recluster_new <- to_recluster_new[stringr::str_detect(to_recluster_new,pattern = "Cluster",negate = T)]
  to_recluster_new <- to_recluster_new[2:length(to_recluster_new)]
  
  #clusters.to.recluster <- unique(srat@meta.data[rownames(srat@meta.data) %in% to_recluster_new,]$seurat_clusters)
  
  tmp_meta <- srat@meta.data[!rownames(srat@meta.data) %in% to_recluster_new,]
  
  # To save the already defined clusters
  df.cell.clustering[rownames(tmp_meta),"cl_1"] <- as.vector(tmp_meta$seurat_clusters)
  
  if (!sum(is.na(df.cell.clustering$cl_1)) == length(to_recluster_new)) {
    print("Problems!")
  }
  
  saveRDS(df.cell.clustering,paste(out_dir_root,"cluster_df_",cond,".RDS",sep = ""))
  write.csv(to_recluster_new,paste(out_dir_root,"to_recluster_",cond,"tot.csv",sep = ""),row.names = F)
  
  print("First part finished! Starting the iterative part.")
  
  round <- 0
  singletons <- 0
  #n_cells <- length(obj@clusters)
  to_recluster_old <- NA
  number.cls.old <- 0 #Just to start from a number higher than any possible real situation
  
  while (sum(is.na(df.cell.clustering$cl_1)) > 50 & round < 25) {
    #data.raw = Read10X(get(paste0("raw.",cond)))
    gc()
    round = round+1
    
    
    # Cells wih NA as cluster have to be clustered
    to_recluster_old <- rownames(df.cell.clustering)[is.na(df.cell.clustering$cl_1)]
    
    print(paste("Length array cells to re-cluster:",length(to_recluster_new),"round:",round, sep = " "))
    
    #Clustering using Seurat
    seurat.obj <- CreateSeuratObject(counts = (obj@raw[,colnames(obj@raw) %in% to_recluster_new]), 
                                   project = "reclustering", min.cells = 1, min.features = 2)
    seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", 
                                scale.factor = 10000)
    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", 
                                     nfeatures = 2000)
    all.genes <- rownames(seurat.obj)
    seurat.obj <- ScaleData(seurat.obj, features = all.genes)
    seurat.obj <- RunPCA(seurat.obj,npcs = min(c(50,(nrow(seurat.obj@meta.data)-1))), 
                         features = VariableFeatures(object = seurat.obj))
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:min(c(25,(nrow(seurat.obj@meta.data)-1))))
    seurat.obj <- FindClusters(seurat.obj, resolution = 0.5,algorithm = 2)
    cl.resolution <- 0.5
    ##################
    
    # The next lines are necessary to make cluster smaller while the number of residual cells decrease and to
    # stop clustering if the algorithm gives too many singletons.
    while(number.cls.old >= length(unique(seurat.obj$seurat_clusters)) & cl.resolution <= 2){
      cl.resolution <- cl.resolution + 0.5
      print(paste0("Same or lower cluster number: reclustering at higher resolution - resolution ",cl.resolution))
      seurat.obj <- FindClusters(seurat.obj, resolution = cl.resolution,algorithm = 2)
      }
    
    
   # if(length(unique(seurat.obj@meta.data$seurat_clusters)) <= 10 & dim(seurat.obj@meta.data)[1] > 50){
  #    print("We try to cluster with 2 of resolution")
   #   seurat.obj.tmp1 <- FindClusters(seurat.obj, resolution = 2, algorithm = 2,group.singletons = F)
      
  #    sum_cl <- t(summary(seurat.obj.tmp1@meta.data$seurat_clusters))
   #   if(!"singleton" %in% colnames(sum_cl)){
  #      print("No singletons (1) !")
   #     seurat.obj <- seurat.obj.tmp1
  #    }else if( (sum_cl[1,"singleton"]< sum(sum_cl[1,]/2) )) {
   #     seurat.obj <- seurat.obj.tmp1
  #    }
    #}
    
    seurat.obj <- RunUMAP(seurat.obj,umap.method = "uwot",metric = "cosine", 
                          dims = 1:min(c(25,(nrow(seurat.obj@meta.data)-1))))
    
    out_dir <- paste(out_dir_root,"round.",round,"/", sep = "")
    
    if(!file.exists(out_dir)){
      dir.create(file.path(out_dir))
    #  dir.create(file.path(paste(out_dir,"cleaning",sep = "" )))
    }
    print(paste0("PDF UMAP! ", out_dir, "pdf_umap.pdf"))
    pdf(paste0(out_dir,"pdf_umap.pdf"))
      plot(DimPlot(seurat.obj, reduction = "umap",label = TRUE)+
             annotate(geom="text", x=0, y=30, 
                      label=paste0("Cell number: ",length(to_recluster_old),"\nCl. resolution: ",
                                   cl.resolution),color="black"))
    dev.off()
    
    #Next lines check for each new cell cluster if it is homogeneous,
    #in homog.clusters are saved the homogeneous clusters while in 
    #to_reclusters_new are kept the other cells.

    #to_recluster = cluster_homogeneity(seurat.obj,data.raw,out_dir,cond = cond,cores = cores)
    to_recluster_new <- NA
    for (cl in unique(seurat.obj@meta.data$seurat_clusters)[!unique(seurat.obj@meta.data$seurat_clusters) == "singleton"]) {
      cells <- rownames(seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters == cl,])
      if (length(cells) < 10 ) {
        to_recluster_new <- c(to_recluster_new,cells)
      }else{
        to_recluster_new <- c(to_recluster_new,cluster_homogeneity_check(obj = obj,cells = cells,
                                                    out_dir = out_dir,
                                                    cores = cores,
                                                    code = cl) 
        )
      }
    }
    homog.clusters <- to_recluster_new[stringr::str_detect(to_recluster_new,pattern = "Cluster")]
    number.cls.old <- length(unique(seurat.obj$seurat_clusters))- (length(homog.clusters)-1)
    to_recluster_new <- to_recluster_new[stringr::str_detect(to_recluster_new,pattern = "Cluster",negate = T)]
    to_recluster_new <- to_recluster_new[!is.na(to_recluster_new)]
    
    #clusters.to.recluster <- unique(seurat.obj@meta.data[rownames(seurat.obj@meta.data) %in% to_recluster_new,]$seurat_clusters)
    
    if (length(to_recluster_new) == length(to_recluster_old) & cl.resolution < 2.5) {
        next  
    }
    
    
    #Next line store the number of cluster that are ok to be saved.
    #The problem is: we need to change to not duplicate in the obj@cluster array!
    cl_number_to_store = unique(seurat.obj@meta.data[!rownames(seurat.obj@meta.data) %in% to_recluster_new,]$seurat_clusters)
    
    #cluster number not already used
    cl_number_to_store_ready = cl_number_to_store[!cl_number_to_store %in% df.cell.clustering$cl_1]
    cl_number_to_store_ready = as.vector(cl_number_to_store_ready)
    
    # singletons are saved as not clustered so need to be added to the dim(seurat.obj@meta.data)[1] in the next round
    cl_cells_ready =  seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters %in% cl_number_to_store_ready,]
    
    singletons = singletons+dim(cl_cells_ready[cl_cells_ready$seurat_clusters == "singleton",])[1]
    
    #obj@clusters[rownames(cl_cells_ready)] = as.numeric(as.vector(cl_cells_ready$seurat_clusters))
    df.cell.clustering[rownames(cl_cells_ready),"cl_1"] <- as.numeric(as.vector(cl_cells_ready$seurat_clusters))
    
    
    #cluster to change value because already present
    cl_to_change = as.vector(cl_number_to_store[!cl_number_to_store %in% cl_number_to_store_ready])
    cl_cells_to_change = seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters %in% cl_to_change,]
    #cl_not_used = c(min(unique(obj@clusters),na.rm = T):max(unique(obj@clusters),na.rm = T))[!c(min(unique(obj@clusters),na.rm = T):max(unique(obj@clusters),na.rm = T)) %in% unique(obj@clusters)]
    cl_not_used = c(min(unique(df.cell.clustering$cl_1),na.rm = T):
                      max(unique(df.cell.clustering$cl_1),na.rm = T))[!c(min(df.cell.clustering$cl_1,na.rm = T):
                                                                           max(unique(df.cell.clustering$cl_1),na.rm = T)) %in% unique(df.cell.clustering$cl_1)]
    k=1
    while (!is.na(cl_to_change[k])) {
      cells_cluster = rownames(seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters == cl_to_change[k],])
      if ( !is.na(cl_not_used[k])) {
        #obj@clusters[cells_cluster] = cl_not_used[k]
        df.cell.clustering[cells_cluster,"cl_1"] <- cl_not_used[k]
      }else{
        #obj@clusters[cells_cluster] = (max(unique(as.numeric(as.vector(obj@clusters))),na.rm = T)+1)
        df.cell.clustering[cells_cluster,"cl_1"] <- (max(unique(as.numeric(as.vector(df.cell.clustering$cl_1))),na.rm = T)+1)
        
      }
      k=k+1
    }
    
    if (!sum(is.na(df.cell.clustering$cl_1)) == length(to_recluster_new)) {
      print("Some problems in cells reclustering")
      break
    }
    
    if (length(to_recluster_new) <= 50 | cl.resolution >= 2.5) {
      print(paste0("NO new possible clusters! Cell left: ",(length(to_recluster_new)-1), " They will be removed!"))
      obj@raw = obj@raw[,!colnames(obj@raw) %in% to_recluster_new]
      #obj@clusters = obj@clusters[!names(obj@clusters) %in% to_recluster_new]
      df.cell.clustering$names <- rownames(df.cell.clustering)
      df.cell.clustering <- df.cell.clustering[!rownames(df.cell.clustering) %in% to_recluster_new,]
      if ((!all(colnames(obj@raw) %in% rownames(df.cell.clustering)) & 
           !all(rownames(df.cell.clustering) %in% colnames(obj@raw)))) {
        print("Problem: different cells in raw and clusters!")
      }
      saveRDS(df.cell.clustering,paste(out_dir_root,"df.clusters_",cond,".cotan.RDS",sep = ""))
      if (identical(colnames(obj@raw), rownames(df.cell.clustering))) {
        obj@clusters <- df.cell.clustering$cl_1
        names(obj@clusters) <- rownames(df.cell.clustering)
      }else{
        errorCondition("Problems between obj@raw and clustered cells!")
        break
      }
      saveRDS(obj,paste(out_dir_root,"obj_",cond,".cotan.RDS",sep = ""))
    }
    
    #print("Saving left cells.")  
    #print(paste0("Left to recluster: ", (length(to_recluster_new)-1), " cells."))
    #write.csv(to_recluster_new,paste(out_dir_root,"to_recluster_",round,"tot.csv",sep = ""),row.names = F)
    #saveRDS(obj,paste(out_dir_root,"obj_",cond,".cotan.RDS",sep = ""))
    
  } # End while
  
  if (any(is.na(df.cell.clustering$cl_1))) {
    print("Problems! Some NA left!")
    #break
  }
  
  obj@meta <- rbind(obj@meta,c("n. cells left out by clustering:", (length(to_recluster_new)-1)))
  
  if ( (obj@n_cells - length(to_recluster_new)) !=  dim(obj@raw)[2]) {
    print("Problems with the cell number! Check!")
    #break
  }
  
  obj@n_cells <- dim(obj@raw)[2]
  
  #--------------------------
  cotan <- obj@clusters
  cotan <- as.data.frame(cotan)
  meta.data <- merge(srat@meta.data,cotan,by = "row.names",all.x = T)
  rownames(meta.data) <- meta.data$Row.names
  meta.data <- meta.data[rownames(srat@meta.data),]
  if(!identical(rownames(meta.data),rownames(srat@meta.data))){
    print("Problems in the cell codes.")
    break
  }
  srat@meta.data <- meta.data
  
  srat <- subset(srat, subset = cotan >= 0 )

  
  if (!dim(srat)[2] == dim(obj@raw)[2]) {
    print("Problem! Dimensions!")
    break
  }
  
  print("Cluster, UMAP and Saving the Seurat dataset")
  
  srat <- FindNeighbors(srat, dims = 1:25)
  srat <- FindClusters(srat, resolution = 0.5,algorithm = 2)
  srat <- RunUMAP(srat,umap.method = "uwot",metric = "cosine", 
                        dims = 1:25)
  
  
  saveRDS(srat,paste(out_dir_root,"Seurat_obj_",cond,"_with_cotan_clusters.RDS",sep = ""))
  rm(srat)
  gc()
  #Re estimation of all parameters after cells removing
  if((length(to_recluster_new)) > 0){
    print("Estimate again parameters after cells dropping")
    ttm <- clean(obj)
    print("Claening step done")
    obj <- ttm$object
    rm(ttm)
    gc()
    obj <- cotan_analysis(obj,cores = cores)
    print("Analysis step done")
    gc()
    saveRDS(obj,paste(out_dir_root,"obj_",cond,".cotan.RDS",sep = ""))
    print("Coex estimation last step started")
    obj <- get.coex(obj)
    print("Coex estimation last step done")
    gc()
  }
  
  #saveRDS(obj,paste(out_dir_root,"obj_",cond,".cotan.RDS",sep = ""))
  return(obj)
  
 }
)
