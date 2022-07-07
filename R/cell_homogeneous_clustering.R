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
#' @param dataset_type Type of data set in input: "Seurat" if it is a Seurat data set saved as RDS,"COTAN" 
#' if it is a COTAN object, "DF"if it is a data frame saved as csv. 
#' @param GEO GEO or other dataset official code
#' @param sc.method single cell method used fot the experiment 
#' @import Seurat
#' @return the scCOTAN object and it saves, in the out_dir, the Seurat elaborated data file
#' @export 
#'
#' @examples
setGeneric("cell_homogeneous_clustering", function(cond,out_dir,in_dir,cores=1, dataset_name,dataset_type,
                                                   GEO, sc.method )
  standardGeneric("cell_homogeneous_clustering"))
#' @rdname cell_homogeneous_clustering
setMethod("cell_homogeneous_clustering","Seurat",
 function(cond,out_dir,in_dir,cores, dataset_name,dataset_type){
  
  out_dir_root <- paste0(out_dir,cond,"/")
  if(!file.exists(out_dir_root)){
    dir.create(file.path(out_dir_root))
    
  }
  dir.create(file.path(paste(out_dir_root,"cleaning",sep = "" )))
  
  #To load the seurat object crated
  
  if(dataset_type == "Seurat"){
    srat <- readRDS(paste0(in_dir,dataset_name))
    
    data.raw <- srat@assays$RNA@counts
    
    obj <- new("scCOTAN",raw = data.raw[,rownames(srat@meta.data)] )
    obj = initRaw(obj,GEO = GEO ,sc.method= sc.method,cond = cond)
    
  }else if(dataset_type == "COTAN"){
    obj <- readRDS(paste0(in_dir,dataset_name))
    
    srat <- CreateSeuratObject(counts = as.data.frame(obj@raw), project = cond, 
                               min.cells = 3, min.features = 200)
    srat <- NormalizeData(srat)
    srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(srat)
    srat <- ScaleData(srat, features = all.genes)
    srat <- RunPCA(srat, features = VariableFeatures(object = srat))
      
  }else if(dataset_type == "DF"){
    DF <- read.csv(paste0(in_dir,dataset_name))
    
    srat <- CreateSeuratObject(counts = as.data.frame(DF), project = cond, 
                               min.cells = 3, min.features = 200)
    srat <- NormalizeData(srat)
    srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(srat)
    srat <- ScaleData(srat, features = all.genes)
    srat <- RunPCA(srat, features = VariableFeatures(object = srat))
    
    obj <- new("scCOTAN",raw = as.data.frame(DF) )
    obj = initRaw(obj,GEO = GEO ,sc.method= sc.method,cond = cond)
  }
  
  out_dir <- out_dir_root
  
  srat <- FindNeighbors(srat, dims = 1:25)
  srat <- FindClusters(srat, resolution = 1)
  
  gc()
  to_recluster <- cluster_homogeneity(srat,data.raw, out_dir_root, cond = cond,cores = cores)
  ##########
  
  #tmp_meta = srat@meta.data[names(obj@clusters[!names(obj@clusters) %in% to_recluster]),]
  tmp_meta <- srat@meta.data[!rownames(srat@meta.data) %in% to_recluster,]
  
  # To save the already defined clusters
  obj@clusters[rownames(tmp_meta)] <- as.vector(tmp_meta$seurat_clusters)
  
  saveRDS(obj,paste(out_dir_root,"obj_",cond,".cotan.RDS",sep = ""))
  write.csv(to_recluster,paste(out_dir_root,"to_recluster_",cond,"tot.csv",sep = ""),row.names = F)
  
  print("First part finished! Starting the iterative part.")
  
  round <- 0
  singletons <- 0
  n_cells <- length(obj@clusters)
  to_recluster_old <- NA
  while (length(to_recluster) > 25 & (length(to_recluster_old) != length(to_recluster)) ) {
    #data.raw = Read10X(get(paste0("raw.",cond)))
    gc()
    round = round+1
    #to check
    
    print(paste("Length array cells to re-cluster:",length(to_recluster[!is.na(to_recluster)]),"round:",round, sep = " "))
    #data.raw = as.data.frame(obj@raw)
    seurat.obj <- CreateSeuratObject(counts = (data.raw[,colnames(data.raw) %in% to_recluster]), project = "reclustering", min.cells = 3, min.features = 200)
    #print(dim(seurat.obj))
    seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    
    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seurat.obj)
    seurat.obj <- ScaleData(seurat.obj, features = all.genes)
    
    seurat.obj <- RunPCA(seurat.obj,npcs = min(c(50,(nrow(seurat.obj@meta.data)-1))), features = VariableFeatures(object = seurat.obj))
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:min(c(50,(nrow(seurat.obj@meta.data)-1))))
    
    seurat.obj <- FindClusters(seurat.obj, resolution = 1.1,algorithm = 2)
    
    # The next lines are necessary to make cluster smaller while the numeber of residual cells decrease and to
    # stop clustering if the algorithm gives too many singletons.
    if(length(unique(seurat.obj@meta.data$seurat_clusters)) <= 10 & dim(seurat.obj@meta.data)[1] > 50){
      print("We try to cluster with 1.5 of resolution")
      seurat.obj.tmp1 <- FindClusters(seurat.obj, resolution = 1.5, algorithm = 2,group.singletons = F)
      
      sum_cl =t(summary(seurat.obj.tmp1@meta.data$seurat_clusters))
      if(!"singleton" %in% colnames(sum_cl)){
        print("No singletons (1) !")
        seurat.obj =seurat.obj.tmp1
      }else if( (sum_cl[1,"singleton"]< sum(sum_cl[1,]/2) )) {
        seurat.obj =seurat.obj.tmp1
      }
      
    }
    
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:min(c(50,(nrow(seurat.obj@meta.data)-1))))
    
    out_dir =paste(out_dir_root,"round.",round,"/", sep = "")
    
    if(!file.exists(out_dir)){
      dir.create(file.path(out_dir))
      dir.create(file.path(paste(out_dir,"cleaning",sep = "" )))
    }
    pdf(paste0(out_dir,"pdf_umap.pdf"))
    print(DimPlot(seurat.obj, reduction = "umap",label = TRUE))
    dev.off()
    #obj2@clusters =
    to_recluster_old = to_recluster
    to_recluster = cluster_homogeneity(seurat.obj,data.raw,out_dir,cond = cond,cores = cores)
    
    #The first time length(to_recluster_old)==length(to_recluster) we try to make the cluster smaller
    if (length(to_recluster_old)==length(to_recluster)) {
      print("We try to cluster with 2.5 of resolution")
      seurat.obj.tmp2 <- FindClusters(seurat.obj, resolution = 3, group.singletons = F)
      sum_cl =t(summary(seurat.obj.tmp2@meta.data$seurat_clusters))
      if(!"singleton" %in% colnames(sum_cl)){
        print("No singletons (2) !")
        seurat.obj =seurat.obj.tmp2
      }else if (sum_cl[1,"singleton"]< sum(sum_cl[1,]/2)) {
        seurat.obj =seurat.obj.tmp2
      }
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:min(c(50,(nrow(seurat.obj@meta.data)-1))))
      
      out_dir =paste(out_dir_root,"round.",round,"/", sep = "")
      
      if(!file.exists(out_dir)){
        dir.create(file.path(out_dir))
        dir.create(file.path(paste(out_dir,"cleaning",sep = "" )))
      }
      pdf(paste(out_dir,"pdf_umap.pdf"))
      print(DimPlot(seurat.obj, reduction = "umap",label = TRUE))
      dev.off()
      #obj2@clusters =
      to_recluster_old = to_recluster
      to_recluster = cluster_homogeneity(seurat.obj,data.raw,out_dir,cond = cond,cores = 11)
    }
    
    if (length(to_recluster_old)==length(to_recluster)) {
      print(paste0("NO new possible clusters! Cell left:",length(to_recluster), " They will be removed!"))
      obj@raw = obj@raw[,!colnames(obj@raw) %in% to_recluster]
      obj@clusters = obj@clusters[names(obj@clusters) %in% colnames(obj@raw)]
      if (any(is.na(obj@clusters))) {
        print("Problems! Some NA left!")
      }
      
      saveRDS(obj,paste(out_dir_root,"obj_",cond,"_in_cotan.RDS",sep = ""))
      write.csv(to_recluster,paste(out_dir_root,"to_recluster_",cond,"tot.csv",sep = ""),row.names = F)
      break
      
    }
    # Next line store the number of cluster that are ok to be saved.
    #The problem is: we need to change to not duplicate in the obj@cluster array!
    cl_number_to_store = unique(seurat.obj@meta.data[!rownames(seurat.obj@meta.data) %in% to_recluster,]$seurat_clusters)
    
    #cluster number not already used
    cl_number_to_store_ready = cl_number_to_store[!cl_number_to_store %in% obj@clusters]
    #cl_number_to_store_ready = as.numeric(as.vector(cl_number_to_store_ready))
    cl_number_to_store_ready = as.vector(cl_number_to_store_ready)
    
    # singletons are saved as not clustered so need to be added to the dim(seurat.obj@meta.data)[1] in the next round
    cl_cells_ready =  seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters %in% cl_number_to_store_ready,]
    
    singletons = singletons+dim(cl_cells_ready[cl_cells_ready$seurat_clusters == "singleton",])[1]
    
    obj@clusters[rownames(cl_cells_ready)] = as.numeric(as.vector(cl_cells_ready$seurat_clusters))
    #cluster to change value because already present
    cl_to_change = as.vector(cl_number_to_store[!cl_number_to_store %in% cl_number_to_store_ready])
    cl_cells_to_change = seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters %in% cl_to_change,]
    cl_not_used = c(min(unique(obj@clusters),na.rm = T):max(unique(obj@clusters),na.rm = T))[!c(min(unique(obj@clusters),na.rm = T):max(unique(obj@clusters),na.rm = T)) %in% unique(obj@clusters)]
    #change_matrix = data.frame("cl_to_change"=cl_to_change,"not_used"=cl_not_used)
    
    
    k=1
    while (!is.na(cl_to_change[k])) {
      cells_cluster = rownames(seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters == cl_to_change[k],])
      if ( !is.na(cl_not_used[k])) {
        obj@clusters[cells_cluster] = cl_not_used[k]
      }else{
        obj@clusters[cells_cluster] = (max(unique(as.numeric(as.vector(obj@clusters))),na.rm = T)+1)
      }
      k=k+1
      
      #out_dir = out_dir_root
    }
    
    
    
    print("Saving the COTAN dataset")  
    saveRDS(obj,paste(out_dir_root,"obj_",cond,"in_cotan.RDS",sep = ""))
    print(paste0("Left to recluster: ", (length(to_recluster)-1), " cells."))
    write.csv(to_recluster,paste(out_dir_root,"to_recluster_",cond,"tot.csv",sep = ""),row.names = F)
    
    #if(length(obj@clusters) != (n_cells + length(to_recluster_old) - length(to_recluster)) ){
    #  print("Problems: changed the number of cells! ")
    #  break
    #}
    
  }
  
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
  
  
  print("Saving the Seurat dataset")
  saveRDS(srat,paste(out_dir_root,"Seurat_obj_",cond,"_with_cotan_clusters.RDS",sep = ""))
  return(obj)
  
 }
)
