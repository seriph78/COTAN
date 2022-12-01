#' cell_homogeneous_clustering
#'
#'This function imports a Seurat object from in_dir with a file name specified by data set_name.
#'From this, and using COTAN functions, it splits the data set until all clusters are homogeneous by GDI
#'(calling the cluster_homogeneity function).
#' @param cond a string specifying the experiment condition
#' @param out_dir an existing directory for the analysis output. In this directory is saved also the final Seurat object.
#' @param cores number of cores (NB for windows system no more that 1 can be used)
#' @import Seurat
#' @importFrom stringr str_detect
#' @return the scCOTAN object and it saves, in the out_dir, the Seurat elaborated data file
#' @export
#'
#' @examples
setGeneric("cell_homogeneous_clustering", function(obj,cond,out_dir,cores=1)
  standardGeneric("cell_homogeneous_clustering"))
#' @rdname cell_homogeneous_clustering
setMethod("cell_homogeneous_clustering",
          "COTAN",
 function(obj,cond,out_dir,cores){

  out_dir_cond <- paste0(out_dir,"/",cond,"/")
  if(!file.exists(out_dir_cond)){
    dir.create(file.path(out_dir_cond))
  }
  #Step 1
  cell.cluster.vector <- rep(NA, length = getNumCells(obj))
  names(cell.cluster.vector) <- getCells(obj)
  
  srat <- CreateSeuratObject(counts = as.data.frame(obj@raw), project = cond,
                             min.cells = 3, min.features = 4)
  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(srat)
  srat <- ScaleData(srat, features = all.genes)
  srat <- RunPCA(srat, features = VariableFeatures(object = srat))

  srat <- FindNeighbors(srat, dims = 1:25)
  srat <- FindClusters(srat, resolution = 0.5,algorithm = 2)

  srat <- RunUMAP(srat,umap.method = "uwot",metric = "cosine",
                        dims = 1:min(c(50,(nrow(srat@meta.data)-1))))
  print(paste0("PDF UMAP! ", out_dir_cond, "pdf_umap.pdf"))
  pdf(paste0(out_dir_cond,"pdf_umap.pdf"))
  plot(DimPlot(srat, reduction = "umap",label = TRUE,group.by = "orig.ident"))
  plot(DimPlot(srat, reduction = "umap",label = TRUE))
  dev.off()
  gc()

  # Step 2
  clustering.number <- length(obj@clustersCoex) + 1
  clustering.name <- paste0("CL_",clustering.number)
  
  #obj@metaCells[clustering.name] <- as.data.frame(matrix(nrow = getNumCells(obj),ncol = 1))
  #rownames(obj@clustersCoex[[clustering.number]]) <- getCells(obj)
  #colnames(obj@clustersCoex[[clustering.number]]) <- "cl_1"

  to_recluster_new <- NA
  for (cl in unique(srat@meta.data$seurat_clusters)) {
  print(cl)
    cells <- rownames(srat@meta.data[srat@meta.data$seurat_clusters == cl,])
    to_recluster_new <- c(to_recluster_new,cluster_homogeneity_check(obj = obj,
                                                                     cells = cells,
                                                                     out_dir = out_dir_cond,
                                                                     cores = cores,
                                                                     code = cl)
    )
  }

  homog.clusters <- to_recluster_new[stringr::str_detect(to_recluster_new,pattern = "Cluster")]
  print(paste(homog.clusters, collapse = " "))
  to_recluster_new <- to_recluster_new[stringr::str_detect(to_recluster_new,pattern = "Cluster",negate = T)]
  to_recluster_new <- to_recluster_new[2:length(to_recluster_new)]

  #clusters.to.recluster <- unique(srat@meta.data[rownames(srat@meta.data) %in% to_recluster_new,]$seurat_clusters)

  tmp_meta <- srat@meta.data[!rownames(srat@meta.data) %in% to_recluster_new,]

  # To save the already defined clusters
  cell.cluster.vector[rownames(tmp_meta)] <- tmp_meta$seurat_clusters
  
  if (!sum(is.na(cell.cluster.vector)) == length(to_recluster_new)) {
    print("Problems!")
  }

  print("First part finished! Starting the iterative part.")
  round <- 0
  singletons <- 0
  to_recluster_old <- names(cell.cluster.vector[is.na(cell.cluster.vector)])
  number.cls.old <- 0 #Just to start from a number higher than any possible real situation

  while (sum(is.na(cell.cluster.vector)) > 50 ) { #& round < 25
    #data.raw = Read10X(get(paste0("raw.",cond)))
    gc()
    round = round+1

    #to_recluster_new <- names(cell.cluster.vector[is.na(cell.cluster.vector)])

    print(paste("Length array cells to re-cluster:", length(to_recluster_old),
                "round:", round))

    #Clustering using Seurat
#    seurat.obj <- CreateSeuratObject(counts = (obj@raw[,getCells(obj) %in% to_recluster_new]),
#                                   project = "reclustering", min.cells = 1, min.features = 2)
    seurat.obj <- CreateSeuratObject(counts = (obj@raw[,getCells(obj) %in% to_recluster_old]),
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

    seurat.obj <- RunUMAP(seurat.obj,umap.method = "uwot",metric = "cosine",
                          dims = 1:min(c(25,(nrow(seurat.obj@meta.data)-1))))

    out_dir_round <- paste0(out_dir_cond, "round.", round, "/")

    if(!file.exists(out_dir_round)){
      dir.create(file.path(out_dir_round))
    #  dir.create(file.path(paste0(out_dir_round, "cleaning")))
    }
    print(paste0("PDF UMAP! ", out_dir_round, "pdf_umap.pdf"))
    pdf(paste0(out_dir_round,"pdf_umap.pdf"))
      plot(DimPlot(seurat.obj, reduction = "umap",label = TRUE)+
             annotate(geom="text", x=0, y=30,
                      label=paste0("Cell number: ",length(to_recluster_old),"\nCl. resolution: ",
                                   cl.resolution),color="black"))
    dev.off()

    #Next lines check for each new cell cluster if it is homogeneous,
    #in homog.clusters are saved the homogeneous clusters while in
    #to_reclusters_new are kept the other cells.

    to_recluster_new <- NA
    for (cl in unique(seurat.obj@meta.data$seurat_clusters)[!unique(seurat.obj@meta.data$seurat_clusters) == "singleton"]) {
      cells <- rownames(seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters == cl,])
      if (length(cells) < 10 ) {
        to_recluster_new <- c(to_recluster_new,cells)
      }else{
        to_recluster_new <- c(to_recluster_new,cluster_homogeneity_check(obj = obj,cells = cells,
                                                    out_dir = out_dir_round,
                                                    cores = cores,
                                                    code = cl)
        )
      }
    }
    homog.clusters <- to_recluster_new[stringr::str_detect(to_recluster_new,pattern = "Cluster")]
    number.cls.old <- length(unique(seurat.obj$seurat_clusters)) - (length(homog.clusters)-1)
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
    cl_number_to_store_ready = cl_number_to_store[!cl_number_to_store %in% cell.cluster.vector]
    cl_number_to_store_ready = as.vector(cl_number_to_store_ready)

    cl_cells_ready =  seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters %in% 
                                             cl_number_to_store_ready,]
    
    # singletons are saved as not clustered so need to be added to the 
    singletons = singletons+dim(cl_cells_ready[cl_cells_ready$seurat_clusters == "singleton",])[1]

    #obj@clusters[rownames(cl_cells_ready)] = as.numeric(as.vector(cl_cells_ready$seurat_clusters))
    cell.cluster.vector[rownames(cl_cells_ready)] <- as.numeric(as.vector(cl_cells_ready$seurat_clusters))

    #cluster to change value because already present
    cl_to_change = as.vector(cl_number_to_store[!cl_number_to_store %in% cl_number_to_store_ready])
    cl_cells_to_change = seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters %in% cl_to_change,]
    #cl_not_used = c(min(unique(obj@clusters),na.rm = T):max(unique(obj@clusters),na.rm = T))[!c(min(unique(obj@clusters),na.rm = T):max(unique(obj@clusters),na.rm = T)) %in% unique(obj@clusters)]
    cl_not_used = which(!c(1:max(cell.cluster.vector, na.rm = T)) %in% cell.cluster.vector)
    
    k=1
    for (cl.used in cl_to_change) {
      cells_cluster = rownames(seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters == cl.used,])
      if (k <= length(cl_not_used)) {
        cell.cluster.vector[cells_cluster] <- cl_not_used[k]
        k <- k+1
      }else{
        cell.cluster.vector[cells_cluster] <- max(cell.cluster.vector,na.rm = T)+1
      }
    }
    
    
    if (!sum(is.na(cell.cluster.vector)) == length(to_recluster_new)) {
      print("Some problems in cells reclustering")
      break
    }

    if (length(to_recluster_new) <= 50 | cl.resolution >= 2.5) {
      print(paste0("NO new possible clusters! Cell left: ",length(to_recluster_new)))
      #obj@raw = obj@raw[,!getCells(obj) %in% to_recluster_new]
      #obj@clusters = obj@clusters[!names(obj@clusters) %in% to_recluster_new]
    #  df.cell.clustering$names <- rownames(df.cell.clustering)
      #obj <- dropGenesCells(obj,cells = to_recluster_new)
      #df.cell.clustering <- df.cell.clustering[!rownames(df.cell.clustering) %in% to_recluster_new,]
      #if ((!all(getCells(obj) %in% rownames(df.cell.clustering)) &
       #    !all(rownames(df.cell.clustering) %in% getCells(obj)))) {
      #  print("Problem: different cells in raw and clusters!")
      #}
      #saveRDS(df.cell.clustering,
       #       paste0(out_dir_round, "df.clusters_", cond, ".cotan.RDS"))
      #if (identical(getCells(obj), rownames(df.cell.clustering))) {
      #  obj@clusters <- df.cell.clustering$cl_1
      #  names(obj@clusters) <- rownames(df.cell.clustering)
      #}else{
      #  errorCondition("Problems between obj@raw and clustered cells!")
      #  break
      #}
      #saveRDS(obj, paste0(out_dir, "obj_", cond, ".cotan.RDS"))
    }

  } # End while

  if (!sum(is.na(cell.cluster.vector)) == length(to_recluster_new)) {
    print("Problems! Some NA left!")
    #break
  }

  obj <- addElementToMetaDataset(obj,tag = "n. cells left out by clustering:",length(to_recluster_new))

  #--------------------------
# Da mettere in una funzione esterna!
  #cotan <- as.data.frame(cotan)
  meta.data <- merge(srat@meta.data,
                     as.matrix(cell.cluster.vector),
                     by = "row.names",all.x = T)
  colnames(meta.data)[dim(meta.data)[2]] <- paste0("cotan_",clustering.name)
  rownames(meta.data) <- meta.data$Row.names
  meta.data <- meta.data[rownames(srat@meta.data),]
  if(!identical(rownames(meta.data),rownames(srat@meta.data))){
    print("Problems in the cell codes.")
    break
  }
  srat@meta.data <- meta.data

  #srat <- subset(srat, subset = cotan >= 0 )


  if (!dim(srat)[2] == dim(obj@raw)[2]) {
    print("Problem! Dimensions!")
    break
  }

  print("Cluster, UMAP and Saving the Seurat dataset")

  srat <- FindNeighbors(srat, dims = 1:25)
  srat <- FindClusters(srat, resolution = 0.5,algorithm = 2)
  srat <- RunUMAP(srat,umap.method = "uwot",metric = "cosine",
                        dims = 1:25)


  saveRDS(srat,
          paste0(out_dir, "Seurat_obj_", cond, "_with_cotan_clusters.RDS"))
  rm(srat)
  #gc()
  #Re estimation of all parameters after cells removing


  #if ((length(to_recluster_new)) > 0) {
  #  print("Estimate again parameters after cells dropping")
  #  obj <- as(clean(obj, calcExtraData = FALSE)[["objCOTAN"]], "scCOTAN")
  #  print("Cleaning step done")

  #  saveRDS(obj, paste0(out_dir, "obj_", cond, ".cotan.RDS"))
  #}

  gc()

  #obj <- as(estimateDispersion(obj, cores = cores), "scCOTAN")
  #print("Analysis step done")
  #gc()
  #saveRDS(obj, paste0(out_dir, "obj_", cond, ".cotan.RDS"))
  #gc()
  #print("Coex estimation last step started")
  #obj <- calculateCoex(obj)
  #print("Coex estimation last step done")
  #gc()

  #obj <- as(obj, "scCOTAN")
  #obj <- addClusterization(obj,clusters = cell.cluster.vector, clusterizationName = clustering.name)
  
  #saveRDS(obj, paste0(out_dir, "obj_", cond, ".cotan.RDS"))
  return(cell.cluster.vector)

 }
)
