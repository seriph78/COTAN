#Conditions is a named list containing the condition name and in each arrays, the cell regarding that condition
# default null: just one condition

#' clusters.summary.plot
#'
#' @param objCOTAN a `COTAN` object
#' @param condition name of the column in the MetadataCells dataframe containing 
#' the condition. Default NULL
#' @param cl.name name of the clustering to consider (without the initial "CL_")
#'
#' @return a list of a dataframe with the data and a ggplot object
#' @export
#'
#' @examples
clusters.summary.plot <- function(objCOTAN, condition = NULL, cl.name = "merge_round"){
  MetadataCells <- getMetadataCells(objCOTAN)
  last.clusterization <- paste0("CL_",cl.name)
  
  
  df <- as.data.frame(table(MetadataCells[,c(last.clusterization,condition)]))
  #rownames(df) <- df$Var1
  if (dim(df)[2] == 2) {
    no.cond <- TRUE
    colnames(df) <- c("Cluster","CellNumber")
  }else if(dim(df)[2] == 3){
    no.cond <- FALSE
    colnames(df) <- c("Cluster",condition,"CellNumber")
  }else{
    errorCondition("Problem!")
  }
  
  
  df$MeanUDE <- NA
  df$MedianUDE <- NA
  df$ExpGenes25 <- NA
  df$ExpGenes <- NA
  
  #UDE <- setNames(rep(NA,length(Freq)), as.character(unique(MetadataCells[,last.clusterization])))
  
  #expressed.genes <- setNames(rep(NA,length(Freq)), as.character(unique(MetadataCells[,last.clusterization])))
  
  if (no.cond == TRUE) {
    for (cl in unique(MetadataCells[,last.clusterization])) {
      df[df$Cluster == cl,]$MeanUDE <- mean(MetadataCells[MetadataCells[last.clusterization] == cl,]$nu)
      df[df$Cluster == cl,]$MedianUDE <- median(MetadataCells[MetadataCells[last.clusterization] == cl,]$nu)
      
      cl.dim <- dim(MetadataCells[MetadataCells[last.clusterization] == as.character(cl),])[1]
      cells.array <- rownames(MetadataCells[MetadataCells[last.clusterization] == as.character(cl),])
      
      
      df[df$Cluster == cl,]$ExpGenes25 <- sum(rowSums(getZeroOneProj(objCOTAN)[,cells.array]) > (cl.dim/100)*25)
      df[df$Cluster == cl,]$ExpGenes <- sum(rowSums(getZeroOneProj(objCOTAN)[,cells.array]) > 0)
      
      
    }
      
  }else{
    for (cond in unique(df[,condition])) {
      for (cl in unique(MetadataCells[,last.clusterization])) {
        df[df$Cluster == cl & df[,condition] == cond,]$MeanUDE <- 
          mean(MetadataCells[MetadataCells[,last.clusterization] == cl & MetadataCells[,condition] == cond,]$nu)
        df[df$Cluster == cl & df[,condition] == cond,]$MedianUDE <- 
          median(MetadataCells[MetadataCells[,last.clusterization] == cl & MetadataCells[,condition] == cond,]$nu)
        
        cl.dim <- dim(MetadataCells[MetadataCells[last.clusterization] == as.character(cl) & MetadataCells[,condition] == cond,])[1]
        cells.array <- rownames(MetadataCells[MetadataCells[last.clusterization] == as.character(cl) & MetadataCells[,condition] == cond,])
        
        
        df[df$Cluster == cl & df[,condition] == cond,]$ExpGenes25 <- sum(rowSums(getZeroOneProj(objCOTAN)[,cells.array]) > (cl.dim/100)*25)
        df[df$Cluster == cl & df[,condition] == cond,]$ExpGenes <- sum(rowSums(getZeroOneProj(objCOTAN)[,cells.array]) > 0)
        
        
      }
    }
    
    
  }
  
  df$MeanUDE <-   round(df$MeanUDE,digits = 2)
  df$MedianUDE <-   round(df$MedianUDE,digits = 2)
  
  df$CellPercentage <- df$CellNumber/sum(df$CellNumber)*100
  df$CellPercentage <-   round(df$CellPercentage,digits = 2)
  
  Df <- df %>% 
    gather(keys, values, CellNumber:CellPercentage)
  
  
  statistics.plot <- ggplot(Df, aes(Cluster,values, fill=cond)) +
    geom_bar(position="dodge", stat="identity", color="black") +
    geom_text(aes(x=Cluster, y=values, label=values, vjust=0.5, hjust=-0.1), position = position_dodge(width=1))+
    facet_wrap(~ keys,ncol = 6,scales = "free") +
    scale_y_continuous(expand = expansion(mult = c(.05, .4)))+
    coord_flip()+
    theme_classic()+
    scale_fill_brewer(palette = "Accent")
  
  results <- list("Data"=df,"Plot"=statistics.plot)
  return(results)  
}






#' clusters.tree.plot
#'
#' @param objCOTAN 
#' @param k numeric scalar (OR a vector) with the number of clusters the tree should be cut into.
#'
#' @return the dendrogram
#' @export
#' @import RColorBrewer
#'
#' @examples
clusters.tree.plot <- function(objCOTAN, k){
  cluster_data <- getClusterizationData(objCOTAN)[["coex"]]
  
  ######## This is the best: cosine dissimilarity
  Matrix <- as.matrix(t(cluster_data))
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  tree <- hclust(D_sim,method = "ward.D2")
  
  ############
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  ########
  dend <- as.dendrogram(tree)
  #colnames(df) <- str_remove_all(colnames(df), pattern = "cl.")
  cut = cutree(tree, k = k)
  dend =branches_color(dend,k=k,col=col_vector[1:k],groupLabels = T)
  dend =color_labels(dend,k=k) 
  dend = set(dend = dend, "branches_lwd", 4)
  
  
return(dend)
  
}

