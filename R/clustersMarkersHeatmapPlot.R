#' clustersMarkersHeatmapPlot
#'
#' @param objCOTAN 
#' @param clName 
#' @param markers a list of gene markers
#' @param k number of clusters/color for the dendrogram (just for coloring)
#' @param conditionsList a list of dataframes coming from the 
#' clustersSummaryPlot() function 
#'
#' @import dendextend
#' @import ComplexHeatmap
#' @return
#' @export
#'
#' @examples
clustersMarkersHeatmapPlot <- function(objCOTAN,clName = NULL, markers,k=3,
                               conditionsList = NULL){
  
  if (!is.null(conditionsList)) {
    cond1 <- conditionsList[1] # a data frame coming from clustersSummaryPlot()
    
    if (is.numeric(cond1$Cluster)) {
      cond1$Cluster <- sprintf("%02d",cond1$Cluster)
    }
    
    cond1 <- cond1[,c("Cluster","cond1","CellNumber")] %>%
      pivot_wider(names_from = cond1, values_from = CellNumber)
    
    cond1[,2:3] <- cond1[,2:3]/rowSums(cond1[,2:3])
    cond1 <- as.data.frame(cond1)
    rownames(cond1) <- cond1$Cluster
    cond1 <- cond1[rownames(df.t),]
    cond1[is.na(cond1)] <- 0
    
    
    cond1_col <- c("F"="deeppink1","M"= "darkturquoise") # qui invece di F ed M ci vogliono le due o più condizioni possibili
    hb = rowAnnotation(cond1 = anno_barplot(cond1[,2:3],width = unit(3, "cm"), 
                                      gp = gpar(fill = cond1_col, col = "black"),
                                            align_to = "right",
                                      labels_gp = gpar(fontsize = 12)),
                       annotation_name_rot = 0)
    
  }
  
  
  if(is.null(clName)){
    clName <- getClusterizationName(objCOTAN)
  }
  expression.cl <- clustersDeltaExpression(objCOTAN,clName = clName)
  score.df <- geneSetEnrichment(groupMarkers = markers,
                                clustersCoex = expression.cl)
  
  df.t <- t(score.df[,1:(ncol(score.df)-2)])
  dend <- clustersTreePlot(objCOTAN, k=k)
  dend <- dend$dend
  
  dend = set(dend = dend, "branches_lwd", 2)
  
  rownames(df.t) <- sprintf("%02d",as.numeric(rownames(df.t)))
  
  cls.info <- clustersSummaryPlot(objCOTAN,clName = "merge_round")
  cls.info$data$Cluster <- sprintf("%02d",cls.info$data$Cluster)
  
  rownames(cls.info$data) <- cls.info$data$Cluster
  cls.info$data <- cls.info$data[rownames(df.t),]
  Freq <- cls.info$data$CellNumber
  names(Freq) <- cls.info$data$Cluster
  
  Freq2 <- paste0(round(cls.info$data$CellPercentage,digits = 1),"%")
  names(Freq2) <- cls.info$data$Cluster
  
  
  ha = rowAnnotation(cell.number = anno_numeric(Freq,
                                                bg_gp = gpar(fill = "orange", 
                                                             col = "black"),
                                                #labels_gp = gpar(fontsize = 10)
  ), annotation_name_rot = 0)
  
  ha2 = rowAnnotation(cell.number = anno_text(Freq2,
                                              gp = gpar(fontsize = 10)), 
                      annotation_name_rot = 0)
  
 
  #condition_col <- c("F"="red4","R"="sienna3","H"="seagreen")
  #hc = rowAnnotation(condition = anno_barplot(condition[,2:3],width =unit(3, "cm"), gp = gpar(fill = condition_col, col = "black"),align_to = "right",labels_gp = gpar(fontsize = 12)),annotation_name_rot = 0)
  
  
  lgd_list = list(
    Legend(labels = c("Female", "Male"), title = "cond1", 
          legend_gp = gpar(fill = cond1_col))#,
  #  Legend(labels = c("Flare", "Remission","Healthy"), title = "Condition", 
  #        legend_gp = gpar(fill = condition_col))
  )
  
  col_fun = colorRamp2(c(0, 1), c( "lightblue", "red"))
  

  final.heatmap <- Heatmap(df.t, rect_gp = gpar(col = "white", lwd = 1),
                           cluster_rows = dend,
                           cluster_columns = FALSE,
                           col = col_fun,
                           width =  unit(28, "cm"),
                           row_dend_width = unit(8, "cm"),
                           #height = unit(6, "cm"),
                           column_names_gp = gpar(fontsize = 11),
                           row_names_gp = gpar(fontsize = 11),
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(sprintf("%.1f", df.t[i, j]), x, y, 
                                       gp = gpar(fontsize = 9))},
                           right_annotation = c(ha,ha2)
                                                   left_annotation = c(hb,hc)#,
                           #                       column_title = 
                                                   #paste0("Gene marker set expression in ", sample)
  )
  
  
  final.heatmap <- draw(final.heatmap, annotation_legend_list = lgd_list)
  
  return(list("heatmapPlot"=final.heatmap,"dataScore"= score.df))
  
  }

