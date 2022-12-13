
#' geneSetEnrichment
#' 
#' Function to get a cumulative score of enrichment in a cluster over a gene set
#'
#' @param objCOTAN a COTAN objCOTANect
#' @param genes a named list of genes
#' @param expression.cl the dataframe for the increased or decreased expression 
#' of every gene in the cluster compared to the whole background   
#'
#' @return a dataframe
#' @export
#'
#' @examples
setMethod(
  "geneSetEnrichment",
  "COTAN",
   function(objCOTAN, expression.cl, genes) {
    
    df <- as.data.frame(matrix(nrow = 1,ncol = ncol(getClusterizationData(objCOTAN)[[1]])+2))
    rownames(df) <- names(genes)
    colnames(df) <- c(colnames(getClusterizationData(objCOTAN)[[1]]),"N. detected","N. total")
    teta <- -1/0.1 * (log(0.25))
    #not_ass_clusters <- NA
    for (m in names(genes)) {
      n.genes <- sum(getGenes(objCOTAN) %in% genes[[m]])
      print(paste0("In ",m , "there are ",n.genes, " detected over ",length(genes[[m]]), " genes"))
      df[m,"N. detected"] <- n.genes
      df[m,"N. total"] <- length(genes[[m]])
      for (ro in colnames(getClusterizationData(objCOTAN)[[1]])) {
        #pv <- p_value[unlist(genes[[m]]),ro]
        #co <- objCOTAN@cluster_data[unlist(genes[[m]]),ro]
        ex <- expression.cl[rownames(expression.cl) %in% genes[[m]],ro]
        ex[ex < 0 & !is.na(ex)] <- 0
        
        ex <- 1-exp(- teta * ex)
        df[m,ro] <- sum(ex,na.rm = T)/n.genes
      }
    }
    
    df <- round(df,digits = 1)
    return(df)
  }

)
