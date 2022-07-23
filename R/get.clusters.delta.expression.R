
#' get.clusters.delta.expression
#'
#'This function estimate the change in gene expression inside the cluster for each genes
#'compared to the average situation in the data set
#'
#' @param obj a COTAN object with defined clusters (field obj@clusters filled)
#'
#' @return a dataframe 
#' @export
#'
#' @examples
setGeneric("get.clusters.delta.expression", function(obj)
  standardGeneric("get.clusters.delta.expression"))
#' @rdname get.clusters.delta.expression
setMethod("get.clusters.delta.expression","scCOTAN",
 function(obj){
    clusters.names = unique(obj@clusters)[!is.na(unique(obj@clusters))]
    list.clusters = list(names(obj@clusters[obj@clusters %in% clusters.names[1]]))
    names(list.clusters)=clusters.names[1]
    for (c in c(2:length(clusters.names))) {
      tmp = list(names(obj@clusters[obj@clusters %in% clusters.names[c]]))
      names(tmp)= clusters.names[c]
      list.clusters = c(list.clusters,tmp)
    }
    
    increased.expression.tot <- data.frame()
    
    cells <- obj@raw
    mu_estimator <- mu_est(obj)
    
    hk <- obj@hk
    
    mu_estimator <- mu_estimator[!rownames(mu_estimator) %in% hk,]
    cells <- cells[!rownames(cells) %in% hk, ]
    M = fun_pzero(obj@a,mu_estimator[,colnames(cells)]) # matrix of 0 probabilities
    N = 1-M
    
    for (condition in names(list.clusters)) {
      gc()
      print(paste("cluster",condition, sep=" "))
      
      cells_set <-  unlist(list.clusters[condition])
      
      stopifnot("ERROR. Some cells are not present!"=  all(cells_set %in% colnames(cells)) )
      
      
      yes_in <- rowSums(cells[,colnames(cells) %in% cells_set] > 0)
      
      estimator_yes_in = rowSums(N[,colnames(cells) %in% cells_set])
      
      sum.UDE.cluster <- sum(obj@nu[cells_set])
      
      increse.expression <- (yes_in - estimator_yes_in)/sum.UDE.cluster
      
      increse.expression <- as.data.frame(increse.expression)
      
      colnames(increse.expression) <- paste("cl.",condition,sep = "")
      
      if (dim(increased.expression.tot)[1] == 0) {
        increased.expression.tot <- as.data.frame(matrix(nrow = length(obj@coex$genes)))
        rownames(increased.expression.tot) <- obj@coex$genes
      }
      
      increased.expression.tot  <- cbind(increased.expression.tot ,increse.expression )
      
    }
    
    increased.expression.tot<- round(increased.expression.tot, digits = 3)
    
    return(increased.expression.tot)
    
  }
)