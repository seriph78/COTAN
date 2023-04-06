
#' FindClustersMarkers
#'
#' Function that takes an `COTAN object, a list of markers and produce a
#' `data.frame` with the `n` most positively enriched and the `n` most
#' negatively enriched markers for cluster and tell whether they are in the
#' markers list or not. It also gives the `pValue` and the adjusted p-values
#' (defaults to "bonferroni")
#'
#' @param objCOTAN
#' @param n
#' @param marks
#' @param expression.cl
#' @param pval.df
#' @param deltaExp
#'
#' @returns
#' @export
#'
#' @examples
FindClustersMarkers <- function(objCOTAN, n = 10, marks, expression.cl=NULL,pval.df=NULL,deltaExp=NULL,
                                method = "bonferroni"){
  marks <- unlist(markers)
  if (is.null(expression.cl)|is.null(pval.df)) {
    DEA <- DEAOnClusters(objCOTAN)
    expression.cl <- DEA$coex
    pval.df <- DEA$`p-value`
  }

  if(is.null(deltaExp)){
    deltaExp <- clustersDeltaExpression(objCOTAN)
  }

  df <- as.data.frame(matrix(data = NA,nrow = 0,ncol = 7))
  colnames(df) <- c("CL","Gene","Score","pVal","adjPVal","DEA","IsMarker")
  for (cl in unique(colnames(expression.cl))) {
    for (type in c("min","max")) {
      tmp.df <- as.data.frame(matrix(data = NA,nrow = n,ncol = 7))
      colnames(tmp.df) <- c("CL","Gene","Score","pVal","adjPVal","DEA","IsMarker")

      if (type == "max") { # Get the maximum score for each cluster
        exp <- expression.cl[order(expression.cl[,cl], decreasing = TRUE),][1:n,]
      }else if((type == "min")){ # Get the minimum score for each cluster
        exp <- expression.cl[order(expression.cl[,cl], decreasing = FALSE),][1:n,]
      }

      tmp.df$Gene <- rownames(exp)
      tmp.df$CL <- cl
      tmp.df$Score <- exp[,cl]
      tmp.df$IsMarker <- 0
      if(any(tmp.df$Gene %in% marks)){
        tmp.df[tmp.df$Gene %in% marks,]$IsMarker <- 1
      }
      tmp.df$pVal <- pval.df[tmp.df$Gene,cl]
      tmp.df$adjPVal <- p.adjust(tmp.df$pVal, n = dim(expression.cl)[1],method = method)
      tmp.df$DEA <- deltaExp[tmp.df$Gene,cl]
      df <- rbind(df, tmp.df)

    }

  }

  return(df)

}
