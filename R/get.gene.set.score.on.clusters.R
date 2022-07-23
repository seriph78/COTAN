#' get.gene.set.score.on.clusters
#'
#'This function compute a wheigted parcentage of the expressed genes for each field 
#'in the gene list (\eqn{\frac{1}{n}\sum_i(1-e^{-\theta X_i})} where the \eqn{X_i} 
#'are the values from the function get.clusters.delta.expression and \eqn{\theta = -1/0.1 \cdot ln(0.25)})
#'
#' @param obj a cotan object
#' @param genes.list a list with genes 
#' @param out_dir where to save the output
#' @param cond prefix for the saved csv
#'
#' @return a dataframe
#' @export
#'
#' @examples
setGeneric("get.gene.set.score.on.clusters", function(obj,genes.list, out_dir,cond) 
  standardGeneric("get.gene.set.score.on.clusters"))
#' @rdname get.gene.set.score.on.clusters
setMethod("get.gene.set.score.on.clusters","scCOTAN",
          function(obj,genes.list, out_dir,cond) {
            # cluster assignment with expression increment
            expression.cl <- get.clusters.delta.expression(obj)
            df <- as.data.frame(matrix(nrow = length(names(markers)),ncol = ncol(obj@cluster_data)))
            rownames(df) <- names(markers)
            colnames(df) <- colnames(obj@cluster_data)
            teta <- -1/0.1 * (log(0.25))
            #not_ass_clusters <- NA
            for (ro in colnames(df)) {
              for (m in names(markers)) {
                #pv <- p_value[unlist(markers[[m]]),ro]
                #co <- obj@cluster_data[unlist(markers[[m]]),ro]
                ex <- expression.cl[unlist(markers[[m]]),ro]
                ex[ex < 0 & !is.na(ex)] <- 0
                
                ex <- 1-exp(- teta * ex)
                n.markers <- sum(unlist(markers[[m]]) %in% rownames(obj@raw))
                df[m,ro] <- sum(ex,na.rm = T)/n.markers
              }
            }
            write.csv(df,file = paste(out_dir,cond,"_gene_set_scores.csv", sep = ""))
            
            return(df)
                     }
)