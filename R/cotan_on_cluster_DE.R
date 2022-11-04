#' DEA_on_clusters
#' Next method is used to do the differential expression analysis using COTAN contingency tables
#' on the cluster checked with the previous method
#'
#' @param obj A COTAN object
#' @param cells_list a list with all the cells divided for each group. Name of each element in the list will
#' be the condition o cell population or cluster number, the arraz will contain all the cell codes (as in the
#' raw matrix column names).
#'
#' @return a list with two objects: the first is a scCOTAN with the new 
#' object with also the correlation matrix for the genes in each cluster (obj@cluster_data), 
#' the second is the p-values matrix. 
#' @export
#'
#' @examples
setGeneric("DEA_on_clusters", function(obj, cells_list = NULL) 
  standardGeneric("DEA_on_clusters"))
#' @rdname DEA_on_clusters
setMethod("DEA_on_clusters","scCOTAN",
          function(obj, cells_list = NULL) {
            
            if (is.null(cells_list)) {
              clusters.names = unique(obj@clusters)[!is.na(unique(obj@clusters))]
              cells_list = list(names(obj@clusters[obj@clusters %in% clusters.names[1]]))
              names(cells_list)=clusters.names[1]
              for (c in c(2:length(clusters.names))) {
                tmp = list(names(obj@clusters[obj@clusters %in% clusters.names[c]]))
                names(tmp)= clusters.names[c]
                cells_list = c(cells_list,tmp)
              }  
            }
            
              cells<-obj@raw
              #cells[cells > 0] <- 1
              
              mu_estimator <- estimateMu(obj)
              
              hk <- obj@hk
              
              mu_estimator <- mu_estimator[!rownames(mu_estimator) %in% hk,]
              cells <- cells[!rownames(cells) %in% hk, ]
              M <- funProbZero(obj@a,mu_estimator[,colnames(cells)]) # matrix of 0 probabilities
              N <- 1-M
              
              cluster_data <- data.frame()
              cluster_pval <- data.frame()
              # cells_set  > set of cell code corresponding to cluster
              for (condition in names(cells_list)) {
                gc()
                  print(paste("cluster",condition, sep<-" "))
                  
                  cells_set <- unlist(cells_list[condition])

                  stopifnot("ERROR. Some cells are not present!"<-  all(cells_set %in% colnames(cells)) )

                  # Cells matrix : formed by row data matrix changed to 0-1 matrix
                  

                  yes_in <- rowSums(cells[,colnames(cells) %in% cells_set]>0)
                  yes_out <- rowSums(cells[,!colnames(cells) %in% cells_set]>0)
                  #no_in <- rowSums(1 - cells[,colnames(cells) %in% cells_set] )
                  #no_out <- rowSums(1 - cells[,!colnames(cells) %in% cells_set])
                  no_in <- rowSums(cells[,colnames(cells) %in% cells_set] == 0)
                  no_out <- rowSums(cells[,!colnames(cells) %in% cells_set] == 0)



                  estimator_no_in <- rowSums(M[,colnames(cells) %in% cells_set])
                  estimator_no_out <- rowSums(M[,!colnames(cells) %in% cells_set])

                  
                  estimator_yes_in <- rowSums(N[,colnames(cells) %in% cells_set])
                  estimator_yes_out <- rowSums(N[,!colnames(cells) %in% cells_set])

                  if( (sum(((yes_in+no_in)-(estimator_no_in+estimator_yes_in))**2)+
                       sum(((yes_out+no_out)-(estimator_yes_out+estimator_no_out))**2)) > 0.0001){
                      print("Problems with estimators!")
                  }


                  exp_yes <- rowSums(cells>0)
                  exp_no <-obj@n_cells - exp_yes
                  if( sum(yes_in+yes_out-exp_yes) != 0 |  sum(no_in+no_out-exp_no) != 0 ){
                      print("Problems with observed counts!")
                  }


                  new_estimator_no_in <- estimator_no_in
                  new_estimator_no_in[new_estimator_no_in < 1] <- 1
                  new_estimator_no_out <- estimator_no_out
                  new_estimator_no_out[new_estimator_no_out < 1] <- 1

                  new_estimator_yes_in <- estimator_yes_in
                  new_estimator_yes_in[new_estimator_yes_in < 1] <- 1
                  new_estimator_yes_out <- estimator_yes_out
                  new_estimator_yes_out[new_estimator_yes_out < 1] <- 1

                  dif_no_in <- (as.matrix(no_in) - estimator_no_in)**2/new_estimator_no_in
                  dif_no_out <- (as.matrix(no_out) - estimator_no_out)**2/new_estimator_no_out
                  dif_yes_in <- (as.matrix(yes_in) - estimator_yes_in)**2/new_estimator_yes_in
                  dif_yes_out <- (as.matrix(yes_out) - estimator_yes_out)**2/new_estimator_yes_out

                  S <- dif_no_in + dif_no_out + dif_yes_in + dif_yes_out
                  
                  rm(dif_yes_out,dif_no_in,dif_no_out,dif_yes_in)
                  gc()

                  if(any(is.na(S))){
                      print(paste("Errore: some Na in matrix S", which(is.na(S),arr.ind = T),sep = " "))
                      break()
                  }


                  print("Calculating p values")
                  p_value <- as.data.frame(pchisq(as.matrix(S), df<-1, lower.tail<-F))

                  gc()

                  coex <- ((as.matrix(yes_in) - as.matrix(estimator_yes_in))/as.matrix(new_estimator_yes_in)) +
                      ((as.matrix(no_out) - as.matrix(estimator_no_out))/as.matrix(new_estimator_no_out)) -
                      ((as.matrix(yes_out) - as.matrix(estimator_yes_out))/as.matrix(new_estimator_yes_out)) -
                      ((as.matrix(no_in) - as.matrix(estimator_no_in))/as.matrix(new_estimator_no_in))
                  
                  rm(yes_in, yes_out,no_out,no_in)
                  rm(estimator_yes_out,estimator_yes_in, estimator_no_out, estimator_no_in)
                  gc()

                  coex <- coex / sqrt(1/new_estimator_yes_in + 1/new_estimator_no_in + 1/new_estimator_yes_out + 1/new_estimator_no_out)

                  coex <- coex/sqrt(obj@n_cells)
                  coex <- as.data.frame(coex)

                  colnames(coex) <- paste("cl.",condition,sep = "")
                  colnames(p_value) <- paste("cl.",condition,sep = "")
                  #to insert gene names when the dataframe is empty

                  if (dim(cluster_data)[1] == 0) {
                      cluster_data <- as.data.frame(matrix(nrow = length(obj@coex$genes)))
                      rownames(cluster_data) <- obj@coex$genes
                      cluster_pval <- as.data.frame(matrix(nrow = length(obj@coex$genes)))
                      rownames(cluster_pval) <- obj@coex$genes
                  }

                  cluster_data <- cbind(cluster_data,coex)
                  cluster_pval <- cbind(cluster_pval,p_value)

              }
              obj@cluster_data <- cluster_data[,2:ncol(cluster_data)]

              return(list(obj,cluster_pval))
          }
)
