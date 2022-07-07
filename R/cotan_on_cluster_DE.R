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
setGeneric("DEA_on_clusters", function(obj, cells_list) 
  standardGeneric("DEA_on_clusters"))
#' @rdname DEA_on_clusters
setMethod("DEA_on_clusters","scCOTAN",
          function(obj, cells_list) {

              #---------------------------------------------
              fun_pzero <- function(a,mu){
                  (a <= 0)*(exp(-(1+abs(a))*mu)) + (a > 0)*(1+abs(a)*mu)^(-1/abs(a))
              }

              fun_pzero_posi <- function(r,mu){ (1+r*mu)^(-1/r) }

              fun_pzero_nega0 <- function(r,mu){ (exp(-(1-r)*mu))}

              fun_dif_mu_zeros <- function(h,x,somma_zeri){
                  if (h > 0) {
                      sum(fun_pzero_posi(h,mu_estimator[x,])) - somma_zeri#/somma_zeri
                  }else{
                      sum(fun_pzero_nega0(h,mu_estimator[x,])) - somma_zeri#/somma_zeri
                  }
              }

              fun_my_opt <- function(x){
                  somma_zeri = sum(cells[x,] == 0)
                  a1 = 0
                  u1 = fun_dif_mu_zeros(a1,x,somma_zeri)
                  a2 = a1
                  u2 = u1
                  if (u1 > 0) {
                      a1 = a1 - 1
                      u1 = fun_dif_mu_zeros(a1,x,somma_zeri)
                      while (u1 > 0) {
                          a2 = a1
                          u2 = u1
                          a1 = 2 * a1
                          u1 = fun_dif_mu_zeros(a1,x,somma_zeri)
                      }
                  }else{
                      a2 = 1
                      u2 = fun_dif_mu_zeros(a2,x,somma_zeri)
                      while (u2 < 0) {
                          a1 = a2
                          u1 =u2
                          a2 = 2 * a2
                          u2 = fun_dif_mu_zeros(a2,x,somma_zeri)
                      }
                  }
                  a = (a1+a2)/2
                  u = fun_dif_mu_zeros(a,x,somma_zeri)
                  while (abs(u)>0.001) {
                      if(u>0){
                          a2 = a
                          u2 = u
                      }else{
                          a1 = a
                          u1 = u
                      }
                      a = (a1 + a2)/2
                      u = fun_dif_mu_zeros(a,x,somma_zeri)
                  }
                  r = data.frame(a,u)
                  rownames(r) = x
                  return(r)
              }


              cluster_data = data.frame()
              cluster_pval = data.frame()
              # cells_set  > set of cell code corresponding to cluster
              for (condition in names(cells_list)) {
                gc()
                  print(paste("cluster",condition, sep=" "))
                  cells=as.matrix(obj@raw)
                  cells_set = unlist(cells_list[condition])

                  stopifnot("ERROR. Some cells are not present!"=  all(cells_set %in% colnames(cells)) )

                  # Cells matrix : formed by row data matrix changed to 0-1 matrix
                  cells[cells > 0] <- 1

                  mu_estimator = mu_est(obj)

                  hk = obj@hk

                  mu_estimator = mu_estimator[!rownames(mu_estimator) %in% hk,]
                  cells = cells[!rownames(cells) %in% hk, ]

                  yes_in = rowSums(cells[,colnames(cells) %in% cells_set])
                  yes_out = rowSums(cells[,!colnames(cells) %in% cells_set])
                  no_in = rowSums(1 - cells[,colnames(cells) %in% cells_set])
                  no_out = rowSums(1 - cells[,!colnames(cells) %in% cells_set])



                  M = fun_pzero(obj@a,mu_estimator[,colnames(cells)]) # matrix of 0 probabilities
                  estimator_no_in = rowSums(M[,colnames(cells) %in% cells_set])
                  estimator_no_out = rowSums(M[,!colnames(cells) %in% cells_set])

                  N = 1-M
                  estimator_yes_in = rowSums(N[,colnames(cells) %in% cells_set])
                  estimator_yes_out = rowSums(N[,!colnames(cells) %in% cells_set])

                  if( (sum(((yes_in+no_in)-(estimator_no_in+estimator_yes_in))**2)+
                       sum(((yes_out+no_out)-(estimator_yes_out+estimator_no_out))**2)) > 0.0001){
                      print("Problems with estimators!")
                  }


                  exp_yes = rowSums(cells)
                  exp_no =obj@n_cells - exp_yes
                  if( sum(yes_in+yes_out-exp_yes) != 0 |  sum(no_in+no_out-exp_no) != 0 ){
                      print("Problems with observed counts!")
                  }


                  new_estimator_no_in = estimator_no_in
                  new_estimator_no_in[new_estimator_no_in < 1] <- 1
                  new_estimator_no_out = estimator_no_out
                  new_estimator_no_out[new_estimator_no_out < 1] <- 1

                  new_estimator_yes_in = estimator_yes_in
                  new_estimator_yes_in[new_estimator_yes_in < 1] <- 1
                  new_estimator_yes_out = estimator_yes_out
                  new_estimator_yes_out[new_estimator_yes_out < 1] <- 1

                  dif_no_in = (as.matrix(no_in) - estimator_no_in)**2/new_estimator_no_in
                  dif_no_out = (as.matrix(no_out) - estimator_no_out)**2/new_estimator_no_out
                  dif_yes_in = (as.matrix(yes_in) - estimator_yes_in)**2/new_estimator_yes_in
                  dif_yes_out = (as.matrix(yes_out) - estimator_yes_out)**2/new_estimator_yes_out

                  S = dif_no_in + dif_no_out + dif_yes_in + dif_yes_out
                  
                  rm(dif_yes_out,dif_no_in,dif_no_out,dif_yes_in)
                  gc()

                  if(any(is.na(S))){
                      print(paste("Errore: some Na in matrix S", which(is.na(S),arr.ind = T),sep = " "))
                      break()
                  }


                  print("Calculating p values")
                  p_value = as.data.frame(pchisq(as.matrix(S), df=1, lower.tail=F))

                  gc()

                  coex = ((as.matrix(yes_in) - as.matrix(estimator_yes_in))/as.matrix(new_estimator_yes_in)) +
                      ((as.matrix(no_out) - as.matrix(estimator_no_out))/as.matrix(new_estimator_no_out)) -
                      ((as.matrix(yes_out) - as.matrix(estimator_yes_out))/as.matrix(new_estimator_yes_out)) -
                      ((as.matrix(no_in) - as.matrix(estimator_no_in))/as.matrix(new_estimator_no_in))
                  
                  rm(yes_in, yes_out,no_out,no_in)
                  rm(estimator_yes_out,estimator_yes_in, estimator_no_out, estimator_no_in)
                  gc()

                  coex = coex / sqrt(1/new_estimator_yes_in + 1/new_estimator_no_in + 1/new_estimator_yes_out + 1/new_estimator_no_out)

                  coex = coex/sqrt(obj@n_cells)
                  coex = as.data.frame(coex)

                  colnames(coex)=paste("cl.",condition,sep = "")
                  colnames(p_value)=paste("cl.",condition,sep = "")
                  #to insert gene names when the dataframe is empty

                  if (dim(cluster_data)[1] == 0) {
                      cluster_data <- as.data.frame(matrix(nrow = length(obj@coex$genes)))
                      rownames(cluster_data) <- obj@coex$genes
                      cluster_pval <- as.data.frame(matrix(nrow = length(obj@coex$genes)))
                      rownames(cluster_pval) <- obj@coex$genes
                  }

                  cluster_data = cbind(cluster_data,coex)
                  cluster_pval = cbind(cluster_pval,p_value)

              }
              obj@cluster_data = cluster_data[,2:ncol(cluster_data)]

              return(list(obj,cluster_pval))
          }
)
