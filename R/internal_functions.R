getPyPath <- function() system.file("inst","python", "python_PCA.py", package="COTAN",mustWork = T)


setGeneric("fun_pzero", function(a,mu) standardGeneric("fun_pzero"))
setMethod("fun_pzero","numeric",
          function(a,mu) {
              #---------------------------------------------
              #fun_pzero <- function(a,mu){
                  #a= as.numeric(a[,2])
                  #print(a)
                  (a <= 0)*(exp(-(1+abs(a))*mu)) + (a > 0)*(1+abs(a)*mu)^(-1/abs(a))
              }
)

setGeneric("fun_pzero_posi", function(r,mu) standardGeneric("fun_pzero_posi"))
setMethod("fun_pzero_posi","numeric",
              #fun_pzero_posi <-
              function(r,mu){ (1+r*mu)^(-1/r) }
          )

setGeneric("fun_pzero_nega0", function(r,mu) standardGeneric("fun_pzero_nega0"))
setMethod("fun_pzero_nega0","numeric",
              #fun_pzero_nega0 <-
              function(r,mu){ (exp(-(1-r)*mu))}
)

setGeneric("fun_dif_mu_zeros", function(h,x,somma_zeri) standardGeneric("fun_dif_mu_zeros"))
setMethod("fun_dif_mu_zeros","numeric",
              #fun_dif_mu_zeros <-
              function(h,x,somma_zeri){
                  if (h > 0) {
                      sum(fun_pzero_posi(h,mu_estimator[x,])) - somma_zeri#/somma_zeri
                  }else{
                      sum(fun_pzero_nega0(h,mu_estimator[x,])) - somma_zeri#/somma_zeri
                  }
              }
)

setGeneric("fun_my_opt", function(x) standardGeneric("fun_my_opt"))
setMethod("fun_my_opt","character",
              #fun_my_opt <-
              function(x){
                  somma_zeri = sum(cells[x,] == 0)
                  #somma_zeri = rowSums(cells[x,] == 0)
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
)






# internal functions

mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                 axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                 axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                 axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))




#' spMat
#' Internal function to convert the matrix in a sparce trinagolar matrix
#' @param object A COTAN object
#' @importFrom Matrix forceSymmetric
#' @return the entire COTAN object
#' @noRd
#'
setGeneric("spMat", function(object) standardGeneric("spMat"))
setMethod("spMat","matrix",
          function(object){
              #object =as.matrix(object)
              object[upper.tri(object,diag = F)] = 0
              object = as(object, "triangularMatrix")
              #---------------------------

              # To symmetrize
              object =forceSymmetric(object,uplo="L")
              return(object)
          }
)


#' fun_linear
#' Internal function to estimeate the cell efficiency
#' @param object COTAN object
#' @param mean_type string: "restricted" or "normal". Default is normal. The other was wsed only for tests.
#'
#' @return
#' @import
#' basilisk
#' reticulate
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colMeans
#' @importFrom parallel mclapply
#'
#'
setGeneric("fun_linear", function(object, mean_type) standardGeneric("fun_linear"))
setMethod("fun_linear","scCOTAN",
          #fun_linear =
          function(object, mean_type) {
              #file.py <- system.file("inst","python", "python_PCA.py", package="COTAN",mustWork = T)

              print("Start estimation mu with linear method")
              print(dim(object@raw))
              genes_means = Matrix::rowMeans(object@raw, dims = 1, na.rm = T)

              max.genes = names(sort(genes_means,decreasing = T)
                                [round(length(genes_means)/2,digits = 0 ):round(length(genes_means)/4*3,digits = 0 )])

              if (mean_type == "restricted") {
                  print("cells mean type: restricted")
                  cells_means = Matrix::colMeans(object@raw[rownames(object@raw) %in% max.genes,], dims = 1, na.rm = T)
              }else{
                  print("cells mean type: normal")
                  cells_means = Matrix::colMeans(object@raw, dims = 1, na.rm = T)
              }


              means = mean(as.matrix(object@raw),na.rm = T )
              print("pippo")
              mu_estimator = (genes_means %*% t(cells_means)) /means

              rownames(mu_estimator) = names(genes_means)

              #nu_est = colMeans(mu_estimator)/means
              object@nu = colMeans(mu_estimator)/means
              object@lambda = genes_means

              #print(paste("End estimation; time",end_time - start_time, sep = " " ))
              gc()
              # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
              #start_time <- Sys.time()

              to_clust <- t(t(as.matrix(object@raw)) * (1/as.vector(object@nu)))
              #to_clust <- as(as.matrix(to_clust), "sparseMatrix")


              #end_time <- Sys.time()
              #print(paste("to clust; time",end_time - start_time, sep = " " ))
              #object@raw.norm = to_clust
              t_to_clust = t(to_clust)

              #start_time <- Sys.time()
              #to import using Basilisk
              proc <- basiliskStart(my_env_cotan)
              on.exit(basiliskStop(proc))

              t_to_clust = as.matrix(t_to_clust)

              pca_cells <- basiliskRun(proc, function(arg1) {
                  reticulate::source_python(getPyPath())
                  output <- python_PCA(arg1)

                  # The return value MUST be a pure R object, i.e., no reticulate
                  # Python objects, no pointers to shared memory.
                  output
              }, arg1=t_to_clust)

              #pca_cells = python_PCA(as.matrix(t_to_clust))
              rownames(pca_cells)=rownames(t_to_clust)
              #end_time <- Sys.time()
              #print(paste("pca; time",end_time - start_time, sep = " " ))
              #dist_cells = dist(t_to_clust, method = "euclidean")
              #---- Mhalanobis distance
              ppp = pca_cells
              ppp = scale(ppp)
              dist_cells = dist(ppp, method = "euclidean") # mhalanobis
              colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "")
              pca_cells = as.data.frame(pca_cells)

              output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "object"=object)#, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )
              #detach("package:reticulate", unload = TRUE)
              return(output)


          }

)



setGeneric("fun_linear_iter", function(object,mean_t) standardGeneric("fun_linear_iter"))
setMethod("fun_linear_iter","scCOTAN",
          function(object, mean_t) {
              # <- function(cells, raw){
              #library(reticulate)

              #source_python(paste(surcedir,"python_PCA.py",sep = ""))
              cells=as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              print("Start estimation mu with average method with iter")

              object@raw = as.matrix(object@raw[rownames(cells),colnames(cells)])
              new_object = object@raw

              genes_on_off = (matrix(ncol = length(colnames(new_object)), nrow = length(rownames(new_object)),data = TRUE))
              rownames(genes_on_off) = rownames(new_object)
              colnames(genes_on_off) = colnames(new_object)
              lambda_i = rowMeans(new_object, dims = 1, na.rm = T)
              free_deg = 0
              old_free_deg = 1
              start_time <- Sys.time()
              it =1
              while (free_deg != old_free_deg) {
                  start_time_it <- Sys.time()
                  old_free_deg = free_deg

                  cells_means = colMeans(new_object, dims = 1, na.rm = T)
                  genes_means = rowMeans(new_object, dims = 1, na.rm = T)
                  means = mean(new_object,na.rm = T )

                  end_time <- Sys.time()
                  print(paste("after means evaluation",end_time - start_time, sep = " "))
                  mu_estimator = (genes_means %*% t(cells_means)) /means
                  end_time <- Sys.time()
                  print(paste("after putting together",end_time - start_time, sep = " "))
                  rownames(mu_estimator) = names(genes_means)
                  colnames(mu_estimator) = names(cells_means)
                  #print(mu_estimator[1:3,1:3])
                  mu_estimator = as.data.frame(mu_estimator)
                  #print(mu_estimator[1:3,1:3])
                  genes_on_off = genes_on_off & !(exp(-mu_estimator) < 10**-3 & object@raw[rownames(cells),colnames(cells)] == 0 ) #& object@raw[rownames(cells),] == 0 )

                  free_deg = sum(genes_on_off == FALSE, na.rm = T)
                  print(paste("free deg",free_deg, sep = " "))

                  new_object[!genes_on_off] = mu_estimator[!genes_on_off]
                  nu_est = colMeans(mu_estimator)/means

                  end_time <- Sys.time()
                  print(paste("after iteration",it,"time",end_time - start_time_it, sep = " "))
                  it = it + 1
                  gc()
              }
              end_time <- Sys.time()
              print(paste("End estimation; time",end_time - start_time, sep = " " ))
              gc()
              # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
              start_time <- Sys.time()

              to_clust <- t(t(as.matrix(object@raw)) * (1/as.vector(nu_est)))

              end_time <- Sys.time()
              print(paste("to clust; time",end_time - start_time, sep = " " ))

              t_to_clust = t(to_clust)



              t_to_clust = as.matrix(t_to_clust)
              print(getwd())

              pca_cells <- basiliskRun(proc, function(arg1) {
                  reticulate::source_python("inst/python/python_PCA.py")
                  output <- python_PCA(arg1)

                  # The return value MUST be a pure R object, i.e., no reticulate
                  # Python objects, no pointers to shared memory.
                  output
              }, arg1=t_to_clust)
              rownames(pca_cells)=rownames(t_to_clust)

              rownames(pca_cells)=rownames(t_to_clust)
              end_time <- Sys.time()
              print(paste("pca; time",end_time - start_time, sep = " " ))
              #---- Mhalanobis distance
              ppp = pca_cells
              ppp = scale(ppp)
              dist_cells = dist(ppp, method = "euclidean") # mhalanobis
              colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "")
              pca_cells = as.data.frame(pca_cells)
              #output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "nu_est"=nu_est, "lambda_i"=lambda_i)#, "si_si"=si_si, "si_no"=si_no,"no_no"=no_no,"no_si"=no_si,"genes_on_off"=genes_on_off )

              object@nu = nu_est
              object@lambda = lambda_i

              output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "object"=object,"nu.data"=nu_est,"lamba.data"= lambda_i)

              #detach("package:reticulate", unload = TRUE)
              #detach("package:propagate", unload = TRUE)
              return(output)
          }
)

setGeneric("mu_est", function(object) standardGeneric("mu_est"))
setMethod("mu_est","scCOTAN",
          function(object) {
              print("mu estimator creation")

              #cells_means = colMeans(as.matrix(object@raw), dims = 1, na.rm = T)
              #genes_means = rowMeans(as.matrix(object@raw), dims = 1, na.rm = T)
              #means = mean(as.matrix(object@raw),na.rm = T )

              #mu_estimator = (genes_means %*% t(cells_means)) /means
              mu_estimator = object@lambda %*% t(object@nu)

              rownames(mu_estimator) = rownames(object@raw)
              #rownames(mu_estimator) = rownames(object@raw.norm)
              #colnames(mu_estimator) = colnames(object@raw)

              return(mu_estimator)
          }
)

setGeneric("get.S", function(object) standardGeneric("get.S"))
setMethod("get.S","scCOTAN",
          function(object) {
              print("function to generate S ")

              S = (object@coex)^2 * object@n_cells

              return(S)
          }
)
