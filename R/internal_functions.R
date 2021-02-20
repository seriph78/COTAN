#' scCOTAN-class
#' Define my COTAN structure
#' @slot raw ANY. To store the raw data matrix
#' @slot raw.norm ANY. To store the raw data matrix divided for the cell efficiency estimated (nu)
#' @slot coex ANY. The coex matrix (sparce)
#' @slot nu vector.
#' @slot lambda vector.
#' @slot a vector.
#' @slot hk vector.
#' @slot n_cells numeric.
#' @slot meta data.frame.
#' @slot yes_yes ANY.
#' @slot clusters vector.
#' @slot cluster_data data.frame.
#'
#' @return
#' @export
#'
#' @examples obj = scCOTAN(raw = raw) # with raw the raw data matrix with cells on colums and genes on rows
scCOTAN <- setClass("scCOTAN", slots =
                        c(raw="ANY",raw.norm="ANY", coex="ANY",
                          nu="vector",lambda="vector",a="vector", hk = "vector", n_cells = "numeric",
                          meta="data.frame",yes_yes="ANY", clusters="vector",cluster_data="data.frame"))



# internal functions

mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                 axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                 axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                 axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))




#' spMat
#' Internal function to convert the matrix in a sparce trinagolar matrix
#' @param object
#'
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
#' Matrix
#' basilisk
#' @importFrom parallel mclapply
#'
#'
setGeneric("fun_linear", function(object, mean_type) standardGeneric("fun_linear"))
setMethod("fun_linear","scCOTAN",
          #fun_linear =
          function(object, mean_type) {


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
              file.py <- system.file("inst","python", "python_PCA.py", package="COTAN")
              pca_cells <- basiliskRun(proc, function(arg1) {
                  reticulate::source_python(file.py)
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


#' fun_sqrt
#' Internal function to estimeate the cell efficiency with the square root
#'
#' @param object It is COTAN object
#'
#' @return
#' @import
#' Matrix
#' reticulate
#'
setGeneric("fun_sqrt", function(object) standardGeneric("fun_sqrt"))
setMethod("fun_sqrt","scCOTAN",
          function(object){
              library(propagate)
              c1 = 6.5
              c2 = 60
              ## psi functions
              psi <- function(x){
                  val_psi <- x**2 + psi_1(x) * psi_2(x) + psi_3(x)
                  return(val_psi)
              }

              psi_1 <- function(x){
                  val_psi_1 = 1/4*x^7 + 3/32*x^5 + x - x^2/sqrt(2)
                  return(val_psi_1)
              }

              psi_2 <- function(x){
                  val_psi_2 = 1 / (x^7 + 1)
                  return(val_psi_2)
              }

              psi_3 <- function(x){
                  val_psi_3 = c2 * x^8 * exp(-c1*x)
                  return(val_psi_3)
              }

              ## derivates of psi
              first_der_psi <- function(x){
                  f_psi = 2 * x + first_der_psi_1(x) * psi_2(x) + psi_1(x) * first_der_psi_2(x) + first_der_psi_3(x)
                  return(f_psi)
              }

              sec_der_psi <- function(x){
                  s_psi = 2 + sec_der_psi_1(x) * psi_2(x) + 2 * first_der_psi_1(x) * first_der_psi_2(x) + psi_1(x) * sec_der_psi_2(x) + sec_der_psi_3(x)
                  return(s_psi)
              }

              first_der_psi_1 <- function(x){
                  f_psi_1 = 7/4 * x^6 + 15/32 * x^4 + 1 - sqrt(2) * x
                  return(f_psi_1)
              }

              sec_der_psi_1 <- function(x){
                  s_psi_1 = 21/2 * x^5 + 15/8 * x^3 - sqrt(2)
                  return(s_psi_1)
              }

              first_der_psi_2 <- function(x){
                  f_psi_2 = -(7 * x^6)/(x^7 + 1)^2
                  return(f_psi_2)
              }

              sec_der_psi_2 <- function(x){
                  s_psi_2 = -(14 * x^5 * (3 - 4*x^7))/(x^7 + 1)^3
                  return(s_psi_2)
              }

              first_der_psi_3 <- function(x){
                  f_psi_3 = c2*(8*x^7 - c1*x^8)*exp(-c1*x)
                  return(f_psi_3)
              }

              sec_der_psi_3 <- function(x){
                  s_psi_3 = c2 * (56 * x^6 - 16 * c1 * x^7 + c1^2 * x^8) * exp(-c1*x)
                  return(s_psi_3)
              }

              # tau
              tau <- function(x){
                  xtau = psi(x) - x^2
                  return(xtau)
              }

              # lambda
              lambda <-function(xi_star, mat_X){
                  v = rowVarsC(mat_X)
                  l_i = psi(xi_star) + 1/2 * sec_der_psi(xi_star)*(v - tau(xi_star))
                  return(l_i)
              }

              # nu
              R_star_j <- function(x_star_j, mat_X){
                  v = colVarsC(mat_X)
                  Rsj = psi(x_star_j) + 1/2 * sec_der_psi(x_star_j)*(v - tau(x_star_j))
                  return(Rsj)
              }

              nu <- function(x_star_j,mat_X){
                  m = length(colnames(mat_X))
                  Rsj = R_star_j(x_star_j, mat_X)
                  nu_j = Rsj/(sum(Rsj)/m)
                  return(nu_j)
              }


              #library(reticulate)
              #setwd("../seriph/COTAN-tool/")
              #source_python(paste(surcedir,"python_PCA.py",sep = ""))

              #print("Start estimation mu")
              # Estimators computation
              raw = object@raw

              new_raw = sqrt(raw)

              print("lambda i")
              lambda_i = lambda(rowMeans(new_raw), new_raw)
              lamba.data = data.frame("means"= rowMeans(new_raw), "var"=rowVarsC(new_raw))
              print("nu est")
              nu_est = nu(colMeans(new_raw),new_raw)
              nu.data = data.frame("means"= colMeans(new_raw), "var"=colVarsC(new_raw))
              print("mu estimator")
              mu_estimator = outer(lambda_i,nu_est, "*")

              object@nu = nu_est
              object@lambda = lambda_i


              gc()
              # To insert an explorative analysis and check for strage cells (as blood) and cells with a too low efficiency (nu est)
              to_clust <- t(t(as.matrix(raw)) * (1/as.vector(nu_est))) #raw counts divided for cell efficiency
              colnames(to_clust) = colnames(raw)
              t_to_clust = t(to_clust)
              #start_time <- Sys.time()
              t_to_clust = as.matrix(t_to_clust)
              proc <- basiliskStart(my_env_cotan)
              on.exit(basiliskStop(proc))
              t_to_clust = as.matrix(t_to_clust)
              print(getwd())
              pca_cells <- basiliskRun(proc, function(arg1) {
                  reticulate::source_python("../../inst/python/python_PCA.py")
                  output <- python_PCA(arg1)

                  # The return value MUST be a pure R object, i.e., no reticulate
                  # Python objects, no pointers to shared memory.
                  output
              }, arg1=t_to_clust)
              rownames(pca_cells)=rownames(t_to_clust)
              #end_time <- Sys.time()
              #print(paste("pca; time",end_time - start_time, sep = " " ))

              #---- Mhalanobis distance
              ppp = pca_cells
              ppp = scale(ppp)
              dist_cells = dist(ppp, method = "euclidean") # mhalanobis
              colnames(pca_cells) = paste("PC",c(1:ncol(pca_cells)), sep = "")
              pca_cells = as.data.frame(pca_cells)
              output = list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells, "t_to_clust"=t_to_clust,"mu_estimator"=mu_estimator, "object"=object,"nu.data"=nu.data,"lamba.data"= lamba.data)

              #detach("package:reticulate", unload = TRUE)
              detach("package:propagate", unload = TRUE)
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


