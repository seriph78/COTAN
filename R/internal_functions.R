
setGeneric("get.zero_one.cells", function(object) standardGeneric("get.zero_one.cells"))
setMethod("get.zero_one.cells","scCOTAN",
          function(object) {
              cells.0.1 = as.data.frame(as.matrix(object@raw))
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells.0.1[cells.0.1 > 0] <- 1
              cells.0.1[cells.0.1 <= 0] <- 0
              # We want to discard genes having less than 3 not 0 counts over 1000 cells
              cells.0.1 <- cells.0.1[rowSums(cells.0.1) > round((length(colnames(object@raw))/1000*3),
                                                               digits = 0),]
              return(cells.0.1)
          }
)


setGeneric("fun_pzero", function(a,mu) standardGeneric("fun_pzero"))
setMethod("fun_pzero","numeric",
          function(a,mu) {
                  (a <= 0)*(exp(-(1+abs(a))*mu)) + (a > 0)*(1+abs(a)*mu)^(-1/abs(a))
              }
)

setGeneric("fun_pzero_posi", function(r,mu) standardGeneric("fun_pzero_posi"))
setMethod("fun_pzero_posi","numeric",
              function(r,mu){ (1+r*mu)^(-1/r) }
          )

setGeneric("fun_pzero_nega0", function(r,mu) standardGeneric("fun_pzero_nega0"))
setMethod("fun_pzero_nega0","numeric",
              function(r,mu){ (exp(-(1-r)*mu))}
)

setGeneric("fun_dif_mu_zeros", function(h,x,somma_zeri,mu_estimator)
    standardGeneric("fun_dif_mu_zeros"))
setMethod("fun_dif_mu_zeros","numeric",
              function(h,x,somma_zeri,mu_estimator){
                  if (h > 0) {
                      sum(fun_pzero_posi(h,mu_estimator[x,])) - somma_zeri
                  }else{
                      sum(fun_pzero_nega0(h,mu_estimator[x,])) - somma_zeri
                  }
              }
)

setGeneric("fun_my_opt", function(x,ce.mat,mu_est) standardGeneric("fun_my_opt"))
setMethod("fun_my_opt","character",
              function(x,ce.mat,mu_est){
                  somma_zeri <- sum(ce.mat[x,] == 0)

                  a1 <- 0
                  u1 <- fun_dif_mu_zeros(a1,x,somma_zeri,mu_est)
                  a2 <- a1
                  u2 <- u1
                  if (u1 > 0) {
                      a1 <- a1 - 1
                      u1 <- fun_dif_mu_zeros(a1,x,somma_zeri,mu_est)
                      while (u1 > 0) {
                          a2 <- a1
                          u2 <- u1
                          a1 <- 2 * a1
                          u1 <- fun_dif_mu_zeros(a1,x,somma_zeri,mu_est)
                      }
                  }else{
                      a2 <- 1
                      u2 <- fun_dif_mu_zeros(a2,x,somma_zeri,mu_est)
                      while (u2 < 0) {
                          a1 <- a2
                          u1 <-u2
                          a2 <- 2 * a2
                          u2 <- fun_dif_mu_zeros(a2,x,somma_zeri,mu_est)
                      }
                  }
                  a <- (a1+a2)/2
                  u <- fun_dif_mu_zeros(a,x,somma_zeri,mu_est)
                  while (abs(u)>0.001) {
                      if(u>0){
                          a2 <- a
                          u2 <- u
                      }else{
                          a1 <- a
                          u1 <- u
                      }
                      a <- (a1 + a2)/2
                      u <- fun_dif_mu_zeros(a,x,somma_zeri,mu_est)
                  }
                  r <- data.frame(a,u)
                  rownames(r) <- x
                  return(r)
              }
)


#' spMat
#'
#' Internal function to convert the matrix in a sparce trinagolar matrix
#' @param object A COTAN object
#' @importFrom Matrix forceSymmetric
#' @return the entire COTAN object
#' @noRd
#'
setGeneric("spMat", function(object) standardGeneric("spMat"))
setMethod("spMat","matrix",
          function(object){

              object[upper.tri(object,diag = FALSE)] <- 0
              object <- as(object, "triangularMatrix")
              return(object)
          }
)


#' fun_linear
#'
#' Internal function to estimate the cell efficiency
#' @param object COTAN object
#' @return a list of object (dist_cells, to_clust, pca_cells,
#' t_to_clust, mu_estimator, object)
#'
#' @importFrom Matrix rowMeans
#' @importFrom irlba prcomp_irlba
#' @importFrom Matrix colMeans
#' @importFrom stats dist
#'
setGeneric("fun_linear", function(object) standardGeneric("fun_linear"))
setMethod("fun_linear","scCOTAN",

          function(object) {
              
              print("Start estimation mu with linear method")
              print(dim(object@raw))
              genes_means <- Matrix::rowMeans(object@raw, dims = 1, na.rm = TRUE)

              max.genes <- names(sort(genes_means,decreasing = TRUE)
                        [round(length(genes_means)/2,digits = 0 ):round(length(genes_means)/4*3,
                                                                                digits = 0 )])

              cells_means <- Matrix::colMeans(object@raw, dims = 1, na.rm = TRUE)

             means <- mean(as.matrix(object@raw),na.rm = TRUE )
              mu_estimator <- (genes_means %*% t(cells_means)) /means

              rownames(mu_estimator) <- names(genes_means)

              object@nu <- colMeans(mu_estimator)/means
              object@lambda <- genes_means
              gc()
              # To insert an explorative analysis and check for strage cells (as blood)
              #and cells with a too low efficiency (nu est)


              to_clust <- Matrix::t(Matrix::t(object@raw) * (1/as.vector(object@nu)))
              
              print("Start PCA")

              p1 <- irlba::prcomp_irlba(Matrix::t(to_clust), n=5)
              
              pca_cells <- p1$x
              rm(p1) 
              gc()

             rownames(pca_cells) <- colnames(to_clust)

              ppp <- pca_cells
              ppp <- scale(ppp)
              dist_cells <- stats::dist(ppp, method = "euclidean") # mhalanobis

              rm(ppp)
              gc()
              pca_cells <- as.data.frame(pca_cells)

              output <- list("dist_cells"=dist_cells, "to_clust"=to_clust,"pca_cells"=pca_cells,
                           "mu_estimator"=mu_estimator, "object"=object)

              return(output)
          }

)


setGeneric("mu_est", function(object) standardGeneric("mu_est"))
setMethod("mu_est","scCOTAN",
          function(object) {
              print("mu estimator creation")
              mu_estimator <- object@lambda %*% t(object@nu)
              rownames(mu_estimator) <- rownames(object@raw)
              return(mu_estimator)
          }
)

setGeneric("get.S", function(object) standardGeneric("get.S"))
setMethod("get.S","scCOTAN",
          function(object) {
              print("function to generate S ")
              #S <- (object@coex)^2 * object@n_cells
              S <- (object@coex$values)^2 * object@n_cells
              S <- list("genes"=object@coex$genes,"values"=S)
              return(S)
          }
)


setGeneric("hk_genes", function(object) standardGeneric("hk_genes"))
setMethod("hk_genes","scCOTAN",
          function(object) {
              print("save effective constitutive genes")
              cells=as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              hk <- names(which(rowSums(cells) == length(colnames(cells))))
              object@hk <- hk

              return(object)
          }
)


setGeneric("obs_ct", function(object) standardGeneric("obs_ct"))
setMethod("obs_ct","scCOTAN",
          function(object) {
              cells <- as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              #cells[cells <= 0] <- 0
              cells <- as.matrix(cells)
              print("Generating contingency tables for observed data")
              somma <- rowSums(cells)
              somma <- as.matrix(somma)
              si_any <- do.call("cbind", replicate(length(rownames(somma)), somma, simplify = FALSE))

              colnames(si_any) = rownames(si_any)
              if (is.null(object@yes_yes)) {
                  object = obs_yes_yes(object)
              }

              si_si <- as.matrix(object@yes_yes)
              si_no <- si_any - si_si

              si_any <- t(si_any)
              no_si <- si_any - si_si

              no_no <- length(colnames(cells)) - (si_si + no_si + si_no)
              out <- list("no_yes"=no_si,"yes_no"=si_no,"no_no"=no_no,"object"=object)
              return(out)
          }
)

setGeneric("obs_yes_yes", function(object) standardGeneric("obs_yes_yes"))
setMethod("obs_yes_yes","scCOTAN",
          function(object) {
              print("creation of observed yes/yes contingency table")

              cells <- as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              #cells[cells <= 0] <- 0

              si_si <- as.matrix(cells) %*% t(as.matrix(cells))
              object@yes_yes <- as(as.matrix(si_si), "sparseMatrix")
              return(object)
          }
)

setGeneric("expected_ct", function(object) standardGeneric("expected_ct"))
setMethod("expected_ct","scCOTAN",
          function(object) {


              cells <- as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              mu_estimator <- mu_est(object)
              mu_estimator <- mu_estimator[!rownames(mu_estimator) %in% object@hk,]
              print("expected contingency tables creation")

              M <- fun_pzero(object@a,mu_estimator[,colnames(cells)])
              N <- 1-M
              n_zero_esti <- rowSums(M)
              # estimated number of zeros for each genes
              n_zero_obs <- rowSums(cells[!rownames(cells) %in% object@hk,] == 0)
              # observed number of zeros for each genes
              dist_zeros <- sqrt(sum((n_zero_esti - n_zero_obs)^2))

              print(paste("The distance between estimated n of zeros and observed number of zero is",
                          dist_zeros,"over", length(rownames(M)), sep = " "))

              stopifnot("Errore: some Na in matrix M " = !any(is.na(M)))

              gc()
              estimator_no_no <- M %*% t(M)
              estimator_no_si <- M %*% t(N)
              estimator_si_no <- t(estimator_no_si)
              estimator_si_si <- N %*% t(N)

              print("Done")

              out <- list("estimator_no_no"=estimator_no_no, "estimator_no_yes"=estimator_no_si,
                         "estimator_yes_no"=estimator_si_no,"estimator_yes_yes"=estimator_si_si)
              return(out)
          }
)


setGeneric("get.G", function(object) standardGeneric("get.G"))
setMethod("get.G","scCOTAN",
          function(object) {
              print("function to generate G ")
              hk <- object@hk
              ll <- obs_ct(object)

              object <- ll$object

              ll$no_yes<- ll$no_yes[!rownames(ll$no_yes) %in% hk,!colnames(ll$no_yes) %in% hk]
              ll$no_no <- ll$no_no[!rownames(ll$no_no) %in% hk,!colnames(ll$no_no) %in% hk]
              si_si <- object@yes_yes[!rownames(object@yes_yes) %in% hk,!colnames(object@yes_yes) %in% hk]
              ll$yes_no <- ll$yes_no[!rownames(ll$yes_no) %in% hk,!colnames(ll$yes_no) %in% hk]

              est <- expected_ct(object)
              for (i in est) {
                  stopifnot("Some expected values are 0!" = !any(i == 0 ))
              }

              print("G estimation")

              t1 <- as.matrix(si_si)    * log( as.matrix(si_si)         /
                                                  as.matrix(est$estimator_yes_yes))
              t1[which(as.matrix(si_si) == 0)] <- 0

              t2 <-  as.matrix(ll$no_no)  * log( as.matrix(ll$no_no) /
                                                       as.matrix(est$estimator_no_no))
              t2[which(as.matrix(ll$no_no) == 0)] <- 0

              t3 <- as.matrix(ll$yes_no) * log( as.matrix(ll$yes_no)/
                                                       as.matrix(est$estimator_yes_no))
              t3[which(as.matrix(ll$yes_no) == 0)] <- 0
              t4 <-  as.matrix(ll$no_yes) * log( as.matrix(ll$no_yes)/
                                                       as.matrix(est$estimator_no_yes))
              t4[which(as.matrix(ll$no_yes) == 0)] <- 0
              G <- 2 * (t1 + t2 + t3 + t4)
              return(G)
          }
)


