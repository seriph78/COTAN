#' cosine.dissimilarity
#'
#' @param mat a matrix
#'
#' @return The dissimilarity matrix between column data
#' @export
#' @rdname cosine.dissimilarity
#' @examples
#' mat <- matrix(c(1:25),nrow = 5, ncol = 5)
#' colnames(mat) <- paste0("col.",c(1:5))
#' rownames(mat) <- paste0("row.",c(1:5))
#' dis <-  cosine.dissimilarity(mat)
setGeneric("cosine.dissimilarity", function(mat)
    standardGeneric("cosine.dissimilarity"))
#' @rdname cosine.dissimilarity
setMethod("cosine.dissimilarity","ANY",

          function(mat) {
              Matrix <- as.matrix(t(mat))
              sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
              sim <- sim %*% t(sim)
              D_sim <- as.dist(1 - sim)
              return(D_sim)
          }
)

