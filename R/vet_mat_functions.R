setGeneric("vec2mat_rfast", function(x) standardGeneric("vec2mat_rfast"))
setMethod("vec2mat_rfast","numeric",
vec2mat_rfast <- function(x){
    p <- sqrt(1 + 8 * length(x))/ 2 - 0.5
    m <- matrix(0, p, p)
    m <- Rfast::lower_tri.assign(m, diag = TRUE,v = x)
    m <- Rfast::upper_tri.assign(m,v = Rfast::upper_tri(Rfast::transpose(m)))
    temp.names <- stringr::str_split(names(x),pattern = "[|]",simplify = T)
    rownames(m) <- unique(temp.names[,1])
    colnames(m) <- unique(temp.names[,2])
    m
}
)

setGeneric("mat2vec_rfast", function(mat) standardGeneric("mat2vec_rfast"))
setMethod("mat2vec_rfast","numeric",
          mat2vec_rfast <- function(mat){
              v <- Rfast::lower_tri(mat,diag = T)
              names.v <- c()
              for (d in c(1:ncol(mat))) {
                  print(colnames(mat)[d])
                  temp <- paste0(rownames(mat)[d:nrow(mat)],"|",colnames(mat)[d])
                  names.v <- c(names.v,temp)
              }
              names.v
              names(v) <- names.v
              v
          }
)

