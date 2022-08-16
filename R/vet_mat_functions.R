#' vec2mat_rfast
#'
#' @param x a list formed by two arrays: "genes" with the gene names and "values" with all unique values.
#' @param genes a vector with all wanted genes or "all". By default is equal to  "all" in this way it recreate the entire coex dataframe.
#'
#' @return a dataframe
#' @export
#' @importFrom Rfast lower_tri.assign
#' @importFrom Rfast upper_tri.assign
#' @importFrom Rfast transpose
#' @rdname vec2mat_rfast
#' @examples
#' v <- list("genes" = paste0("gene.", c(1:10)), "values" = c(1:55))
#' genes <- c("gene.3", "gene.4", "gene.7")
#' df <- vec2mat_rfast(v, genes)
#' df2 <- vec2mat_rfast(v)
setGeneric("vec2mat_rfast", function(x, genes = "all") standardGeneric("vec2mat_rfast"))
#' @rdname vec2mat_rfast
setMethod(
  "vec2mat_rfast", "list",
  function(x, genes = "all") {
    if (genes[1] == "all") {
      p <- sqrt(1 + 8 * length(x$values)) / 2 - 0.5
      m <- matrix(0, p, p)
      m <- Rfast::lower_tri.assign(m, diag = TRUE, v = x$values)
      m <- Rfast::upper_tri.assign(m, v = Rfast::upper_tri(Rfast::transpose(m)))
      rownames(m) <- x$genes
      colnames(m) <- x$genes
    } else {
      m <- matrix(0, nrow = length(x$genes), ncol = length(genes))
      rownames(m) <- x$genes
      colnames(m) <- genes

      for (pos.gene in match(genes, x$genes)) {
        temp.array <- x$values[pos.gene]
        p <- 1
        l <- length(x$genes)
        s <- pos.gene
        while (p < (pos.gene)) {
          l <- l - 1
          s <- s + l
          temp.array <- c(temp.array, x$values[s])
          p <- p + 1
        }


        # linear part
        start.reading.position <- 1
        i <- 1
        while (i < (pos.gene)) {
          start.reading.position <- start.reading.position + (length(x$genes) - (i - 1))
          i <- i + 1
        }
        start.reading.position <- start.reading.position + 1
        start.reading.position


        end.reading.position <- 0
        for (i in c(0:(pos.gene - 1))) {
          end.reading.position <- end.reading.position + (length(x$genes) - i)
        }
        end.reading.position

        temp.array <- c(temp.array, x$values[(start.reading.position):end.reading.position])
        m[, x$genes[pos.gene]] <- temp.array
      }
    }

    return(m)
  }
)

#' mat2vec_rfast
#'
#' @param mat a symmetric matrix with all genes as row and column names
#'
#' @return a list formed by two arrays: "genes" with the gene names and "values" with all unique values.
#' @export
#' @importFrom Rfast lower_tri.assign
#' @importFrom Rfast upper_tri.assign
#' @importFrom Rfast transpose
#' @rdname mat2vec_rfast
#' @examples
#' mat <- matrix(0, nrow = 10, ncol = 10)
#' mat <- Rfast::lower_tri.assign(mat, c(1:55), diag = TRUE)
#' mat <- Rfast::upper_tri.assign(mat, v = Rfast::upper_tri(Rfast::transpose(mat)))
#' v <- mat2vec_rfast(mat)
setGeneric("mat2vec_rfast", function(mat) standardGeneric("mat2vec_rfast"))
#' @rdname mat2vec_rfast
setMethod(
  "mat2vec_rfast", "matrix",
  function(mat) {
    if (!dim(mat)[1] == dim(mat)[2]) stop("The matrix is not simmetric!")
    v <- Rfast::lower_tri(mat, diag = TRUE)
    names.v <- rownames(mat)
    v <- list("genes" = rownames(mat), "values" = v)
    return(v)
  }
)
