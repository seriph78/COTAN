#' get.pval
#'
#' This function computes the p-values for genes in the COTAN object. It can be used genome-wide
#' or setting some specific genes of interest. By default it computes the p-values using the S
#' statistics (\eqn{\chi^{2}})
#' @param object a COTAN object
#' @param gene.set.col an array of genes. It will be put in columns.
#' If left empty the function will do it genome-wide.
#' @param gene.set.row an array of genes. It will be put in rows.
#' If left empty the function will do it genome-wide.
#' @param type_stat By default it computes the S (\eqn{\chi^{2}})
#'
#' @return a p-value matrix
#' @export
#' @rdname get.pval
#' @importFrom Matrix forceSymmetric
#' @examples
#'
#' data("ERCC.cotan")
#' ERCC.cotan <- get.pval(ERCC.cotan, type_stat = "S")
#'
setGeneric("get.pval", function(object, gene.set.col = c(), gene.set.row = c(), type_stat = "S") {
  standardGeneric("get.pval")
})
#' @rdname get.pval
setMethod(
  "get.pval", "scCOTAN",
  function(object, gene.set.col = c(), gene.set.row = c(), type_stat = "S") {
    print(gene.set.col)
    if (!is.null(gene.set.row)) {
      # a set for rows, not Genome Wide
      cond.row <- "on a set of genes on rows"
      stopifnot("can't have genome wide on columns and not rows! Use a
                           subset on gene.set.col, not on rows." = !is.null(gene.set.col))
      cond.col <- "on a set of genes on columns"
    } else {
      cond.row <- "genome wide on rows"
      if (is.null(gene.set.col)) {
        cond.col <- "genome wide on columns"
      } else {
        cond.col <- "on a set of genes on columns"
      }
    }
    print(paste("Get p-values", cond.col, cond.row, sep = " "))
    if (type_stat == "S") {
      print("Using function S")
      S <- get.S(object)
    } else if (type_stat == "G") {
      print("Using function G")
      S <- get.G(object)
    }


    if (cond.col == "on a set of genes on columns") {
      S <- vec2mat_rfast(S, genes = S$genes[S$genes %in% gene.set.col])
      if (cond.row == "on a set of genes on rows") {
        S <- S[rownames(S) %in% gene.set.row, ]
      }
    } else {
      S <- vec2mat_rfast(S)
    }
    p_value <- pchisq(as.matrix(S), df = 1, lower.tail = FALSE)
    return(p_value)
  }
)
