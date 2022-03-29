#' extract.coex
#'
#' This function extract a complete or a prtial coex matrix from the COTAN object.
#'
#' @param object A COTAN object
#' @param genes A vector of gene names. These will be on the columns of the final data frame. By defaults the function will use all genes.
#'
#' @return the coex values data frame
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' coex <- extract.coex(ERCC.cotan)
setGeneric("extract.coex", function(object, genes = "all") standardGeneric("extract.coex"))
#' @rdname extract.coex
setMethod("extract.coex","scCOTAN",
          function(object,genes = "all") {
              vect <- object@coex
              coex <- vec2mat_rfast(vect,genes = genes)
              return(coex)
          }
)


