#' get.subset
#'
#' This function estract the row data for a subset of cluesters.
#'
#' @param object A COTAN object
#' @param cluster.names A cluster names array
#'
#' @return the subset raw count dataframe
#' @export
#'
#' @examples
#' data("ERCC.cotan")
setGeneric("get.subset", function(object, cluster.names) standardGeneric("get.subset"))
#' @rdname get.subset
setMethod("get.subset","scCOTAN",
          function(object,cluster.names) {
              meta = object@raw
              meta = meta[,colnames(meta) %in% names(object@clusters[object@clusters %in% cluster.names])]
              return(as.data.frame(meta))
          }
)
