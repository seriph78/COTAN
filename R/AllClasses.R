#' Definition of COTAN class
#' @slot raw raw UMI count matrix 𝑛×𝑚 (gene number × cell number)
#' @slot rawNorm raw UMI count matrix divided by UDE, 𝑛×𝑚
#' @slot coex correlation of COTAN between genes, 𝑛×𝑛
#' @slot nu vector that stores the estimated UDE, size 𝑚
#' @slot lambda vector to store the average for the gene expression, size 𝑛
#' @slot dispetrion vector to store all
#' the negative binomial dispersion factors, size 𝑛.
#' @slot hKGenes house-keeping genes. It is a vector to store the name 
#' of the genes with positive UMI count in every single cell of the sample
#' @slot metaDataset data.frame
#' @slot metaCells data.frame
#' @slot clustersCoex coex
setClass("COTAN",
         slots = c(
           raw          = "dgCMatrix",
           rawNorm      = "dgCMatrix",
           coex         = "ANY",
           cellsCoex    = "ANY",
           nu           = "vector",
           lambda       = "vector",
           dispersion   = "vector",
           hkGenes      = "vector",
           metaDataset  = "data.frame",
           metaCells    = "data.frame",
           clustersCoex = "list"
         ),
         prototype=list(
           raw          = as(matrix(0, 0, 0), "dgCMatrix"),
           rawNorm      = as(matrix(0, 0, 0), "dgCMatrix"),
           coex         = as(matrix(0, 0, 0), "dgCMatrix"),
           cellsCoex    = as(matrix(0, 0, 0), "dgCMatrix"),
           nu           = vector(mode = "numeric"),
           lambda       = vector(mode = "numeric"),
           dispersion   = vector(mode = "numeric"),
           hkGenes      = vector(mode = "character"),
           metaDataset  = data.frame(),
           metaCells    = data.frame(),
           clustersCoex = vector(mode = "list")
         ),
         validity=function(object) {
           if ((!is_empty(object@rawNorm)) && (!identical(dim(object@rawNorm), dim(object@raw)))) {
             stop(paste("'rawNorm'[", nrow(object@rawNorm), "," , ncol(object@rawNorm),
                        "] must have the same sizes as ",
                        "'raw'[", nrow(object@rawNorm), "," , ncol(object@rawNorm),
                        "] when not empty."))
           }
           if ((!is_empty(object@nu)) && (length(object@nu) != ncol(object@raw))) {
             stop(paste("'nu'[", length(object@nu), "] must have size equal",
                        " to the number of columns of 'raw'[", dim(object@raw), "] when not empty."))
           }
           if ((!is_empty(object@lambda)) && (length(object@lambda) != nrow(object@raw))) {
             stop(paste("'lambda'[", length(object@lambda), "] must have size equal",
                        " to the number of rows of 'raw'[", dim(object@raw), "] when not empty."))
           }
           if ((!is_empty(object@dispersion) &&
                ((length(object@dispersion) + length(object@hkGenes)) != nrow(object@raw)))) {
             stop(paste("'dispersion'[", length(object@dispersion), "] size plus",
                        " 'hkGenes'[", length(object@hkGenes), "] size must be equal",
                        " to the number of rows of 'raw'[", dim(object@raw), "] when not empty."))
           }
           # TODO: check remaining slots
           return(TRUE)
         }
)

# constructor of the COTAN CLASS
#' @export 
COTAN <- function(raw = "ANY") {
  raw <- as(as.matrix(raw), "dgCMatrix")
  new("COTAN", raw = raw)
}


#' scCOTAN-class (for legacy usage)
#'
#' Define my COTAN structure
#' @slot raw ANY. To store the raw data matrix
#' @slot raw.norm ANY. To store the raw data matrix divided for the cell
#' efficiency estimated (nu)
#' @slot coex ANY. The coex matrix (sparce)
#' @slot nu vector.
#' @slot lambda vector.
#' @slot a vector.
#' @slot hk vector.
#' @slot meta data.frame.
#' @slot yes_yes ANY. Unused and deprecated
#' @slot clusters vector.
#' @slot cluster_data data.frame.
#'
#' @return the object class
#' @export
#' @examples
#'
#' data("ERCCraw")
#' obj <- new("scCOTAN", raw = data)
#'
setClass("scCOTAN",
         slots = c(
           raw = "ANY",
           raw.norm = "ANY",
           coex = "ANY",
           nu = "vector",
           lambda = "vector",
           a = "vector",
           hk = "vector",
           meta = "data.frame",
           yes_yes = "ANY",
           clusters = "vector",
           cluster_data = "data.frame"
         )
) -> scCOTAN



## Automatically convert an object from class "scCOTAN" into "COTAN"
setIs("scCOTAN",
      "COTAN",
      coerce = function(from) {
        if (!is_empty(from@yes_yes)) {
            warning("scCOTAN as COTAN: non-empty yes_yes member found: will be discarded",
                    call. = FALSE)
        }

        if (is_empty(from@raw)) {
          raw = as(matrix(0, 0, 0), "dgCMatrix")
        }
        else if(is(from@raw, "dgCMatrix")) {
          raw = from@raw
        }
        else {
          raw = as(as.matrix(from@raw), "dgCMatrix")
        }

        if (is_empty(from@raw.norm)) {
          rawNorm = as(matrix(0, 0, 0), "dgCMatrix")
        }
        else if(is(from@raw.norm, "dgCMatrix")) {
          rawNorm = from@raw.norm
        }
        else {
          rawNorm = as(as.matrix(from@raw.norm), "dgCMatrix")
        }
        
        if (!is_empty(from@clusters)) {
          metaCells <- data.frame(clusters = from@clusters,
                                  row.names = colnames(from@raw))
        }
        else {
          metaCells <- data.frame()
        }
        
        if (!is_empty(from@cluster_data)) {
          clustersCoex <- list(cluster_data = from@cluster_data)
        }
        else {
          clustersCoex <- vector(mode = "list")
        }

        new("COTAN",
            raw          = raw,
            rawNorm      = rawNorm,
            coex         = from@coex,
            cellsCoex    = as(matrix(0, 0, 0), "dgCMatrix"),
            nu           = from@nu,
            lambda       = from@lambda,
            dispersion   = from@a,
            hkGenes      = from@hk,
            metaDataset  = from@meta,
            metaCells    = metaCells,
            clustersCoex = clustersCoex )
      },
      # 'from' arg-name is convention: it is actually a destination!
      replace = function(from, value) {
        if(!is_empty(value@yes_yes)) {
          warning("scCOTAN<- as COTAN<-: non-empty yes_yes member found: will be discarded",
                  call. = FALSE)
        }

        if (is_empty(value@raw)) {
          raw = as(matrix(0, 0, 0), "dgCMatrix")
        }
        else if(is(value@raw, "dgCMatrix")) {
          raw = value@raw
        }
        else {
          raw = as(as.matrix(value@raw), "dgCMatrix")
        }

        if (is_empty(value@raw.norm)) {
          rawNorm = as(matrix(0, 0, 0), "dgCMatrix")
        }
        else if(is(value@raw.norm, "dgCMatrix")) {
          rawNorm = value@raw.norm
        }
        else {
          rawNorm = as(as.matrix(value@raw.norm), "dgCMatrix")
        }

        if (!is_empty(value@clusters)) {
          metaCells <- data.frame(clusters = value@clusters,
                                  row.names = colnames(value@raw))
        }
        else {
          metaCells <- data.frame()
        }

        if (!is_empty(value@cluster_data)) {
          clustersCoex <- list(cluster_data = value@cluster_data)
        }
        else {
          clustersCoex <- vector(mode = "list")
        }

        from@raw          <- raw
        from@rawNorm      <- rawNorm
        from@coex         <- value@coex
        from@cellCoex     <- as(matrix(0, 0, 0), "dgCMatrix")
        from@nu           <- value@nu
        from@lambda       <- value@lambda
        from@dispersion   <- value@a
        from@hkGenes      <- value@hk
        from@metaDataset  <- value@meta
        from@metaCells    <- metaCells
        from@clustersCoex <- clustersCoex
        from})

## Explicitly convert an object from class "COTAN" into "scCOTAN"
setAs("COTAN",
      "scCOTAN",
      function(from) {
        assertthat::assert_that(!is.null(from@metaCells),
                                !is.null(from@clustersCoex),
                                msg = "Unexpected COTAN null members")

        if (!is_empty(from@metaCells[['clusters']])) {
          clusters <- from@metaCells[['clusters']]
          names(clusters) <- rownames(from@metaCells)
        }
        else {
          clusters <- vector()
        }

        if (!is_empty(from@clustersCoex[['cluster_data']])) {
          cluster_data <- from@clustersCoex[['cluster_data']]
        }
        else {
          cluster_data <- data.frame()
        }

        new("scCOTAN",
            raw          = from@raw,
            raw.norm     = from@rawNorm,
            coex         = from@coex,
            nu           = from@nu,
            lambda       = from@lambda,
            a            = from@dispersion,
            hk           = from@hkGenes,
            meta         = from@metaDataset,
            clusters     = clusters,
            cluster_data = cluster_data )
      },
      # 'from' arg-name is convention: it is actually a destination!
      replace = function(from, value) {
        assertthat::assert_that(!is.null(from@metaCells),
                                !is.null(from@clustersCoex),
                                msg = "Unexpected scCOTAN null members")

        if (!is_empty(value@metaCells[['clusters']])) {
          clusters <- value@metaCells[['clusters']]
        }
        else {
          clusters <- vector()
        }

        if (!is_empty(value@clustersCoex[['cluster_data']])) {
          cluster_data <- value@clustersCoex[['cluster_data']]
        }
        else {
          cluster_data <- data.frame()
        }

        from@raw          <- value@raw
        from@raw.norm     <- value@rawNorm
        from@coex         <- value@coex
        from@nu           <- value@nu
        from@lambda       <- value@lambda
        from@a            <- value@dispersion
        from@hk           <- value@hkGenes
        from@meta         <- value@metaDataset
        from@clusters     <- clusters
        from@cluster_data <- cluster_data
        from})
