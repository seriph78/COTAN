#' Definition of COTAN class
#' @slot raw raw UMI count matrix 𝑛×𝑚 (gene number × cell number)
#' @slot coex correlation of COTAN between genes, 𝑛×𝑛
#' @slot nu vector that stores the estimated UDE, size 𝑚
#' @slot lambda vector to store the average for the gene expression, size 𝑛
#' @slot dispersion vector to store all
#' the negative binomial dispersion factors, size 𝑛.
#' @slot hKGenes house-keeping genes. It is a vector to store the name
#' of the genes with positive UMI count in every single cell of the sample
#' @slot metaDataset data.frame
#' @slot metaCells data.frame
#' @slot clustersCoex coex
#'
#' @importFrom rlang is_empty
#' @importClassesFrom Matrix dgCMatrix
#'
setClass(
  "COTAN",
  slots = c(
    raw          = "dgCMatrix",
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
  prototype = list(
    raw          = as(matrix(0, 0, 0), "dgCMatrix"),
    coex         = vector(mode = "list"),
    cellsCoex    = vector(mode = "list"),
    nu           = vector(mode = "numeric"),
    lambda       = vector(mode = "numeric"),
    dispersion   = vector(mode = "numeric"),
    hkGenes      = vector(mode = "character"),
    metaDataset  = data.frame(),
    metaCells    = data.frame(),
    clustersCoex = vector(mode = "list")
  ),
  validity = function(object) {
    if (!is_empty(object@raw) && is_empty(object@lambda)) {
      if (anyNA(object@raw)) {
        stop("Input raw data contains NA!")
      }
      if (isFALSE(all.equal(object@raw, round(object@raw), tolerance = 0))) {
        stop("Input raw data contains not integer numbers!")
      }
      if (any(object@raw < 0)) {
        stop("Input raw data must contain only non negative integers!")
      }
    }
    if (!is_empty(object@nu) && length(object@nu) != ncol(object@raw)) {
      stop(paste0("'nu'[", length(object@nu), "] must have size equal",
                  " to the number of columns [", ncol(object@raw),
                  "] of 'raw' when not empty."))
    }
    if (!is_empty(object@lambda) && length(object@lambda) != nrow(object@raw)) {
      stop(paste0("'lambda'[", length(object@lambda), "] must have size equal",
                  " to the number of rows [", nrow(object@raw),
                  "]  of 'raw' when not empty."))
    }
    # if (!is_empty(object@dispersion) && length(object@dispersion) != nrow(object@raw)) {
    #   stop(paste0("'dispersion'[", length(object@dispersion), "] must have size equal",
    #               " to the number of rows [", nrow(object@raw),
    #               "] of 'raw' when not empty."))
    # }
    # if (!is_empty(object@hkGenes) && length(object@hkGenes) != nrow(object@raw)) {
    #   stop(paste0("'hkGenes'[", length(object@hkGenes), "] must have size equal",
    #               " to the number of rows [", nrow(object@raw),
    #               "]  of 'raw' when not empty."))
    # }
    # metaDataset has no required fields as of now
    if (!is_empty(object@metaCells) && nrow(object@metaCells) != ncol(object@raw)) {
      stop(paste0("The number of rows [", nrow(object@metaCells),
                  "] of 'metaCells' must be the same",
                  " as the number of cols [", ncol(object@raw),
                  "]  of 'raw' when not empty."))
    }
    for (clusterName in names(object@clustersCoex)) {
      if (!clusterName %in% colnames(object@metaCells)) {
        stop(paste0("The cluster name '", clusterName, "' found in 'clustersCoex'",
                    " must be one of the column names of 'metaCells'"))
      }
      if (!isa((object@clustersCoex)[[clusterName]], "data.frame")) {
        stop(paste0("'clusterCoex' is supposedly compose of data.frames.",
                    " A '", class((object@clustersCoex)[[clusterName]]), "' was given instead" ))
      }
      if (nrow((object@clustersCoex)[[clusterName]]) != nrow(object@raw)) {
        stop(paste0("The number of rows [", nrow((object@clustersCoex)[[clusterName]]),
                    "] of the data.frames in 'clustersCoex' must be the same",
                    " as the number of rows [", nrow(object@raw),
                    "]  of 'raw'."))
      }
    }
    return(TRUE)
  }
) #end class COTAN


#' COTAN
#'
#' constructor of the class COTAN
#' @importFrom methods new
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#'
#' data("raw.dataset")
#' obj <- COTAN(raw = raw.dataset)
#'
#' @rdname COTAN
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
#' @slot yes_yes ANY. Unused and deprecated. Kept for backward compatibility
#' only
#' @slot clusters vector.
#' @slot cluster_data data.frame.
#'
#' @return the object class
#' @export
setClass(
  "scCOTAN",
  slots = c(
    raw = "ANY",
    raw.norm = "ANY",
    coex = "ANY",
    nu = "vector",
    lambda = "vector",
    a = "vector",
    hk = "vector",
    meta = "data.frame",
    yes_yes = "ANY", # Unused and deprecated: kept for backward compatibility only
    clusters = "vector",
    cluster_data = "data.frame"
  )
) -> scCOTAN

# Automatically convert an object from class "scCOTAN" into "COTAN"
#' @importFrom methods setIs
#' @importFrom rlang is_empty
#' @importClassesFrom Matrix dgCMatrix
setIs("scCOTAN",
      "COTAN",
      coerce = function(from) {
        if (!is_empty(from@yes_yes)) {
            warning(paste0("scCOTAN as COTAN: non-empty yes_yes member",
                           " found: will be discarded"),
                    call. = FALSE)
        }

        if (is_empty(from@raw)) {
          raw = as(matrix(0, 0, 0), "dgCMatrix")
        }
        else if(isa(from@raw, "dgCMatrix")) {
          raw = from@raw
        }
        else {
          raw = as(as.matrix(from@raw), "dgCMatrix")
        }

        if (!is_empty(from@clusters)) {
          metaCells <- data.frame(clusters = from@clusters,
                                  row.names = colnames(from@raw))
        }
        else if (!is_empty(from@cluster_data)) {
          # ensure some 'formal' validity
          metaCells <- data.frame(clusters = rep(NA, ncol(from@raw)),
                                  row.names = colnames(from@raw))
        }
        else {
          metaCells <- data.frame()
        }

        if (!is_empty(from@cluster_data)) {
          clustersCoex <- list(clusters = from@cluster_data)
        }
        else {
          clustersCoex <- vector(mode = "list")
        }

        new("COTAN",
            raw          = raw,
            coex         = from@coex,
            cellsCoex    = list(),
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
          warning(paste0("scCOTAN<- as COTAN<-: non-empty yes_yes member",
                         " found: will be discarded"),
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

        if (!is_empty(value@clusters)) {
          metaCells <- data.frame(clusters = value@clusters,
                                  row.names = colnames(value@raw))
        }
        else if (!is_empty(value@cluster_data)) {
          # ensure some 'formal' validity
          metaCells <- data.frame(clusters = rep(NA, ncol(value@raw)),
                                  row.names = colnames(value@raw))
        }
        else {
          metaCells <- data.frame()
        }

        if (!is_empty(value@cluster_data)) {
          clustersCoex <- list(clusters = value@cluster_data)
        }
        else {
          clustersCoex <- vector(mode = "list")
        }

        from@raw          <- raw
        from@coex         <- value@coex
        from@cellCoex     <- list()
        from@nu           <- value@nu
        from@lambda       <- value@lambda
        from@dispersion   <- value@a
        from@hkGenes      <- value@hk
        from@metaDataset  <- value@meta
        from@metaCells    <- metaCells
        from@clustersCoex <- clustersCoex
        from}
      ) # end setIs

# Explicitly convert an object from class "COTAN" into "scCOTAN"
#' @importFrom methods setAs
#' @importFrom rlang is_empty
#' @importClassesFrom Matrix dgCMatrix
setAs("COTAN",
      "scCOTAN",
      function(from) {
        if (!is_empty(from@clustersCoex)) {
          # pick last element as the most relevant!
          clustersName <- names(from@clustersCoex)[length(from@clustersCoex)]
          cluster_data <- from@clustersCoex[[clusterName]]
        }
        else {
          clustersName <- "clusters"
          cluster_data <- data.frame()
        }

        if (!is_empty(from@metaCells[[clustersName]])) {
          clusters <- from@metaCells[[clustersName]]
          names(clusters) <- rownames(from@metaCells)
        }
        else {
          # ensure non-empty vector
          clusters <- rep(NA, ncol(from@raw))
          names(clusters) <- colnames(from@raw)
        }

        new("scCOTAN",
            raw          = from@raw,
            raw.norm     = getNormalizedData(from),
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
        if (!is_empty(value@clustersCoex)) {
          # pick last element as the most relevant!
          clustersName <- names(value@clustersCoex)[length(value@clustersCoex)]
          cluster_data <- value@clustersCoex[[clusterName]]
        }
        else {
          clustersName <- "clusters"
          cluster_data <- data.frame()
        }

        if (!is_empty(value@metaCells[[clustersName]])) {
          clusters <- value@metaCells[[clustersName]]
          names(clusters) <- rownames(value@metaCells)
        }
        else {
          # ensure non-empty vector
          clusters <- rep(NA, ncol(value@raw))
          names(clusters) <- colnames(value@raw)
        }

        from@raw          <- value@raw
        from@raw.norm     <- getNormalizedData(value)
        from@coex         <- value@coex
        from@nu           <- value@nu
        from@lambda       <- value@lambda
        from@a            <- value@dispersion
        from@hk           <- value@hkGenes
        from@meta         <- value@metaDataset
        from@clusters     <- clusters
        from@cluster_data <- cluster_data
        from}
      ) # end setAs
