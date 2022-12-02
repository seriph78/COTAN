
#' @importClassesFrom Matrix dgCMatrix
#'
emptySparseMatrix <- function() {
  return( as(as(as(matrix(0, 0, 0), "dMatrix"), "generalMatrix"), "CsparseMatrix") )
}

#' @importClassesFrom Matrix dspMatrix
#'
emptySymmetricMatrix <- function() {
  return( as(as(as(matrix(0, 0, 0), "dMatrix"), "symmetricMatrix"), "packedMatrix") )
}


#' Definition of COTAN class
#'
#' @slot raw raw UMI count matrix 𝑛×𝑚 (gene number × cell number)
#' @slot genesCoex correlation of COTAN between genes, 𝑛×𝑛
#' @slot cellsCoex correlation of COTAN between cells, 𝑚×𝑚
#' @slot metaDataset data.frame
#' @slot metaCells data.frame
#' @slot clustersCoex a list of coex data.frames for each clustering in the metaCells
#'
#' @importFrom rlang is_empty
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom Matrix dspMatrix
#'
setClass(
  "COTAN",
  slots = c(
    raw          = "dgCMatrix",
    genesCoex    = "dspMatrix",
    cellsCoex    = "dspMatrix",
    metaDataset  = "data.frame",
    metaGenes    = "data.frame",
    metaCells    = "data.frame",
    clustersCoex = "list" # of data.frames
  ),
  prototype = list(
    raw          = emptySparseMatrix(),
    genesCoex    = emptySymmetricMatrix(),
    cellsCoex    = emptySymmetricMatrix(),
    metaDataset  = data.frame(),
    metaGenes    = data.frame(),
    metaCells    = data.frame(),
    clustersCoex = vector(mode = "list")
  ),
  validity = function(object) {
    if (is_empty(object@raw)) {
      stop("'raw' is empty")
    }
    else if (is_empty(object@metaGenes)) { #run the test only at the beginning
      if (anyNA(object@raw)) {
        stop("Input 'raw' data contains NA!")
      }
      if (isFALSE(all.equal(object@raw, round(object@raw), tolerance = 0))) {
        stop("Input 'raw' data contains non integer numbers.")
      }
      if (any(object@raw < 0)) {
        stop("Input 'raw' data must contain only non negative integers.")
      }
    }
    numGenes <- nrow(object@raw)
    numCells <- ncol(object@raw)
    if (!is_empty(object@genesCoex)) {
      if (!isa(object@genesCoex, "dspMatrix")) {
        stop("'genesCoex' must be of type 'dspMatrix.")
      }
      if (nrow(object@genesCoex) != numGenes) {
        stop("'genesCoex' sizes should match the number of genes.")
      }
    }
    if (!is_empty(object@cellsCoex)) {
      if (!isa(object@genesCoex, "dspMatrix")) {
        stop("'genesCoex' must be of type 'dspMatrix.")
      }
      if (nrow(object@genesCoex) != numGenes) {
        stop("'genesCoex' sizes should match the number of genes.")
      }
    }
    # metaDataset has no required fields as of now
    if (!is_empty(object@metaGenes) && nrow(object@metaGenes) != numGenes) {
      stop(paste0("The number of rows [", nrow(object@metaGenes),
                  "] of 'metaGenes' must be the same",
                  " as the number of rows [", numGenes,
                  "]  of 'raw' when not empty."))
    }
    if (!is_empty(object@metaCells) && nrow(object@metaCells) != numCells) {
      stop(paste0("The number of rows [", nrow(object@metaCells),
                  "] of 'metaCells' must be the same",
                  " as the number of cols [", numCells,
                  "]  of 'raw' when not empty."))
    }
    for (name in names(object@clustersCoex)) {
      if (substring(name, 1, 3) != "CL_") {
        stop(paste0("The clusterization name '", name, "' does not start",
                    " with 'CL_' as per conventions"))
      }
      if (!name %in% colnames(object@metaCells)) {
        stop(paste0("The clusterization name '", name, "', found in 'clustersCoex',",
                    " must be one of the column names of 'metaCells'."))
      }
      coexDF <- object@clustersCoex[[name]]
      if (!isa(coexDF, "data.frame")) {
        stop(paste0("'clusterCoex' is supposedly composed of data.frames.",
                    " A '", class(coexDF), "' was given instead." ))
      }
      if (!is_empty(coexDF)) {
        if (nrow(coexDF) != numGenes) {
          stop(paste0("The number of rows [", nrow(object@clustersCoex[[name]]),
                      "] of the non-empty data.frames in 'clustersCoex' must be the same",
                      " as the number of rows [", numGenes, "]  of 'raw'."))
        }
        # TODO: fix merge_cell.clusters that causes the following to fail
        # if (!all(colnames(coexDF) %in% object@metaCells[[name]])) {
        #   badClustersNames <- colnames(coexDF)[!colnames(coexDF) %in% object@metaCells[[name]]]
        #   stop(paste0("The column names of the data.frames in 'clustersCoex' must",
        #               " be a subset of the names in the clusterization.",
        #               " The data.frame for clusterization '", name, "' does not",
        #               " satisfy this condition with the names [",
        #               paste(badClustersNames, collapse = ","), "] vs cluester names [",
        #               paste(unique(object@metaCells[[name]]), collapse = ","), "]."))
        # }
      }
    }
    for (name in colnames(object@metaCells)) {
      if (substring(name, 1, 3) != "CL_") {
        next # not a clusterization name
      }
      if (!name %in% names(object@clustersCoex)) {
        stop(paste0("The clusterization name '", name, "' does not have",
                    " an element in the 'clusterCoex' list"))
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
  raw <- as(as(as(as.matrix(raw), "dMatrix"), "generalMatrix"), "CsparseMatrix")
  new("COTAN", raw = raw)
}


#' scCOTAN-class (for legacy usage)
#'
#' Define my COTAN structure
#' @slot raw ANY. To store the raw data matrix
#' @slot raw.norm ANY. To store the raw data matrix divided for the cell
#' efficiency estimated (nu)
#' @slot coex ANY. The coex matrix
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
#' @importClassesFrom Matrix dspMatrix
#' @importFrom Matrix forceSymmetric
#' @importFrom Matrix pack
#'
setIs("scCOTAN",
      "COTAN",
      coerce = function(from) {
        if (!is_empty(from@yes_yes)) {
            warning(paste0("scCOTAN as COTAN: non-empty yes_yes member",
                           " found: will be discarded"),
                    call. = FALSE)
        }

        if (is_empty(from@raw)) {
          raw = emptySparseMatrix()
        }
        else if(isa(from@raw, "dgCMatrix")) {
          raw = from@raw
        }
        else {
          raw = as(as(as(as.matrix(from@raw), "dMatrix"), "generalMatrix"), "CsparseMatrix")
        }

        genesCoex <- emptySymmetricMatrix()
        if (!is_empty(from@coex)) {
          if (isa(from@coex, "list")) {
            genesCoex <- vec2mat_rfast(from@coex[["values"]])
          }
          else if(isa(from@coex, "dMatrix")) {
            genesCoex <- as.matrix(from@coex)
          }
          else if(isa(from@coex, "matrix")) {
            genesCoex <- from@coex
          }
          else {
            stop(paste0("Unsupported coex type: ", class(from@coex)))
          }
          # now 'coex' is of type 'matrix'
          genesCoex <- pack(forceSymmetric(genesCoex))
        }

        cellsCoex <- emptySymmetricMatrix()

        metaGenes <- data.frame()
        if (!is_empty(from@hk)) {
          metaGenes <- setColumnInDF(metaGenes, rownames(from@raw) %in% from@hk,
                                     "hkGenes", rownames(from@raw))
        }
        if (!is_empty(from@lambda)) {
          metaGenes <- setColumnInDF(metaGenes, from@lambda,
                                     "lambda", rownames(from@raw))
        }
        if (!is_empty(from@a)) {
          metaGenes <- setColumnInDF(metaGenes, from@a,
                                     "dispersion", rownames(from@raw))
        }

        if (is_empty(from@clusters) && !is_empty(from@cluster_data)) {
          stop("Cannot convert to 'COTAN' when 'cluster_data' is not empty and 'clusters' is empty")
        }

        metaCells <- data.frame()
        if (!is_empty(from@clusters)) {
          metaCells <- setColumnInDF(metaCells, from@clusters,
                                     "CL_clusters", colnames(from@raw))
        }
        if (!is_empty(from@nu)) {
          metaCells <- setColumnInDF(metaCells, from@nu,
                                     "nu", colnames(from@raw))
        }

        if (!is_empty(from@cluster_data)) {
          clustersCoex <- list(CL_clusters = from@cluster_data)
        }
        else if (!is_empty(from@clusters)) {
          clustersCoex <- list("CL_clusters" = data.frame())
        }
        else {
          clustersCoex <- vector(mode = "list")
        }

        new("COTAN",
            raw          = raw,
            genesCoex    = genesCoex,
            cellsCoex    = cellsCoex,
            metaDataset  = from@meta,
            metaGenes    = metaGenes,
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
          raw = emptySparseMatrix()
        }
        else if(isa(value@raw, "dgCMatrix")) {
          raw = value@raw
        }
        else {
          raw = as(as(as(as.matrix(value@raw), "dMatrix"), "generalMatrix"), "CsparseMatrix")
        }

        genesCoex <- emptySymmetricMatrix()
        if (!is_empty(value@coex)) {
          if (isa(value@coex, "list")) {
            genesCoex <- vec2mat_rfast(value@coex[["values"]])
          }
          else if(isa(value@coex, "dMatrix")) {
            genesCoex <- as.matrix(value@coex)
          }
          else if(isa(value@coex, "matrix")) {
            genesCoex <- value@coex
          }
          else {
            stop(paste0("Unsupported coex type: ", class(value@coex)))
          }
          # now 'genesCoex' is of type 'matrix'
          genesCoex <- pack(forceSymmetric(genesCoex))
        }

        cellsCoex <- emptySymmetricMatrix()

        metaGenes <- data.frame()
        if (!is_empty(value@hk)) {
          metaGenes <- setColumnInDF(metaGenes, rownames(value@raw) %in% value@hk,
                                     "hkGenes", rownames(value@raw))
        }
        if (!is_empty(value@lambda)) {
          metaGenes <- setColumnInDF(metaGenes, value@lambda,
                                     "lambda", rownames(value@raw))
        }
        if (!is_empty(value@dispersion)) {
          metaGenes <- setColumnInDF(metaGenes, value@dipsersion,
                                     "dispersion", rownames(value@raw))
        }

        if (is_empty(value@clusters) && !is_empty(value@cluster_data)) {
          stop("Cannot convert to 'COTAN' when 'cluster_data' is not empty and 'clusters' is empty")
        }

        metaCells <- data.frame()
        if (!is_empty(value@clusters)) {
          metaCells <- setColumnInDF(metaCells, value@clusters,
                                     "CL_clusters", colnames(value@raw))
        }
        if (!is_empty(value@nu)) {
          metaCells <- setColumnInDF(metaCells, value@nu,
                                     "nu", colnames(value@raw))
        }

        if (!is_empty(value@cluster_data)) {
          clustersCoex <- list(CL_clusters = value@cluster_data)
        }
        else if (!is_empty(value@clusters)) {
          clustersCoex <- list("CL_clusters" = data.frame())
        }
        else {
          clustersCoex <- vector(mode = "list")
        }

        from@raw          <- raw
        from@genesCoex    <- genesCoex
        from@cellsCoex    <- cellsCoex
        from@metaDataset  <- value@meta
        from@metaGenes    <- metaGenes
        from@metaCells    <- metaCells
        from@clustersCoex <- clustersCoex
        from}
      ) # end setIs


# Explicitly convert an object from class "COTAN" into "scCOTAN"

#' @importFrom methods setAs
#' @importFrom rlang is_empty
#'
setAs("COTAN",
      "scCOTAN",
      function(from) {
        lambda <- vector(mode = "numeric")
        if (!is_empty(from@metaGenes[["lambda"]])) {
          lambda <- from@metaGenes[["lambda"]]
          names(lambda) <- rownames(from@metaGenes)
        }

        a <- vector(mode = "numeric")
        if (!is_empty(from@metaGenes[["dispersion"]])) {
          a <- from@metaGenes[["dispersion"]]
          names(a) <- rownames(from@metaGenes)
        }

        nu <- vector(mode = "numeric")
        if (!is_empty(from@metaCells[["nu"]])) {
          nu <- from@metaCells[["nu"]]
          names(nu) <- rownames(from@metaCells)
        }

        hk <- vector(mode = "character")
        if (!is_empty(from@metaGenes[["hkGenes"]])) {
          hk <- rownames(from@raw)[from@metaGenes[["hkGenes"]] != 0]
        }

        raw.norm <- emptySparseMatrix()
        if (!is_empty(nu)) {
          raw.norm <- t(t(from@raw) * (1/nu))
        }

        if (!is_empty(from@clustersCoex)) {
          if (length(from@clustersCoex) > 1) {
            warning(paste0("The 'COTAN' object contains more than one clusterization:",
                           "picking up the last one as likely to be the most relevant"))
          }
          clusterizationName <- names(from@clustersCoex)[length(from@clustersCoex)]
          cluster_data <- from@clustersCoex[[clusterizationName]]
        } else {
          clusterizationName <- "CL_clusters"
          cluster_data <- data.frame()
        }

        if (!is_empty(from@metaCells[[clusterizationName]])) {
          clusters <- from@metaCells[[clusterizationName]]
          names(clusters) <- rownames(from@metaCells)
        } else {
          # ensure non-empty vector
          clusters <- rep(NA, ncol(from@raw))
          names(clusters) <- colnames(from@raw)
        }

        if (!is_empty(from@cellsCoex)) {
          warning("'cellsCoex' is not empty: will be lost in conversion to 'scCOTAN'")
        }

        new("scCOTAN",
            raw          = from@raw,
            raw.norm     = raw.norm,
            coex         = from@genesCoex,
            nu           = nu,
            lambda       = lambda,
            a            = a,
            hk           = hk,
            meta         = from@metaDataset,
            clusters     = clusters,
            cluster_data = cluster_data)
      },
      # 'from' arg-name is convention: it is actually a destination!
      replace = function(from, value) {
        lamda <- vector(mode = "numeric")
        if (!is_empty(value@metaGenes[["lambda"]])) {
          lambda <- value@metaGenes[["lambda"]]
          names(lambda) <- rownames(value@metaGenes)
        }

        a <- vector(mode = "numeric")
        if (!is_empty(value@metaGenes[["dispersion"]])) {
          a <- value@metaGenes[["dispersion"]]
          names(a) <- rownames(value@metaGenes)
        }

        nu <- vector(mode = "numeric")
        if (!is_empty(value@metaCells[["nu"]])) {
          nu <- value@metaCells[["nu"]]
          names(nu) <- rownames(value@metaCells)
        }

        hk <- vector(mode = "character")
        if (!is_empty(value@metaGenes[["hkGenes"]])) {
          hk <- rownames(value@raw)[value@metaGenes[["hkGenes"]] != 0]
        }

        raw.norm <- emptySparseMatrix()
        if (!is_empty(nu)) {
          raw.norm <- t(t(value@raw) * (1/nu))
        }

        if (!is_empty(value@clustersCoex)) {
          if (length(value@clustersCoex) > 1) {
            warning(paste0("The 'COTAN' object contains more than one clusterization:",
                           "picking up the last one as likely to be the most relevant"))
          }
          clusterizationName <- names(value@clustersCoex)[length(value@clustersCoex)]
          cluster_data <- value@clustersCoex[[clusterizationName]]
        }
        else {
          clusterizationName <- "CL_clusters"
          cluster_data <- data.frame()
        }

        if (!is_empty(value@metaCells[[clusterizationName]])) {
          clusters <- value@metaCells[[clusterizationName]]
          names(clusters) <- rownames(value@metaCells)
        }
        else {
          # ensure non-empty vector
          clusters <- rep(NA, ncol(value@raw))
          names(clusters) <- colnames(value@raw)
        }

        if (!is_empty(value@cellsCoex)) {
          warning("'cellsCoex' is not empty: will be lost in conversion to 'scCOTAN'")
        }

        from@raw          <- value@raw
        from@raw.norm     <- raw.norm
        from@coex         <- value@genesCoex
        from@nu           <- nu
        from@lambda       <- lambda
        from@a            <- a
        from@hk           <- hk
        from@meta         <- value@metaDataset
        from@clusters     <- clusters
        from@cluster_data <- cluster_data
        from}
      ) # end setAs
