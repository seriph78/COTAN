#' @title emptySparseMatrix
#'
#' @description Useful to default initialize `COTAN` slots
#'
#' @returns an empty matrix of type `dgCMatrix`
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @noRd
#'
emptySparseMatrix <- function() {
  return(as(as(as(matrix(0.0, 0L, 0L), "dMatrix"),
               "generalMatrix"), "CsparseMatrix"))
}

#' @title emptySymmetricMatrix
#'
#' @description Useful to default initialize `COTAN` slots
#'
#' @returns an empty matrix of type `dspMatrix`
#'
#' @importClassesFrom Matrix dspMatrix
#'
#' @noRd
#'
emptySymmetricMatrix <- function() {
  return(as(as(as(matrix(0.0, 0L, 0L), "dMatrix"),
               "symmetricMatrix"), "packedMatrix"))
}

#---------- COTAN class --------------

#' Definition of the `COTAN` class
#'
#' @slot raw `dgCMatrix` - the raw UMI count matrix \eqn{n \times m} (gene
#'   number × cell number)
#' @slot genesCoex `dspMatrix` - the correlation of `COTAN` between genes,
#'   \eqn{n \times n}
#' @slot cellsCoex `dspMatrix` - the correlation of `COTAN` between cells,
#'   \eqn{m \times m}
#' @slot metaDataset `data.frame`
#' @slot metaCells `data.frame`
#' @slot clustersCoex a `list` of `COEX` `data.frames` for each clustering in
#'   the metaCells
#'
#' @importFrom rlang is_empty
#'
#' @importFrom methods validObject
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
    } else if (is_empty(object@metaGenes)) { #run the test only at the beginning
      if (anyNA(object@raw)) {
        stop("Input 'raw' data contains NA!")
      }
      if (isFALSE(all.equal(object@raw, round(object@raw), tolerance = 0.0))) {
        stop("Input 'raw' data contains non integer numbers.")
      }
      if (any(object@raw < 0.0)) {
        stop("Input 'raw' data must contain only non negative integers.")
      }
    }

    numGenes <- nrow(object@raw)
    numCells <- ncol(object@raw)

    if (length(rownames(object@raw)) != numGenes) {
      stop("Input 'raw' data must have row names")
    }
    if (length(colnames(object@raw)) != numCells) {
      stop("Input 'raw' data must have column names")
    }

    if (!is_empty(object@genesCoex)) {
      if (!isa(object@genesCoex, "dspMatrix")) {
        stop("'genesCoex' must be of type 'dspMatrix.")
      }
      if (nrow(object@genesCoex) != numGenes) {
        stop("'genesCoex' sizes should match the number of genes.")
      }
    }

    if (!is_empty(object@cellsCoex)) {
      if (!isa(object@cellsCoex, "dspMatrix")) {
        stop("'cellsCoex' must be of type 'dspMatrix.")
      }
      if (nrow(object@cellsCoex) != numCells) {
        stop("'cellsCoex' sizes should match the number of cells.")
      }
    }

    # metaDataset has no required fields as of now

    if (!is_empty(object@metaGenes) && nrow(object@metaGenes) != numGenes) {
      stop("The number of rows [", nrow(object@metaGenes), "] of 'metaGenes'",
           " must be the same as the number of rows [", numGenes, "] of 'raw'",
           " when not empty.")
    }

    if (!is_empty(object@metaCells) && nrow(object@metaCells) != numCells) {
      stop("The number of rows [", nrow(object@metaCells), "] of 'metaCells'",
           " must be the same as the number of cols [", numCells, "] of 'raw'",
           " when not empty.")
    }

    for (name in names(object@clustersCoex)) {
      if (!startsWith(name, "CL_")) {
        stop("The clusterization name '", name, "' does not start",
             " with 'CL_' as per conventions")
      }
      if (!name %in% colnames(object@metaCells)) {
        stop("The clusterization name '", name, "', found in 'clustersCoex',",
             " must be one of the column names of 'metaCells'.")
      }
      coexDF <- object@clustersCoex[[name]]
      if (!is_empty(coexDF)) {
        if (!isa(coexDF, "data.frame")) {
          stop("'clusterCoex' is supposedly composed of data.frames.",
               " A '", class(coexDF), "' was given instead.")
        }
        if (!identical(rownames(coexDF), rownames(object@raw))) {
          stop("The row-names of the non-empty data.frames in 'clustersCoex'",
               " must be the same as those of the raw data.")
        }
        if (!all(colnames(coexDF) %in% object@metaCells[[name]])) {
          badClustersTags <- colnames(coexDF)[!colnames(coexDF) %in%
                                                object@metaCells[[name]]]
          stop("The column names of the data.frames in 'clustersCoex' must be",
               " a subset of the tags in the clusterization.",
               " The data.frame for clusterization '", name, "' does not",
               " satisfy this condition with the names [",
               paste(badClustersTags, collapse = ","), "] vs cluster tags [",
               paste(unique(object@metaCells[[name]]), collapse = ","), "].")
        }
      }
    }

    for (name in colnames(object@metaCells)) {
      if (startsWith(name, "CL_")) {
        if (!inherits(object@metaCells[[name]], "factor")) {
          # ensure the clusters are factors
          object@metaCells[[name]] <- factor(object@metaCells[[name]])
        }
        if (!name %in% names(object@clustersCoex)) {
          stop("The clusterization name '", name, "' does not have",
               " an element in the 'clusterCoex' list")
        }
      } else if (startsWith(name, "COND_")) {
        if (!inherits(object@metaCells[[name]], "factor")) {
          # ensure the clusters are factors
          object@metaCells[[name]] <- factor(object@metaCells[[name]])
        }
      }
    }

    return(TRUE)
  }
) #end class COTAN


#'
#' @title `COTAN` shortcuts
#'
#' @description These functions create a [COTAN-class] object and/or also run
#'   all the necessary steps until the genes' `COEX` matrix is calculated.
#'
#' @name COTAN_ObjectCreation
NULL

#' @details Constructor of the class `COTAN`
#'
#' @param raw any object that can be converted to a matrix, but with row (genes)
#'   and column (cells) names
#'
#' @returns a `COTAN` object
#'
#' @importFrom methods new
#' @importFrom assertthat assert_that
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' obj <- COTAN(raw = test.dataset)
#'
#' @rdname COTAN_ObjectCreation
#'
COTAN <- function(raw = "ANY") {
  raw <- as(as(raw, "Matrix"), "sparseMatrix")

  assert_that(!is_empty(rownames(raw)), !is_empty(colnames(raw)),
              msg = "Inputs must have both row and column names!")

  new("COTAN", raw = raw)
}

#'
#' @title Handle legacy `scCOTAN`-class and related symmetric matrix <-> vector
#' conversions
#'
#' @description A class and some functions related to the `V1` version of the
#'   `COTAN` package
#'
#' @name COTAN_Legacy
NULL


#' @details Define the legacy `scCOTAN`-class
#'
#' @slot raw `ANY`. To store the raw data matrix
#' @slot raw.norm `ANY`. To store the raw data matrix divided for the cell
#'   efficiency estimated (nu)
#' @slot coex `ANY`. The coex matrix
#' @slot nu `vector`.
#' @slot lambda `vector`.
#' @slot a `vector`.
#' @slot hk `vector`.
#' @slot n_cells `numeric`.
#' @slot meta `data.frame`.
#' @slot yes_yes `ANY`. Unused and deprecated. Kept for backward compatibility
#'   only
#' @slot clusters `vector`.
#' @slot cluster_data `data.frame`.
#'
#' @returns a `scCOTAN` object
#'
#' @export
#'
#' @rdname COTAN_Legacy
#'
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
    n_cells = "numeric",
    meta = "data.frame",
    yes_yes = "ANY", # Unused and deprecated: for backward compatibility only
    clusters = "vector",
    cluster_data = "data.frame"
  )
) -> scCOTAN


#' @description Helper function to be shared by coerce() and replace()
#'
#' @param a `scCOTAN` object
#'
#' @returns a `list` with all non trivially converted slots of the equivalent
#'   `COTAN` class
#'
#' @importFrom rlang is_empty
#'
#' @importFrom stringr str_remove
#' @importFrom stringr fixed
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom Matrix dspMatrix
#'
#' @importFrom Matrix forceSymmetric
#' @importFrom Matrix pack
#'
#' @noRd
#'
getCOTANSlots <- function(from) {
  if (!is_empty(from@yes_yes)) {
    warning("scCOTAN as COTAN: non-empty yes_yes member found -",
            " Will be discarded!",
            call. = FALSE)
  }

  if (from@n_cells != ncol(from@raw)) {
    warning("scCOTAN as COTAN: n_cells member misaligned",
            " with raw matrix member - Will be ignored!",
            call. = FALSE)
  }

  if (is_empty(from@raw)) {
    raw <- emptySparseMatrix()
  } else if (isa(from@raw, "dgCMatrix")) {
    raw <- from@raw
  } else {
    raw <- as(as(as(as.matrix(from@raw), "dMatrix"),
                 "generalMatrix"), "CsparseMatrix")
  }

  genesCoex <- emptySymmetricMatrix()
  if (!is_empty(from@coex)) {
    if (isa(from@coex, "list")) {
      genesCoex <- vec2mat_rfast(from@coex)
    } else if (isa(from@coex, "dMatrix")) {
      genesCoex <- as.matrix(from@coex)
    } else if (isa(from@coex, "matrix")) {
      genesCoex <- from@coex
    } else {
      warning("scCOTAN as COTAN: unsupported coex type '", class(from@coex),
              "' - genes' coex will be discarded!", call. = FALSE)
      genesCoex <- emptySymmetricMatrix()
    }
    # now 'coex' is of type 'matrix'
    genesCoex <- pack(forceSymmetric(genesCoex))
  }

  cellsCoex <- emptySymmetricMatrix()

  metaGenes <- data.frame()

  if (!is_empty(from@hk)) {
    metaGenes <- setColumnInDF(metaGenes, rownames(from@raw) %in% from@hk,
                               "feGenes", rownames(from@raw))
  }

  if (!is_empty(from@lambda)) {
    metaGenes <- setColumnInDF(metaGenes, from@lambda,
                               "lambda", rownames(from@raw))
  }

  if (!is_empty(from@a)) {
    metaGenes <- setColumnInDF(metaGenes, from@a,
                               "dispersion", rownames(from@raw))
  }

  metaCells <- data.frame()

  if (!is_empty(from@nu)) {
    metaCells <- setColumnInDF(metaCells, from@nu,
                               "nu", colnames(from@raw))
  }

  clustersCoex <- vector(mode = "list")

  hasClusters <- !is_empty(from@clusters) && !all(is.na(from@clusters))
  if (hasClusters) {
    metaCells <- setColumnInDF(metaCells, factor(from@clusters),
                               "CL_clusters", colnames(from@raw))
  }

  {
    clusterData <- from@cluster_data

    if (!hasClusters && !is_empty(clusterData)) {
      warning("scCOTAN as COTAN: cannot have 'cluster_data' along",
              " empty 'clusters' - genes' coex will be discarded!",
              call. = FALSE)
      clusterData <- data.frame()
    }

    if (!is_empty(clusterData) &&
          !all(rownames(clusterData) %in% rownames(from@raw))) {
      warning("scCOTAN as COTAN: 'cluster_data' refers to unknown genes",
              " - clusters' coex will be discarded!", call. = FALSE)
      clusterData <- data.frame()
    }

    if (!is_empty(clusterData)) {
      if (!all(colnames(clusterData) %in% from@clusters)) {
        # It might be possible that the column names have the old extra
        # prefix 'cl.'. It will be remove in such cases.
        colnames(clusterData) <- str_remove(colnames(clusterData),
                                            pattern = fixed("cl."))
      }

      if (!all(colnames(clusterData) %in% from@clusters)) {
        warning("scCOTAN as COTAN: 'cluster_data' refers to unknown",
                " clusters - clusters' coex will be discarded!",
                call. = FALSE)
        clusterData <- data.frame()
      }
    }

    if (!is_empty(clusterData) &&
          !all(rownames(from@raw) %in% union(rownames(clusterData), from@hk))) {
      warning("scCOTAN as COTAN: 'cluster_data' has no information",
              " on some genes that are not fully-expressed",
              " - clusters' coex will be discarded!", call. = FALSE)
      clusterData <- data.frame()
    }

    if (!is_empty(clusterData)) {
      missingGenes <- setdiff(rownames(from@raw), rownames(clusterData))
      if (!is_empty(missingGenes)) {
        # missing genes are fully-expressed thus have all zero coex!
        missingData <- matrix(0.0, nrow = length(missingGenes),
                              ncol = ncol(clusterData),
                              dimnames = c(missingGenes, colnames(clusterData)))
        clusterData <- rbind(clusterData, as.data.frame(missingData))
        # reorder the rows to match original data
        clusterData <- clusterData[rownames(from@raw)]
      }
      clustersCoex <- list(CL_clusters = clusterData)
    } else if (hasClusters) {
      # keep aligned
      clustersCoex <- list("CL_clusters" = data.frame())
    }
  }

  return(list(raw, genesCoex, cellsCoex, metaGenes, metaCells, clustersCoex))
}

#' @details Automatically converts an object from class `scCOTAN` into `COTAN`
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom methods setIs
#'
#' @name scCotan_coerce_to_COTAN
#' @rdname COTAN_Legacy
#'
setIs(
  "scCOTAN",
  "COTAN",
  coerce = function(from) {
    c(raw, genesCoex, cellsCoex,
      metaGenes, metaCells, clustersCoex) %<-% getCOTANSlots(from)

      new("COTAN",
          raw          = raw,
          genesCoex    = genesCoex,
          cellsCoex    = cellsCoex,
          metaDataset  = from@meta,
          metaGenes    = metaGenes,
          metaCells    = metaCells,
          clustersCoex = clustersCoex)
  },
  # 'from' arg-name is convention: it is actually a destination!
  replace = function(from, value) {
    c(raw, genesCoex, cellsCoex,
      metaGenes, metaCells, clustersCoex) %<-% getCOTANSlots(value)

    # Update the existing scCOTAN object in place
    from@raw          <- raw
    from@genesCoex    <- genesCoex
    from@cellsCoex    <- cellsCoex
    from@metaDataset  <- value@meta
    from@metaGenes    <- metaGenes
    from@metaCells    <- metaCells
    from@clustersCoex <- clustersCoex

    # Return the modified object
    return(from)
  })


#'
#' @description Helper function to be shared by coerce() and replace()
#'
#' @param a `COTAN` object
#'
#' @returns a `list` with all non trivially converted slots of the equivalent
#'   `scCOTAN` class
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom Matrix dspMatrix
#'
#' @importFrom Matrix forceSymmetric
#' @importFrom Matrix pack
#'
#' @noRd
#'
getScCOTANSlots <- function(from) {
  lambda <- vector(mode = "numeric")
  if (!is_empty(from@metaGenes[["lambda"]])) {
    lambda <- set_names(from@metaGenes[["lambda"]], rownames(from@metaGenes))
  }

  a <- vector(mode = "numeric")
  if (!is_empty(from@metaGenes[["dispersion"]])) {
    a <- set_names(from@metaGenes[["dispersion"]], rownames(from@metaGenes))
  }

  nu <- vector(mode = "numeric")
  if (!is_empty(from@metaCells[["nu"]])) {
    nu <- set_names(from@metaCells[["nu"]], rownames(from@metaCells))
  }

  hk <- vector(mode = "character")
  if (!is_empty(from@metaGenes[["feGenes"]])) {
    hk <- rownames(from@raw)[from@metaGenes[["feGenes"]]]
  }

  rawNorm <- emptySparseMatrix()
  if (!is_empty(nu)) {
    rawNorm <- t(t(from@raw) * (1.0 / nu))
  }

  if (!is_empty(from@clustersCoex)) {
    if (length(from@clustersCoex) > 1L) {
      warning("COTAN as scCOTAN: more than one clusterization found:",
              "picking up the last one as likely to be the most relevant")
    }
    clName <- names(from@clustersCoex)[length(from@clustersCoex)]
    clusterData <- from@clustersCoex[[clName]]
  } else {
    clName <- "CL_clusters"
    clusterData <- data.frame()
  }

  if (!is_empty(from@metaCells[[clName]])) {
    clusters <- from@metaCells[[clName]]
    clusters <- set_names(levels(clusters)[clusters], rownames(from@metaCells))
  } else {
    # ensure non-empty vector
    clusters <- set_names(rep(NA, ncol(from@raw)), colnames(from@raw))
  }

  if (!is_empty(from@cellsCoex)) {
    warning("COTAN as scCOTAN: 'cellsCoex' is not empty:",
            " will be lost in conversion")
  }

  return(list(rawNorm, nu, lambda, a, hk, clusters, clusterData))
}

#' @details Explicitly converts an object from class `COTAN` into `scCOTAN`
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom methods as
#' @importFrom methods setAs
#'
#' @name COTAN_coerce_to_scCOTAN
#' @rdname COTAN_Legacy
#'
setAs(
  "COTAN",
  "scCOTAN",
  function(from) {
    c(rawNorm, nu, lambda, a,
      hk, clusters, clusterData) %<-% getScCOTANSlots(from)

      new("scCOTAN",
          raw          = from@raw,
          raw.norm     = rawNorm,
          coex         = from@genesCoex,
          nu           = nu,
          lambda       = lambda,
          a            = a,
          hk           = hk,
          n_cells      = ncol(from@raw),
          meta         = from@metaDataset,
          clusters     = clusters,
          cluster_data = clusterData)
  },
  # 'from' arg-name is convention: it is actually a destination!
  replace = function(from, value) {
    # Extract the slots needed for updating the scCOTAN object
    c(rawNorm, nu, lambda, a,
      hk, clusters, clusterData) %<-% getScCOTANSlots(value)

    # Update the existing scCOTAN object in place
    from@raw          <- value@raw
    from@raw.norm     <- rawNorm
    from@coex         <- value@genesCoex
    from@nu           <- nu
    from@lambda       <- lambda
    from@a            <- a
    from@hk           <- hk
    from@n_cells      <- ncol(value@raw)
    from@meta         <- value@metaDataset
    from@clusters     <- clusters
    from@cluster_data <- clusterData

    # Return the modified object
    return(from)
  })


#---------- Transcript Uniformity Check --------------

#' @title Definition of the **Transcript Uniformity Check** classes
#'
#' @description A hierarchy of classes to specify the method for checking
#'   whether a **cluster** has the *Uniform Transcript* property. It also
#'   doubles as result object.
#'
#' @name UT_Check
NULL

#' @details `BaseUniformityCheck` is the base class of the check methods
#'
#' @slot isUniform Logical. Output. The result of the check
#' @slot clusterSize Integer. Output. The number of cells in the cluster
#'
#' @importFrom rlang is_logical
#' @importFrom rlang is_integer
#'
#' @rdname UT_Check
#'
setClass(
  "BaseUniformityCheck",
  slots = c(
    isUniform = "logical",
    clusterSize = "integer"
  ),
  prototype = list(
    isUniform   = FALSE,
    clusterSize = 0L
  ),
  validity = function(object) {
    if (!is_logical(object@isUniform, n = 1)) {
      stop("Given 'isUniform' data must be a Boolean")
    }
    if (!is_integer(object@clusterSize, n = 1) || object@clusterSize < 0) {
      stop("Given 'clusterSize' data must be a positive number")
    }
    return(TRUE)
  }
) # end Base Uniformity Check


#' @details `SimpleGDIUniformityCheck` represents the simplified (and legacy)
#'   mechanism to determine whether a cluster has the *Uniform Transcript*
#'   property
#'
#'   The method is based on checking whether the fraction of the genes' `GDI`
#'   below the given *threshold* is less than the given *ratio*
#'
#' @slot GDIThreshold Numeric. (Simple check only) The level of `GDI` above
#'   which the **cluster** is deemed not uniform (default is 1.43)
#' @slot ratioAboveThreshold Numeric. (Simple check only) The quantile of the
#'   `GDI` distribution used to select the gene whose `GDI` must be below the
#'   threshold (default is 0.01)
#' @slot fractionAboveThreshold Numeric. Output. The fraction of genes whose
#'   `GDI` is above the threshold
#' @slot quantileAtRatio Numeric. Output. The `GDI` distribution quantile
#'   corresponding at the given ratio
#'
#' @importFrom rlang is_double
#'
#' @rdname UT_Check
#'
setClass(
  "SimpleGDIUniformityCheck",
  contains = "BaseUniformityCheck",
  slots = c(
    GDIThreshold           = "numeric",
    ratioAboveThreshold    = "numeric",
    fractionAboveThreshold = "numeric",
    quantileAtRatio        = "numeric"
  ),
  prototype = list(
    GDIThreshold           = 1.43,
    ratioAboveThreshold    = 0.01,
    fractionAboveThreshold = NaN,
    quantileAtRatio        = NaN
  ),
  validity = function(object) {
    if (!is_double(object@GDIThreshold, n = 1, finite = TRUE) ||
        object@GDIThreshold <= 1.0) {
      stop("Input 'GDIThreshold' data must be a finite number above 1.0")
    }
    if (!is_double(object@ratioAboveThreshold, n = 1, finite = TRUE) ||
        (object@ratioAboveThreshold <= 0.0 ||
         object@ratioAboveThreshold >= 1.0)) {
      stop("Input 'ratioAboveThreshold' data must be a finite",
           " number between 0.0 and 1.0")
    }
    if (!is_double(object@fractionAboveThreshold, n = 1, finite = FALSE) ||
        (is.finite(object@fractionAboveThreshold) &&
         (object@fractionAboveThreshold <= 0.0 ||
          object@fractionAboveThreshold >= 1.0))) {
      stop("Given 'fractionAboveThreshold' data must be NaN or a finite",
           " number between 0.0 and 1.0")
    }
    if (!is_double(object@quantileAtRatio, n = 1, finite = FALSE) ||
        (is.finite(object@fractionAboveThreshold) &&
         object@quantileAtRatio <= 0.0)) {
      stop("Given 'quantileAtRatio' data must be NaN or a finite",
           " positive number")
    }
    return(TRUE)
  }
) # end class SimpleGDIUniformityCheck


#' @details `AdvancedGDIUniformityCheck` represents the more precise and
#'   advanced mechanism to determine whether a cluster has the *Uniform
#'   Transcript* property
#'
#'   The method is based on checking the genes' `GDI` against three *thresholds*
#'   A *cluster* is deemed uniform if the ...
#'
#' @slot lowCheckThreshold Numeric. (Advanced check only) The level of `GDI`
#'   above which the
#'   **cluster** is deemed not uniform (default is 1.297)
#' @slot lowCheckQuantile Numeric. (Advanced check only) The quantile of the
#'   `GDI` distribution used to select the gene whose `GDI` must be below the
#'   low threshold (default is 95%)
#' @slot lowCheckRank Integer. (Advanced check only) The rank of the `GDI`
#'   distribution used to select the gene whose `GDI` must be below the low
#'   threshold (default is `NULL`)
#' @slot middleCheckThreshold Numeric. (Advanced check only) The level of `GDI`
#'   below which the
#'   **cluster** is deemed not uniform
#' @slot middleCheckQuantile Numeric. (Advanced check only) The quantile of the
#'   `GDI` distribution used to select the gene whose `GDI` must be below the
#'   middle threshold (default is 98%)
#' @slot middleCheckRank Integer. (Advanced check only) The rank of the `GDI`
#'   distribution used to select the gene whose `GDI` must be below the middle
#'   threshold (default is `NULL`)
#' @slot highCheckThreshold Numeric. (Advanced check only) The level of `GDI`
#'   above which the
#'   **cluster** is deemed not uniform
#' @slot highCheckQuantile Numeric. (Advanced check only) The quantile of the
#'   `GDI` distribution used to select the gene whose `GDI` must be below the
#'   high threshold (default is `NULL`)
#' @slot highCheckRank Integer. (Advanced check only) The rank of the `GDI`
#'   distribution used to select the gene whose `GDI` must be below the high
#'   threshold (default is 2)
#'
#' @importFrom rlang is_integer
#' @importFrom rlang is_double
#'
#' @rdname UT_Check
#'
setClass(
  "AdvancedGDIUniformityCheck",
  contains = "BaseUniformityCheck",
  slots = c(
    lowCheckThreshold    = "numeric",
    lowCheckQuantile     = "numeric",
    lowCheckRank         = "integer",
    middleCheckThreshold = "numeric",
    middleCheckQuantile  = "numeric",
    middleCheckRank      = "integer",
    highCheckThreshold   = "numeric",
    highCheckQuantile    = "numeric",
    highCheckRank        = "integer"
  ),
  prototype = list(
    lowCheckThreshold    = 1.297,
    lowCheckQuantile     = 0.95,
    lowCheckRank         = 0L,
    middleCheckThreshold = 1.400,
    middleCheckQuantile  = NaN,
    middleCheckRank      = 2L,
    highCheckThreshold   = 1.307,
    highCheckQuantile    = 0.98,
    highCheckRank        = 0L
  ),
  validity = function(object) {
    allOK <- function(threshold, quantile, rank) {
      if (!is_double(threshold, n = 1, finite = TRUE) || threshold <= 1.0) {
        stop("Input check `threshold` must be a finite number above 1.0")
      }
      if (!is_double(quantile, n = 1)) {
        stop("Input check `quantile` must be a double")
      }
      if (!is_integer(rank, n = 1)) {
        stop("Input check `rank` must be an integer")
      }
      if (is.nan(quantile) == (rank == 0L)) {
        stop("Input check `quantile` and `rank` are mutually",
             " exclusive: only one is allowed")
      }

      if (!is.nan(quantile)) {
        if (quantile <= 0.0 || quantile >= 1.0) {
          stop("Input 'quantile' must be a number between 0.0 and 1.0")
        }
      } else {
        if(rank < 0L) {
          stop("Input 'rank' must be a positive integer")
        }
      }
      return(TRUE)
    }
    return(
      allOK(object@lowCheckThreshold,
            object@lowCheckQuantile,
            object@lowCheckRank) &&
      allOK(object@middleCheckThreshold,
            object@middleCheckQuantile,
            object@middleCheckRank) &&
      allOK(object@highCheckThreshold,
            object@highCheckQuantile,
            object@highCheckRan))
  }
) # end class AdvancedGDIUniformityCheck
