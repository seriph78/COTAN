# ----------------- execution options --------------------

#' @title Execution options for heavy COTAN computations
#'
#' @description A small parameter object bundling execution-related controls
#' such as multi-core usage, `torch` device selection, and solver batching
#' options. Some fields are only used by specific functions.
#'
#' @slot cores Integer scalar. Requested number of CPU cores.
#' @slot optimizeForSpeed Logical scalar. Whether to try `torch` acceleration.
#' @slot deviceStr Character scalar. Requested `torch` device string.
#' @slot chunkSize Integer scalar. Requested solver batch size.
#'
#' @name ExecutionOptions-class
#'
#' @exportClass ExecutionOptions
#'
#' @rdname ExecutionOptions
#'
setClass(
  "ExecutionOptions",
  slots = c(
    cores = "integer",
    optimizeForSpeed = "logical",
    deviceStr = "character",
    chunkSize = "integer"
  ),
  prototype = list(
    cores = 1L,
    optimizeForSpeed = TRUE,
    deviceStr = "cuda",
    chunkSize = 1024L
  ),
  validity = function(object) {
    if (length(object@cores) != 1L || is.na(object@cores) ||
        object@cores < 0L) {
      return("`cores` must be a non-negative integer scalar")
    }

    if (length(object@optimizeForSpeed) != 1L ||
        is.na(object@optimizeForSpeed)) {
      return("`optimizeForSpeed` must be a non-missing logical scalar")
    }

    if (length(object@deviceStr) != 1L || is.na(object@deviceStr)) {
      return("`deviceStr` must be a character scalar")
    }

    if (length(object@chunkSize) != 1L || is.na(object@chunkSize) ||
        object@chunkSize < 1L) {
      return("`chunkSize` must be a positive integer scalar")
    }

    return(TRUE)
  }
)

#' @title Build execution options
#'
#' @param cores Requested number of CPU cores
#' @param optimizeForSpeed Whether to try `torch` acceleration
#' @param deviceStr Requested `torch` device string
#' @param chunkSize Integer scalar. Requested solver batch size.
#'
#' @returns An object of class `ExecutionOptions`
#'
#' @export
#'
#' @examples
#'   exec <- ExecutionOptions(
#'     cores = 4L,
#'     optimizeForSpeed = TRUE,
#'     deviceStr = "cuda",
#'     chunkSize = 1024L
#'   )
#'
#' @rdname ExecutionOptions
ExecutionOptions <- function(cores = 1L,
                             optimizeForSpeed = TRUE,
                             deviceStr = "cuda",
                             chunkSize = 1024L) {
  methods::new(
    "ExecutionOptions",
    cores = as.integer(cores),
    optimizeForSpeed = as.logical(optimizeForSpeed),
    deviceStr = as.character(deviceStr),
    chunkSize = as.integer(chunkSize)
  )
}


# internal helper: build the object from legacy loose parameters
legacyExecutionOptions <- function(cores = 1L,
                                   optimizeForSpeed = TRUE,
                                   deviceStr = "cuda",
                                   chunkSize = 1024L) {
  ExecutionOptions(
    cores = cores,
    optimizeForSpeed = optimizeForSpeed,
    deviceStr = deviceStr,
    chunkSize = chunkSize
  )
}


# internal helper: apply runtime feasibility checks
resolveExecutionOptions <- function(executionOptions) {
  assertthat::assert_that(
    methods::is(executionOptions, "ExecutionOptions"),
    msg = "`executionOptions` must be an `ExecutionOptions` object"
  )

  effCores <- handleMultiCore(executionOptions@cores)

  torchInfo <- canUseTorch(
    optimizeForSpeed = executionOptions@optimizeForSpeed,
    deviceStr = executionOptions@deviceStr
  )

  return(list(
    "cores" = effCores,
    "useTorch" = torchInfo[["useTorch"]],
    "deviceStr" = torchInfo[["deviceStr"]],
    "chunkSize" = executionOptions@chunkSize
  ))
}


# ----------------- reduction options --------------------

#' @title Reduction options for COTAN dimensionality reduction
#'
#' @description A small parameter object bundling controls used to build the
#' reduced data matrix used by clusterizations and UMAP plots.
#'
#' This object intentionally stores only dimensionality-reduction policy
#' parameters. UMAP layout options and clustering-specific options should stay
#' outside this class.
#'
#' @slot useCoexEigen Logical scalar. Whether to use the first COEX eigenvectors
#'   instead of PCA on a selected gene matrix.
#' @slot dataMethod Character scalar. Data matrix method to use. Empty string is
#'   allowed during the compatibility phase and is resolved by the owning public
#'   function.
#' @slot numComp Integer scalar. Number of reduced components to calculate.
#' @slot genesSel Character vector. Gene-selection method or explicit gene list.
#'   Empty string is allowed during the compatibility phase and is resolved by
#'   the owning public function.
#' @slot numGenes Integer scalar. Number of genes to select when `genesSel`
#'   names a selection method.
#'
#' @name ReductionOptions-class
#'
#' @exportClass ReductionOptions
#'
#' @rdname ReductionOptions
#'
setClass(
  "ReductionOptions",
  slots = c(
    useCoexEigen = "logical",
    dataMethod = "character",
    numComp = "integer",
    genesSel = "character",
    numGenes = "integer"
  ),
  prototype = list(
    useCoexEigen = TRUE,
    dataMethod = "LogLikelihood",
    numComp = 25L,
    genesSel = "HGDI",
    numGenes = 2000L
  ),
  validity = function(object) {
    if (length(object@useCoexEigen) != 1L ||
        is.na(object@useCoexEigen)) {
      return("`useCoexEigen` must be a non-missing logical scalar")
    }

    if (length(object@dataMethod) != 1L || is.na(object@dataMethod)) {
      return("`dataMethod` must be a non-missing character scalar")
    }

    if (length(object@numComp) != 1L || is.na(object@numComp) ||
        object@numComp < 1L) {
      return("`numComp` must be a positive integer scalar")
    }

    if (length(object@genesSel) < 1L || anyNA(object@genesSel)) {
      return("`genesSel` must be a non-missing character vector")
    }

    if (length(object@numGenes) != 1L || is.na(object@numGenes) ||
        object@numGenes < 1L) {
      return("`numGenes` must be a positive integer scalar")
    }

    return(TRUE)
  }
)

#' @title Build reduction options
#'
#' @param useCoexEigen Whether to use the first COEX eigenvectors instead of PCA
#'   on a selected gene matrix.
#' @param dataMethod Data matrix method to use.
#' @param numComp Number of reduced components to calculate.
#' @param genesSel Gene-selection method or explicit gene list.
#' @param numGenes Number of genes to select when `genesSel` names a selection
#'   method.
#'
#' @returns An object of class `ReductionOptions`
#'
#' @export
#'
#' @examples
#'   redOpt <- ReductionOptions(
#'     useCoexEigen = FALSE,
#'     dataMethod = "LogNormalized",
#'     numComp = 25L,
#'     genesSel = "HGDI",
#'     numGenes = 2000L
#'   )
#'
#' @rdname ReductionOptions
#'
ReductionOptions <- function(useCoexEigen = TRUE,
                             dataMethod = "LogLikelihood",
                             numComp = 25L,
                             genesSel = "HGDI",
                             numGenes = 2000L) {
  methods::new(
    "ReductionOptions",
    useCoexEigen = as.logical(useCoexEigen),
    dataMethod = as.character(dataMethod),
    numComp = as.integer(numComp),
    genesSel = as.character(genesSel),
    numGenes = as.integer(numGenes)
  )
}


# internal helper: build the object from legacy loose parameters
legacyReductionOptions <- function(useCoexEigen,
                                   dataMethod,
                                   numComp,
                                   genesSel,
                                   numGenes) {
  ReductionOptions(
    useCoexEigen = useCoexEigen,
    dataMethod = dataMethod,
    numComp = numComp,
    genesSel = genesSel,
    numGenes = numGenes
  )
}


