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

