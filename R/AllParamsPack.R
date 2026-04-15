# ----------------- execution options --------------------

#' @title Execution options for heavy COTAN computations
#'
#' @description A small parameter object bundling execution-related controls
#' such as multi-core usage and the request to use `torch` on a given device.
#'
#' @slot cores Integer scalar. Requested number of CPU cores.
#' @slot optimizeForSpeed Logical scalar. Whether to try `torch` acceleration.
#' @slot deviceStr Character scalar. Requested `torch` device string.
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
    deviceStr = "character"
  ),
  prototype = list(
    cores = 1L,
    optimizeForSpeed = TRUE,
    deviceStr = "cuda"
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

    TRUE
  }
)

#' @title Build execution options
#'
#' @param cores Requested number of CPU cores
#' @param optimizeForSpeed Whether to try `torch` acceleration
#' @param deviceStr Requested `torch` device string
#'
#' @returns An object of class `ExecutionOptions`
#'
#' @export
#'
#' @rdname ExecutionOptions
ExecutionOptions <- function(cores = 1L,
                             optimizeForSpeed = TRUE,
                             deviceStr = "cuda") {
  methods::new(
    "ExecutionOptions",
    cores = as.integer(cores),
    optimizeForSpeed = as.logical(optimizeForSpeed),
    deviceStr = as.character(deviceStr)
  )
}


# internal helper: build the object from legacy loose parameters
legacyExecutionOptions <- function(cores = 1L,
                                   optimizeForSpeed = TRUE,
                                   deviceStr = "cuda") {
  ExecutionOptions(
    cores = cores,
    optimizeForSpeed = optimizeForSpeed,
    deviceStr = deviceStr
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
    "deviceStr" = torchInfo[["deviceStr"]]
  ))
}

