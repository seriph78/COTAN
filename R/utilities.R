#----------------- log functions --------------------

#' @title Logging in the `COTAN` package
#'
#' @description Logging is currently supported for all `COTAN` functions. It is
#'   possible to see the output on the terminal and/or on a log file. The level
#'   of output on terminal is controlled by the  `COTAN.LogLevel` option while
#'   the logging on file is always at its maximum verbosity
#'
#' @name LoggingFunctions
NULL

#' @details `setLoggingLevel()` sets the `COTAN` logging level. It set the
#'   `COTAN.LogLevel` options to one of the following values:
#'    * 0 - Always on log messages
#'    * 1 - Major log messages
#'    * 2 - Minor log messages
#'    * 3 - All log messages
#'
#' @param newLevel the new default logging level. It defaults to 1
#'
#' @returns `setLoggingLevel()` returns the old logging level or default level
#'   if not set yet.
#' @export
#'
#' @examples
#' setLoggingLevel(3) # for debugging purposes only
#'
#' @rdname LoggingFunctions
#'
setLoggingLevel <- function(newLevel = 1L) {
  message("Setting new log level to ", newLevel)
  oldLevel <- options(COTAN.LogLevel = newLevel)
  if (is.null(oldLevel)) {
    oldLevel <- 1L
  }
  return(invisible(oldLevel))
}


#' @details `setLoggingFile()` sets the log file for all `COTAN` output logs. By
#'   default no logging happens on a file (only on the console). Using this
#'   function `COTAN` will use the indicated file to dump the logs produced by
#'   all [logThis()] commands, independently from the log level. It stores the
#'   `connection` created by the call to [bzfile()] in the option:
#'   `COTAN.LogFile`
#'
#' @param logFileName the log file.
#'
#' @export
#'
#' @examples
#' logFile <- file.path(".", "COTAN_Test1.log")
#' setLoggingFile(logFile)
#' logThis("Some log message")
#' setLoggingFile("") # closes the log file
#' file.remove(logFile)
#'
#' @rdname LoggingFunctions
#'
setLoggingFile <- function(logFileName) {
  currentFile <- getOption("COTAN.LogFile")
  if (!is.null(currentFile)) {
    message("Closing previous log file - ", appendLF = FALSE)
    tryCatch({
      flush(currentFile)
      close(currentFile)
    }, error = function(e) {
      message("Connection to previous log file broken, will be discarded")
      options(COTAN.LogFile = NULL)
    })
  }

  message("Setting log file to be: ", logFileName)
  if (isEmptyName(logFileName)) {
    options(COTAN.LogFile = NULL)
  } else {
    options(COTAN.LogFile = file(logFileName, open = "at"))
  }
}


#' @details `logThis()` prints the given message string if the current log level
#'   is greater or equal to the given log level (it always prints its message on
#'   file if active). It uses [message()] to actually print the messages on the
#'   [stderr()] connection, so it is subject to [suppressMessages()]
#'
#' @param msg the message to print
#' @param logLevel the logging level of the current message. It defaults to 2
#' @param appendLF whether to add a new-line character at the end of the message
#'
#' @returns `logThis()` returns TRUE if the message has been printed on the
#'   terminal
#'
#' @export
#'
#' @examples
#' logThis("LogLevel 0 messages will always show, ",
#'         logLevel = 0, appendLF = FALSE)
#' suppressMessages(logThis("unless all messages are suppressed",
#'                          logLevel = 0))
#'
#' @rdname LoggingFunctions
#'
logThis <- function(msg, logLevel = 2L, appendLF = TRUE) {
  # write to log file if any
  currentFile <- getOption("COTAN.LogFile")
  if (!is.null(currentFile)) {
    tryCatch({
      tsMsg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg)
      if (isTRUE(appendLF)) { sep <- "\n" } else { sep <- "" }
      writeLines(tsMsg, currentFile, sep = sep)
      flush(currentFile)
    }, error = function(e) {
      setLoggingFile("")
    })
  }
  # set the logging level global variable
  if (is.null(getOption("COTAN.LogLevel"))) {
    setLoggingLevel() # to default
  }
  currentLevel <- getOption("COTAN.LogLevel")
  showMessage <- currentLevel >= logLevel
  if (showMessage) {
    message(msg, appendLF = appendLF)
  }
  return(invisible(showMessage))
}

#----------------- multi-threading --------------------

#' @title Handling Multi-Core Enviroments
#'
#' @description Check whether session supports multi-core evaluation and
#'   provides an effective upper bound to the number of cores.
#'
#' @details It uses [supportsMulticore()] and [availableCores()] from the
#'   package `parallelly` to actually check whether the session supports
#'   multi-core evaluation.
#'
#' @seealso the help page of [supportsMulticore()] about the flags influencing
#'   the multi-core support; e.g. the usage of `R` option
#'   `parallelly.fork.enable`.
#'
#' @param cores the number of cores asked for
#'
#' @returns the maximum sensible number of cores to use
#'
#' @importFrom parallelly availableCores
#' @importFrom parallelly supportsMulticore
#'
#' @noRd
#'
handleMultiCore <- function(cores) {
  cores <- max(1L, cores)

  if (!supportsMulticore() && cores != 1L) {
    if (is.null(getOption("COTAN.MultiCoreWarning"))) {
      warning("On this system multi-core is not currently supported;",
              " this can happen on some systems like 'windows'")
      warning("In case you might try 'options(parallelly.fork.enable = TRUE)'",
              " to enable multi-core support")
    }
    options(COTAN.MultiCoreWarning = "Published")
    warning("The number of cores used will be set 1!")
    cores <- 1L
  }

  cores <- min(cores, availableCores(omit = 1L))

  logThis(paste("Effective number of cores used:", cores), logLevel = 3L)

  return(cores)
}


#' @title Internal function to handle the [torch] library
#'
#' @description Returns whether the torch library is ready to be used:
#'   it obeys the opt-out flag set via the `COTAN.UseTorch` option
#'
#' @param optimizeForSpeed A Boolean to indicate whether to try to use the
#'   faster torch library
#' @param deviceStr The name of the device to be used by torch
#'
#' @returns A list with 2 elements:
#' * `"useTorch"`: a Boolean indicating whether the torch library can be used
#' * `"deviceStr"`: the updated name of the device to be used: if no `cuda` GPU
#'   is available it will fallback to CPU calculations
#'
#' @noRd
#'
canUseTorch <- function(optimizeForSpeed, deviceStr) {
  warnedAboutTorch <- !is.null(getOption("COTAN.TorchWarning"))

  useTorch <- isTRUE(optimizeForSpeed) &&
    requireNamespace("torch", quietly = TRUE)

  if (useTorch) {
    # if torch is not explicitly opted-in, we avoid using it as
    # there is no clean way to check if it is usable
    useTorchOpt <- getOption("COTAN.UseTorch")
    if (is.null(useTorchOpt)) {
      # default case: explicit opt-out only!
      useTorchOpt <- TRUE
    }
    if (is.character(useTorchOpt)) {
      useTorchOpt <- toupper(useTorchOpt)
      useTorchOpt <- !(useTorchOpt %in% c("FALSE", "F"))
    }
    useTorch <- isTRUE(useTorchOpt)

    if (!useTorch && !warnedAboutTorch) {
      warning("The `torch` library is installed,",
              " but has not been opted in yet")
      warning("In case you might try 'options(COTAN.UseTorch = TRUE)'",
              " to enable it")
      warnedAboutTorch <- TRUE
    }
  }

  if (useTorch) {
    tryCatch({
      if (!torch::torch_is_installed()) {
        stop("The `torch` library is installed but the required",
             " additional libraries are not avalable yet")
      }
      library("torch", character.only = TRUE)
      # Call a simple torch function to check if it's working
      if (is.null(torch::torch_tensor(1))) {
        stop("The `torch` library is installed but not working correctly")
      }
    },
    error = function(err) {
      logThis(paste("While trying to load the `torch` library", err),
              logLevel = 0L)
      if (!warnedAboutTorch) {
        warning("The `torch` library is installed,",
                " but might require further initialization")
        warning("Please look at the `torch` package installation guide",
                " to complete the installation")
      }
      useTorch <<- FALSE
    })
    warnedAboutTorch <- !useTorch
  }

  if (useTorch) {
    # Device configuration - fall-back to cpu if no cuda device is available
    if (substr(deviceStr, 1L, 4L) == "cuda" &&
        !torch::cuda_is_available()) {
      if (!warnedAboutTorch) {
        warning("The `torch` library could not find any `CUDA` device")
        warning("Falling back to CPU calculations")
        warnedAboutTorch <- TRUE
      }
      deviceStr <- "cpu"
    }
  } else {
    if(optimizeForSpeed) {
      if (!warnedAboutTorch) {
        warning("The `torch` library is not installed.")
        warnedAboutTorch <- TRUE
      }
      warning("Falling back to legacy [non-torch] code.")
    }
    deviceStr <- ""
  }

  if (warnedAboutTorch) {
    options(COTAN.TorchWarning = "Published")
  }

  return(list("useTorch" = useTorch, "deviceStr" = deviceStr))
}


#----------------- miscellanea --------------------

#' @title Internal function to handle names subset...
#'
#' @description Returns the given subset or the full list of names if none were
#'   specified
#'
#' @param names The full list of the genes
#' @param subset The names' subset. When empty all names are returned instead!
#'
#' @returns The updated list of names' subset, reordered according to the given
#'   names' list
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @noRd
#'
handleNamesSubsets <- function(names, subset = vector(mode = "character")) {
  if (is_empty(subset)) {
    subset <- names
  } else {
    assert_that(all(subset %in% names),
                msg = "Passed names are not a subset of the full list")

    subset <- names[names %in% subset]
  }
  return(subset)
}


#' @title Internal function to check if a name is empty...
#'
#' @description Returns whether the passed name is not null and has non-zero
#'   characters
#'
#' @param name the name to check
#'
#' @returns whether the passed name is not equivalent to an empty string
#'
#' @importFrom assertthat assert_that
#'
#' @noRd
#'
isEmptyName <- function(name) {
  return(!(length(name) && any(nzchar(name))))
}


#' @title Internal function to convert a named factor to a character vector
#'
#' @description Returns the factor cast to a vector along its names
#'
#' @param a `factor` object
#'
#' @returns a character vector that is preserving the names of the input
#'
#' @importFrom rlang set_names
#'
#' @importFrom assertthat assert_that
#'
#' @noRd
#'
factorToVector <- function(v) {
  assert_that(inherits(v, "factor"), msg = "Passed object is not a factor")
  return(set_names(levels(v)[v], names(v)))
}


#' @title Internal function to have nicer factor labels that have all the same
#'   number of characters
#'
#' @description Returns the factor with the new labels
#'
#' @param an `array` or `factor` object
#'
#' @returns a character vector that is preserving the names of the input
#'
#' @importFrom rlang set_names
#'
#' @importFrom assertthat assert_that
#'
#' @noRd
#'
niceFactorLevels <- function(v) {
  names <- names(v)
  if (inherits(v, "factor")) {
    v <- factorToVector(v)
  }
  nv <- suppressWarnings(as.numeric(v))
  if (!anyNA(nv) && all(as.integer(nv) == nv)) {
    numDigits <- floor(log10(max(nv))) + 1L
    v <- formatC(nv, width = numDigits, flag = "0")
  } else if (is.character(v)) {
    numChars <- max(nchar(v))
    v <- str_pad(v, width = numChars, side = "left", pad = "_")
  }
  names(v) <- names
  return(factor(v))
}


#' @details `getColumnFromDF()` is a function to extract a column from a
#'   `data.frame`, while keeping the `rowNames` as `vector` names
#'
#' @param df the `data.frame`
#' @param colName the name of the new or existing column in the `data.frame`
#'
#' @returns `getColumnFromDF()` returns the column in the `data.frame` as named
#'   `array`, `NULL` if the wanted column is not available
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @rdname HandleMetaData
#'
getColumnFromDF <- function(df, colName) {
  if (is_empty(df)) {
    return(NULL)
  } else {
    retArray <- df[[colName]]
    if (!is_empty(retArray)) {
      names(retArray) <- rownames(df)
    }
    return(retArray)
  }
}


#' @details `setColumnInDF()` is a function to append, if missing, or resets, if
#'   present, a column into a `data.frame`, whether the `data.frame` is empty or
#'   not. The given `rowNames` are used only in the case the `data.frame` has
#'   only the default row numbers, so this function cannot be used to override
#'   row names
#'
#' @param df the `data.frame`
#' @param colToSet the column to add
#' @param colName the name of the new or existing column in the `data.frame`
#' @param rowNames when not empty, if the input `data.frame` has no real row
#'   names, the new row names of the resulting `data.frame`
#'
#' @returns `setColumnInDF()` returns the updated, or the newly created,
#'   `data.frame`
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @rdname HandleMetaData
#'
setColumnInDF <- function(df, colToSet, colName,
                          rowNames = vector(mode = "character")) {
  if (is_empty(df)) {
    df <- set_names(data.frame(colToSet), colName)
  } else {
    if (colName %in% colnames(df)) {
      df[[colName]] <- colToSet
    } else {
      df <- cbind(df, set_names(data.frame(colToSet), colName))
    }
  }

  # default assigned rownames are numbers...
  if (!is_empty(rowNames) && is.integer(attr(df, "row.names"))) {
    rownames(df) <- rowNames
  }

  return(df)
}


#------------------- clusters utilities ----------

#' @title *Clusters* utilities
#'
#' @description Handle *clusterization* <-> *clusters* `list` conversions,
#'   *clusters* grouping and merge
#'
#' @name ClustersList
#'
#' @examples
#' ## create a clusterization
#' clusters <- paste0("",sample(7, 100, replace = TRUE))
#' names(clusters) <- paste0("E_",formatC(1:100,  width = 3, flag = "0"))
#'
#' ## create a clusters list from a clusterization
#' clustersList <- toClustersList(clusters)
#' head(clustersList, 1)
#'
#' ## recreate the clusterization from the cluster list
#' clusters2 <- fromClustersList(clustersList, names(clusters))
#' all.equal(factor(clusters), clusters2)
#'
#' cl1Size <- length(clustersList[["1"]])
#'
#' ## establish the permutation that groups clusters together
#' perm <- groupByClusters(clusters)
#' !is.unsorted(head(names(clusters)[perm],cl1Size))
#' head(clusters[perm], cl1Size)
#'
#' ## it is possible to have the list of the element names different
#' ## from the names in the clusters list
#' selectedNames <- paste0("E_",formatC(11:110,  width = 3, flag = "0"))
#' perm2 <- groupByClustersList(selectedNames, toClustersList(clusters))
#' all.equal(perm2[91:100], c(91:100))
#'
#' ## is is possible to merge a few clusters together
#' clustersMerged <- mergeClusters(clusters, names = c("7", "2"),
#'                                 mergedName = "7__2")
#' sum(table(clusters)[c(2, 7)]) == table(clustersMerged)[["7__2"]]
#'
#' ## it is also possible to do multiple merges at once!
#' ## Note the default new clusters' names
#' clustersMerged2 <-
#'   multiMergeClusters(clusters2, namesList = list(c("2", "7"),
#'                                                  c("1", "3", "5")))
#' table(clustersMerged2)
#'
#'
NULL

#' @details `toClustersList()` given a *clusterization*, creates a `list` of
#'   *clusters* (i.e. for each *cluster*, which elements compose the *cluster*)
#'
#' @param clusters A named `vector` or `factor` that defines the *clusters*
#'
#' @returns `toClustersList()` returns a `list` of clusters
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom assertthat assert_that
#'
#' @rdname ClustersList
#'
toClustersList <- function(clusters) {
  assert_that(!is_empty(names(clusters)),
              msg = "passed clusterization has no names")

  clustersNames <- levels(factor(clusters))

  getCl <- function(cl, clusters) {
    names(clusters)[clusters %in% cl]
  }

  clustersList <- set_names(lapply(clustersNames, getCl, clusters),
                            clustersNames)

  return(clustersList)
}

#' @details `fromClustersList()` given a `list` of *clusters* returns a
#'   *clusterization* (i.e. a named `vector` that for each element indicates to
#'   which cluster it belongs)
#'
#' @param clustersList A named `list` whose elements define the various clusters
#' @param elemNames A `list` of names to which associate a cluster
#' @param throwOnOverlappingClusters When `TRUE`, in case of overlapping
#'   clusters, the function `fromClustersList` and `groupByClustersList` will
#'   throw. This is the default. When FALSE, instead, in case of overlapping
#'   clusters, `fromClustersList` will return the last cluster to which each
#'   element belongs, while `groupByClustersList` will return a vector of
#'   positions that is longer than the given `elemNames`
#'
#' @returns `fromClustersList()` returns a clusterization. If the given
#'   `elemNames` contain values not present in the `clustersList`, those will be
#'   marked as `"-1"`
#'
#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom assertthat assert_that
#'
#' @rdname ClustersList
#'
fromClustersList <- function(clustersList,
                             elemNames = vector(mode = "character"),
                             throwOnOverlappingClusters = TRUE) {
  clustersNames <- names(clustersList)

  assert_that(!is_empty(clustersNames),
              msg = "Passed clusterization has no names")

  if (is_empty(elemNames)) {
    elemNames <- unlist(clustersList, use.names = FALSE)
  }

  clusters <- set_names(rep.int("-1", length(elemNames)), elemNames)

  for (clName in clustersNames) {
    cluster <- clustersList[[clName]]
    cluster <- cluster[cluster %in% elemNames]
    assert_that((!throwOnOverlappingClusters ||
                   all(clusters[cluster] == "-1")),
                msg = "Found overlapping clusters")
    clusters[cluster] <- clName
  }

  return(factor(clusters))
}

#' @details `groupByClusters()` given a *clusterization* returns a permutation,
#'   such that using the permutation on the input the *clusters* are grouped
#'   together
#'
#' @param elemNames A `list` of names to which associate a cluster
#' @param clustersList A named `list` whose elements define the various clusters
#' @param throwOnOverlappingClusters When `TRUE`, in case of overlapping
#'   clusters, the function `fromClustersList` and `groupByClustersList` will
#'   throw. This is the default. When FALSE, instead, in case of overlapping
#'   clusters, `fromClustersList` will return the last cluster to which each
#'   element belongs, while `groupByClustersList` will return a vector of
#'   positions that is longer than the given `elemNames`
#'
#' @returns `groupByClusters()` and `groupByClustersList()` return a permutation
#'   that groups the clusters together. For each cluster the positions are
#'   guaranteed to be in increasing order. In case, all elements not
#'   corresponding to any cluster are grouped together as the last group
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @rdname ClustersList
#'
groupByClustersList <- function(elemNames, clustersList,
                                throwOnOverlappingClusters = TRUE) {
  assert_that(!is_empty(elemNames), msg = "passed no elemNames")
  assert_that(!is_empty(clustersList), msg = "passed no clustersList")

  positions <- vector(mode = "integer")

  for (cluster in clustersList) {
    clPos <- which(elemNames %in% cluster)
    # clPos should be already sorted
    assert_that(!throwOnOverlappingClusters || !any(clPos %in% positions),
                msg = "Found overlapping clusters")
    positions <- append(positions, clPos)
  }

  # add all non-clustered elements as tail group!
  positions <- append(positions, setdiff(seq_along(elemNames), positions))

  return(positions)
}

#' @details `groupByClustersList()` given the elements' names and a `list` of
#'   *clusters* returns a permutation, such that using the permutation on the
#'   given names the *clusters* are grouped together.
#'
#' @param clusters A named `vector` or `factor` that defines the *clusters*.
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @rdname ClustersList
#'
groupByClusters <- function(clusters) {
  return(groupByClustersList(names(clusters), toClustersList(clusters)))
}


#' @details `mergeClusters()` given a *clusterization*, creates a new one where
#'   the given *clusters* are merged.
#'
#' @param clusters A named `vector` or `factor` that defines the *clusters*
#' @param names A list of *clusters* names to be merged
#' @param mergedName The name of the new merged clusters
#'
#' @returns `mergeClusters()` returns a new *clusterization* with the wanted
#'   *clusters* being merged. If less than 2 *cluster* names were passed the
#'   function will emit a warning and return the initial *clusterization*
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @rdname ClustersList
#'
mergeClusters <- function(clusters, names, mergedName = "") {
  if (!inherits(clusters, "factor")) {
    clusters <- factor(clusters)
  }

  effNames <- names[names %in% levels(clusters)]
  if (is_empty(effNames) || length(effNames) < 2L) {
    warning("Passed a list of clusters to merge with less than 2 elements",
            " actually present in the clusterization")
    # nothing to do...
    return(factor(clusters))
  }

  if (isEmptyName(mergedName)) {
    effNames <- sort(effNames)
    mergedName <- paste0(paste(effNames, collapse = "_"), "-merge")
  }

  # cannot assign new levels into a factor
  cclusters <- factorToVector(clusters)
  cclusters[clusters %in% effNames] <- mergedName

  return(factor(cclusters))
}


#' @details `multiMergeClusters()` given a *clusterization*, creates a new one
#'   where the given sets of *clusters* are merged.
#'
#' @param clusters A named `vector` or `factor` that defines the *clusters*
#' @param namesList A `list` of `list`s of *clusters* names to be respectively
#'   merged
#' @param mergedNames The names of the new merged *clusters*
#'
#' @returns `multiMergeClusters()` returns a new *clusterization* with the
#'   wanted *clusters* being merged by consecutive iterations of
#'   [mergeClusters()] on the given `namesList`
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @rdname ClustersList
#'
multiMergeClusters <- function(clusters, namesList, mergedNames = NULL) {
  if (is_empty(namesList)) {
    warning("Passed no clusters lists to merge")
  }

  if (is_empty(mergedNames)) {
    mergedNames <- rep_len("", length(namesList))
  }

  assert_that(length(mergedNames) == length(namesList),
              msg = paste("When given 'mergedNames' must have",
                          "the same length as 'namesList'"))

  for (i in seq_along(namesList)) {
    clusters <- mergeClusters(clusters,
                              names = namesList[[i]],
                              mergedName = mergedNames[[i]])
  }

  return(clusters)
}

#----------------- plot utilities --------------------

#' @details `plotTheme()` returns the appropriate theme for the selected plot
#'   kind. Supported kinds are:  `"common"`, `"pca"`, `"genes"`, `"UDE"`,
#'   `"heatmap"`, `"GDI"`, `"UMAP"`, `"size-plot"`
#'
#' @seealso [ggplot2::theme()] and [ggplot2::ggplot()]
#'
#' @param plotKind a string indicating the plot kind
#' @param textSize axes and strip text size (default=14)
#'
#' @returns `plotTheme()` returns a `ggplot2::theme` object
#
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 unit
#'
#' @importFrom ggthemes theme_tufte
#'
#' @export
#'
#' @examples
#' theme <- plotTheme("pca")
#'
#' @rdname HeatmapPlots
#'
plotTheme <- function(plotKind = "common", textSize = 14L) {
  myDarkBlue <- "#3C5488FF"
  ts <- textSize

  basicTheme <- theme(
    axis.text.x  = element_text(size = ts, angle = 0L, hjust = 0.5, vjust = 0.5,
                                face = "plain", colour = myDarkBlue),
    axis.text.y  = element_text(size = ts, angle = 0L, hjust = 0.0, vjust = 0.5,
                                face = "plain", colour = myDarkBlue),
    axis.title.x = element_text(size = ts, angle = 0L, hjust = 0.5, vjust = 0.0,
                                face = "plain", colour = myDarkBlue),
    axis.title.y = element_text(size = ts, angle =90L, hjust = 0.5, vjust = 0.5,
                                face = "plain", colour = myDarkBlue))

  if (plotKind == "common") {
    return(basicTheme)
  }

  if (plotKind == "pca") {
    return(basicTheme +
           theme(legend.title = element_blank(),
                 legend.text = element_text(size = 12L, face = "italic",
                                            color = myDarkBlue),
                 legend.position = "bottom"))
  }

  if (plotKind == "genes") {
    return(basicTheme +
           theme(plot.title = element_text(size = 20L, hjust = 0.02,
                                           vjust = -10.0, face = "italic",
                                           color = myDarkBlue),
                 plot.subtitle = element_text(vjust = -15.0, hjust = 0.01,
                                              color = "darkred")))
  }

  if (plotKind == "UDE") {
    return(basicTheme +
           theme(plot.title   = element_text(size = 20L, color = myDarkBlue),
                 legend.title = element_text(size = 14L, color = myDarkBlue,
                                             face = "italic"),
                 legend.text  = element_text(size = 11L, color = myDarkBlue),
                 legend.key.width = unit(2.0, "mm"),
                 legend.position = "right"))
  }

  if (plotKind == "heatmap") {
    return(basicTheme +
           theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.spacing = unit(0.0, "lines"),
                 strip.background = element_rect(fill = "#8491B44C"),
                 strip.text.y = element_text(size = ts, colour = myDarkBlue),
                 strip.text.x = element_text(size = ts, angle = 90L,
                                             colour = myDarkBlue),
                 legend.text = element_text(color = myDarkBlue,
                                            face = "italic"),
                 legend.position = "bottom",
                 legend.title = element_blank(),
                 legend.key.height = unit(2.0, "mm")))
  }

  if (plotKind == "GDI") {
    return(basicTheme +
           theme(legend.title = element_blank(),
                 plot.title = element_text(size = ts + 2L,
                                           face = "bold.italic",
                                           color = myDarkBlue),
                 legend.text = element_text(color = myDarkBlue,
                                            face = "italic"),
                 legend.position = "bottom"))
  }

  if (plotKind == "UMAP") {
    return(basicTheme +
             theme(legend.title = element_blank(),
                   plot.title = element_text(size = ts + 2L,
                                             face = "bold.italic",
                                             color = myDarkBlue),
                   legend.text = element_text(color = myDarkBlue,
                                              face = "italic"),
                   legend.position = "right"))
  }

  if (plotKind == "size-plot") {
    return(ggthemes::theme_tufte() +
           theme(legend.position = "none"))
                 # axis.text.x  = element_blank(),
                 # axis.ticks.x = element_blank()) )
  }

  warning("plotTheme: no match found in listed themes for: ", plotKind)
  return(basicTheme)
}

#' @title getColorsVector
#'
#' @description This function returns a list of colors based on the
#'   [brewer.pal()] function
#'
#' @details The colors are taken from the [brewer.pal.info()] sets with
#'   `Set1`, `Set2`, `Set3` placed first.
#'
#' @param numNeededColors The number of returned colors. If omitted it returns
#'   all available colors
#'
#' @returns an array of `RGB` colors of the wanted size
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' colorsVector <- getColorsVector(17)
#'
#' @rdname getColorsVector
#'
getColorsVector <- function(numNeededColors = 0L) {
  qualColPalets <- brewer.pal.info[brewer.pal.info[["category"]] == "qual", ]
  numColPalets <- nrow(qualColPalets)

  qualColPalets <- rbind(qualColPalets[(numColPalets - 2L):numColPalets, ],
                         qualColPalets[1L:(numColPalets - 3L), ])

  colVector <- unlist(Map(brewer.pal, qualColPalets[["maxcolors"]],
                          rownames(qualColPalets)))

  rm(qualColPalets, numColPalets)

  assert_that(numNeededColors <= length(colVector),
              msg = paste("Needed more colors than the number",
                          "of possible supported colors:", length(colVector)))

  if (numNeededColors == 0L) {
    return(colVector)
  }
  return(head(colVector, numNeededColors))
}
