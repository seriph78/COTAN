#----------------- log functions --------------------

#' Logging in the `COTAN` package
#'
#' @description Logging is currently supported for all `COTAN` functions. It is
#'   possible to see the output on the terminal and/or on a log file. The level
#'   of output on terminal is controlled by the  `COTAN.LogLevel` option while
#'   the logging on file is always at its maximum verbosity
#'
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
#' setLoggingFile("./COTAN_Test1.log") # for debugging purposes only
#' logThis("Some log message")
#' setLoggingFile("") # closes the log file
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
  emptyName <- !(length(logFileName) && any(nzchar(logFileName)))
  if (emptyName) {
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
      writeLines(msg, currentFile)
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

#----------------- miscellanea --------------------

#' handleMultiCore
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
              " this can happen on some systems like 'windows'.\n",
              " In case you can try 'options(parallelly.fork.enable = TRUE)'",
              " to enable multi-core support.\n",
              " The number of cores used will be set 1!")
    }
    options(COTAN.MultiCoreWarning = "Published")
    cores <- 1L
  }

  cores <- min(cores, availableCores(omit = 1L))

  logThis(paste("Effective number of cores used:", cores), logLevel = 3L)

  return(cores)
}


#' Internal function to handle names subset...
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
#' @importFrom assertthat assert_that
#'
#' @noRd
#'
handleNamesSubsets <- function(names, subset = c()) {
  if (is_empty(subset)) {
    subset <- names
  } else {
    assert_that(all(subset %in% names),
                msg = "Passed genes are not a subset of the full list")

    subset <- names[names %in% subset]
  }
  return(subset)
}


#' @details `setColumnInDF()` is a function to append, if missing, or resets, if
#'   present, a column into a `data.frame`, whether the `data.frame` is empty or
#'   not. The given `rowNames` are used only in the case the `data.frame` has
#'   only the default row numbers, so this function cannot be used to override
#'   row names
#'
#' @param df the `data.frame`
#' @param colToSet the the column to add
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
setColumnInDF <- function(df, colToSet, colName, rowNames = c()) {
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

#' Handle clusterization <-> clusters list conversions and clusters grouping
#'
#' @description `toClustersList` given a clusterization, creates a `list` of
#'   clusters (i.e. for each cluster, which elements compose the cluster).
#'
#'   `fromClustersList` given a `list` of clusters returns a clusterization (i.e
#'   a named `vector` that for each element indicates to which cluster it
#'   belongs).
#'
#'   `groupByClusters` given a clusterization returns a permutation, such that
#'   using the permutation on the input the clusters are grouped together.
#'
#'   `groupByClustersList` given the elements' names and a `list` of clusters
#'   returns a permutation, such that using the permutation on the given names
#'   the clusters are grouped together.
#'
#' @usage toClustersList(clusters)
#'
#' @param clusters A named `vector` or `factor` that defines the clusters.
#' @param clustersList A named `list` whose elements define the various
#'   clusters.
#' @param elemNames A `list` of names to which associate a cluster.
#' @param throwOnOverlappingClusters When `TRUE`, in case of overlapping
#'   clusters, the function `fromClustersList` and `groupByClustersList` will
#'   throw. This is the default. When FALSE, instead, in case of overlapping
#'   clusters, `fromClustersList` will return the last cluster to which each
#'   element belongs, while `groupByClustersList` will return a vector of
#'   positions that is longer than the given `elemNames`.
#'
#' @returns `toClustersList` returns a `list` of clusters.
#'
#'   `fromClustersList` returns a clusterization. If the given `elemNames`
#'   contain values not present in the `clustersList`, those will be marked as
#'   `"not_clustered"`
#'
#'   `groupByClusters` and `groupByClustersList` return a permutation that
#'   groups the clusters together. For each cluster the positions are guaranteed
#'   to be in increasing order. In case, all elements not corresponding to any
#'   cluster are grouped together as the last group.
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
#' all.equal(clusters, clusters2)
#'
#' cl1Size <- length(clustersList[["1"]])
#'
#' ## establish the permutation that groups clusters together
#' perm <- groupByClusters(clusters)
#' is.unsorted(head(names(clusters)[perm],cl1Size))
#' head(clusters[perm], cl1Size)
#'
#' ## it is possible to have the list of the element names different
#' ## from the names in the clusters list
#' selectedNames <- paste0("E_",formatC(11:110,  width = 3, flag = "0"))
#' perm2 <- groupByClustersList(selectedNames, toClustersList(clusters))
#' all.equal(perm2[91:100], c(91:100))
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
    names(clusters[clusters %in% cl])
  }

  clustersList <- set_names(lapply(clustersNames, getCl, clusters),
                            clustersNames)

  return(clustersList)
}


#' @export
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom assertthat assert_that
#'
#' @rdname ClustersList
#'
fromClustersList <- function(clustersList, elemNames = c(),
                             throwOnOverlappingClusters = TRUE) {
  clustersNames <- names(clustersList)

  assert_that(!is_empty(clustersNames),
              msg = "Passed clusterization has no names")

  if (is_empty(elemNames)) {
    elemNames <- unlist(clustersList, use.names = FALSE)
  }

  clusters <- set_names(rep.int("not_clustered", length(elemNames)), elemNames)

  for (clName in clustersNames) {
    cluster <- clustersList[[clName]]
    cluster <- cluster[cluster %in% elemNames]
    assert_that((!throwOnOverlappingClusters ||
                   all(clusters[cluster] == "not_clustered")),
                msg = "Found overlapping clusters")
    clusters[cluster] <- clName
  }

  return(clusters)
}


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

  positions <- c()

  for (cluster in clustersList) {
    clPos <- which(elemNames %in% cluster)
    # clPos should be already sorted
    assert_that(!throwOnOverlappingClusters || !any(clPos %in% positions),
                msg = "Found overlapping clusters")
    positions <- append(positions, clPos)
  }

  # add all non-clustered elements as tail group!
  positions <- append(positions, setdiff(seq_len(length(elemNames)), positions))

  return(positions)
}

#' @export
#'
#' @importFrom rlang is_empty
#'
#' @rdname ClustersList
#'
groupByClusters <- function(clusters) {
  return(groupByClustersList(names(clusters), toClustersList(clusters)))
}


#------------------- numeric utilities ----------

#' funProbZero
#'
#' @description Private function that gives the probability of a sample gene
#'   count being zero given the given the dispersion and mu
#'
#' @details Using \eqn{d}{`disp`} for `disp` and \eqn{\mu}{`mu`} for `mu`,
#'   it returns:
#'   \eqn{(1 + d \mu)^{-\frac{1}{d}}}{(1 + `disp` * `mu`)^(-1 / `disp`)}
#'   when \eqn{d > 0}{`disp > 0`} and
#'   \eqn{\exp{((d - 1) \mu)}}{exp((`disp` - 1) * `mu`)} otherwise.
#'   The function is continuous in \eqn{d = 0}{`disp = 0`},
#'   increasing in \eqn{d}{`disp`} and decreasing in \eqn{\mu}{`mu`}.
#'   It returns 0 when \eqn{d = -\infty}{`disp = -Inf`} or
#'   \eqn{\mu = \infty}{`mu = Inf`}.
#'   It returns 1 when \eqn{\mu = 0}{`mu = 0`}.
#'
#' @param disp the estimated dispersion (can be a 𝑛-sized vector)
#' @param mu the lambda times nu value  (can be a 𝑛×𝑚 matrix)
#'
#' @returns the probability (matrix) that a count is identically zero
#'
#' @rdname funProbZero
#'
funProbZero <- function(disp, mu) {
  neg <- (disp <= 0.0)
  ad <- abs(disp)
  return(( neg) * (exp(-(1.0 + ad) * mu)) +
         (!neg) * (1.0 + ad * mu)^(-1.0 / ad))
}


#' dispersionBisection
#'
#' @description Private function invoked by [estimateDispersionBisection()] for
#'   the estimation of 'dispersion' slot of a `COTAN` object via a bisection
#'   solver
#'
#' @details The goal is to find a dispersion value that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param sumZeros the number of cells that didn't express the gene
#' @param lambda the lambda parameter
#' @param nu the nu vector
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion value
#'
#' @noRd
#'
dispersionBisection <-
  function(sumZeros,
           lambda,
           nu,
           threshold = 0.001,
           maxIterations = 100L) {
  if (sumZeros == 0L) {
    # cannot match exactly zero prob of zeros with finite values
    return(-Inf)
  }
  mu <- lambda * nu

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  disp1 <- 0.0
  diff1 <- sum(funProbZero(disp1, mu)) - sumZeros
  if (abs(diff1) <= threshold) {
    return(disp1)
  }

  # we assume error is an increasing function of disp
  disp2 <- -1.0 * sign(diff1)
  iter <- 1L
  repeat {
    diff2 <- sum(funProbZero(disp2, mu)) - sumZeros

    if (diff2 * diff1 < 0.0) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached ",
           "while finding the solution straddling intervals")
    }
    iter <- iter + 1L

    disp1 <- disp2 # disp2 is closer to producing 0.0

    disp2 <- 2.0 * disp2 # we double at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  iter <- 1L
  repeat {
    disp <- (disp1 + disp2) / 2.0

    diff <- sum(funProbZero(disp, mu)) - sumZeros

    if (abs(diff) <= threshold) {
      return(disp)
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached ",
           "while finding the solution straddling intervals")
    }
    iter <- iter + 1L

    # drop same sign diff point
    if (diff * diff2 > 0.0) {
      disp2 <- disp
    } else {
      disp1 <- disp
    }
  }
}


#' parallelDispersionBisection
#'
#' @description Private function invoked by [estimateDispersionBisection()] for
#'   the estimation of the `dispersion` slot of a `COTAN` object via a parallel
#'   bisection solver
#'
#' @details The goal is to find a dispersion array that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param genes names of the relevant genes
#' @param sumZeros the number of cells that didn't express the relevant gene
#' @param lambda the lambda vector
#' @param nu the nu vector
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion values
#'
#' @noRd
#'
parallelDispersionBisection <-
  function(genes,
           sumZeros,
           lambda,
           nu,
           threshold = 0.001,
           maxIterations = 100L) {
  sumZeros <- sumZeros[genes]
  lambda <- lambda[genes]

  goodPos <- sumZeros != 0L

  # cannot match exactly zero prob of zeros with finite values
  # so we ignore the rows with no zeros from the solver and return -Inf
  output <- rep(-Inf, length(sumZeros))

  if (sum(goodPos) == 0L) {
    return(output)
  }

  sumZeros <- sumZeros[goodPos]
  lambda <- lambda[goodPos]

  mu <- lambda %o% nu

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  disps1 <- rep(0.0, length(sumZeros))
  diffs1 <- rowSums(funProbZero(disps1, mu)) - sumZeros
  if (all(abs(diffs1) <= threshold)) {
    output[goodPos] <- disps1
    return(output)
  }

  # we assume error is an increasing function of the dispersion
  disps2 <- -1.0 * sign(diffs1)
  diffs2 <- diffs1
  runPos <- rep(TRUE, length(diffs1))
  iter <- 1L
  repeat {
    diffs2[runPos] <- (rowSums(funProbZero(disps2[runPos],
                                           mu[runPos, , drop = FALSE])) -
                       sumZeros[runPos])

    runPos <- (diffs2 * diffs1 >= 0.0)

    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding",
           " the solution straddling intervals")
    }
    iter <- iter + 1L

    disps1[runPos] <- disps2[runPos] # disps2 are closer to producing 0

    disps2[runPos] <- 2.0 * disps2[runPos] # we double at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  runNum <- length(diffs1)
  runPos <- rep(TRUE, runNum)
  disps <- disps1
  diffs <- diffs1
  iter <- 1L
  repeat {
    disps[runPos] <- (disps1[runPos] + disps2[runPos]) / 2.0

    diffs[runPos] <- (rowSums(funProbZero(disps[runPos],
                                          mu[runPos, , drop = FALSE])) -
                      sumZeros[runPos])

    runPos <- abs(diffs) > threshold
    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solutions")
    }
    iter <- iter + 1L

    # drop same sign diff point
    pPos <- runPos & (diffs * diffs2 > 0.0)
    disps2[pPos] <- disps[pPos]

    nPos <- runPos & !pPos
    disps1[nPos] <- disps[nPos]
  }

  output[goodPos] <- disps
  return(output)
}


#' nuBisection
#'
#' @description Private function invoked by [estimateNuBisection()] for the
#'   estimation of `nu` slot of a `COTAN` object via a bisection solver
#'
#' @details The goal is to find a nu value that reduces to zero the difference
#'   between the number of estimated and counted zeros
#'
#' @param sumZeros the number non expressed genes in the cell
#' @param lambda the lambda vector
#' @param dispersion the dispersion vector
#' @param initialGuess the initial guess for nu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the nu value
#'
#' @noRd
#'
nuBisection <-
  function(sumZeros,
           lambda,
           dispersion,
           initialGuess,
           threshold = 0.001,
           maxIterations = 100L) {
  if (sumZeros == 0L) {
    # cannot match exactly zero prob of zeros with finite values
    return(Inf)
  }

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  nu1 <- initialGuess
  diff1 <- sum(funProbZero(dispersion, nu1 * lambda)) - sumZeros
  if (abs(diff1) <= threshold) {
    return(nu1)
  }

  factor <- 2.0 ^ sign(diff1)
  nu2 <- nu1 * factor # we assume error is an decreasing function of nu
  iter <- 1L
  repeat {
    diff2 <- sum(funProbZero(dispersion, nu2 * lambda)) - sumZeros

    if (diff2 * diff1 < 0.0) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while ",
           "finding the solution straddling intervals")
    }
    iter <- iter + 1L

    nu1 <- nu2 # nu2 is closer to producing 0

    nu2 <- nu2 * factor # we double/half at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  iter <- 1L
  repeat {
    nu <- (nu1 + nu2) / 2.0

    diff <- sum(funProbZero(dispersion, nu * lambda)) - sumZeros

    if (abs(diff) <= threshold) {
      return(nu)
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while ",
           "finding the solution straddling intervals")
    }
    iter <- iter + 1L

    # drop same sign diff point
    if (diff * diff2 > 0L) {
      nu2 <- nu
    } else {
      nu1 <- nu
    }
  }
}


#' parallelNuBisection
#'
#' @description Private function invoked by [estimateNuBisection()] for the
#'   estimation of `nu` slot of a `COTAN` object via a parallel bisection solver
#'
#' @details The goal is to find a nu array that reduces to zero the difference
#'   between the number of estimated and counted zeros
#'
#' @param cells names of the relevant cells
#' @param sumZeros the number of genes not expressed in the relevant cell
#' @param lambda the lambda vector
#' @param dispersion the dispersion vector
#' @param initialGuess the initial guess vector for nu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion values
#'
#' @importFrom assertthat assert_that
#'
#' @noRd
#'
parallelNuBisection <-
  function(cells,
           sumZeros,
           lambda,
           dispersion,
           initialGuess,
           threshold = 0.001,
           maxIterations = 100L) {
  sumZeros <- sumZeros[cells]
  initialGuess <- initialGuess[cells]

  assert_that(all(initialGuess > 0.0),
              msg = "initialGuess must hold only positive values")

  goodPos <- sumZeros != 0L

  # cannot match exactly zero prob of zeros with finite values
  output <- rep(Inf, length(initialGuess))

  if (sum(goodPos) == 0L) {
    return(output)
  }

  sumZeros <- sumZeros[goodPos]

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  nus1 <- initialGuess[goodPos]

  diffs1 <- colSums(funProbZero(dispersion, lambda %o% nus1)) - sumZeros
  if (all(abs(diffs1) <= threshold)) {
    output[goodPos] <- nus1
    return(output)
  }

  # we assume error is an increasing function of the dispersion
  factors <- 2.0 ^ sign(diffs1)
  nus2 <- nus1 * factors
  diffs2 <- diffs1
  runPos <- rep(TRUE, length(diffs1))
  iter <- 1L
  repeat {
    diffs2[runPos] <- (colSums(funProbZero(dispersion,
                                           lambda %o% nus2[runPos])) -
                       sumZeros[runPos])

    runPos <- (diffs2 * diffs1 >= 0.0)

    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while ",
           "finding the solution straddling intervals")
    }
    iter <- iter + 1L

    nus1[runPos] <- nus2[runPos] # nus2 are closer to producing 0

    nus2[runPos] <- nus2[runPos] * factors[runPos] # we double/half at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  runNum <- length(diffs1)
  runPos <- rep(TRUE, runNum)
  nus <- nus1
  diffs <- diffs1
  iter <- 1L
  repeat {
    nus[runPos] <- (nus1[runPos] + nus2[runPos]) / 2.0

    diffs[runPos] <- (colSums(funProbZero(dispersion, lambda %o% nus[runPos])) -
                      sumZeros[runPos])

    runPos <- (abs(diffs) > threshold)
    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solutions")
    }
    iter <- iter + 1L

    # drop same sign diff point
    pPos <- runPos & (diffs * diffs2 > 0.0)
    nus2[pPos] <- nus[pPos]

    nPos <- runPos & !pPos
    nus1[nPos] <- nus[nPos]
  }

  output[goodPos] <- nus
  return(output)
}



#' cosineDissimilarity
#'
#' @param m a matrix
#'
#' @returns The dissimilarity matrix between columns' data
#'
#' @export
#'
#' @importFrom Matrix t
#'
#' @importFrom stats as.dist
#'
#' @examples
#' mat <- matrix(c(1:25), nrow = 5, ncol = 5,
#'               dimnames = list(paste0("row.", c(1:5)),
#'                               paste0("col.", c(1:5))))
#' dist <- cosineDissimilarity(mat)
#'
#' @rdname cosineDissimilarity
#'
cosineDissimilarity <- function(m) {
  m <- as.matrix(t(m))
  sim <- m / sqrt(rowSums(m^2.0))
  sim <- tcrossprod(sim)
  return(as.dist(1.0 - sim))
}

#----------------- plot utilities --------------------

#' plotTheme
#'
#' @description This function returns the appropriate theme for the selected
#'   plot kind
#'
#' @details Supported kinds are:
#' * `common`
#' * `pca`
#' * `genes`
#' * `UDE`
#' * `heatmap`
#' * `GDI`
#' * `size-plot`
#'
#' @seealso [ggplot2::theme()] and [ggplot2::ggplot()]
#'
#' @param plotKind a string indicating the plot kind
#' @param textSize axes and strip text size (default=14)
#'
#' @returns a `ggplot2` theme object
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
#' @rdname plotTheme
#'
plotTheme <- function(plotKind = "common", textSize = 14L) {
  myDarkBlue <- "#3C5488FF"
  ts <- textSize

  basicTheme <- theme(
    axis.text.x  = element_text(size = ts, angle = 0L,  hjust = .5, vjust = .5,
                                face = "plain", colour = myDarkBlue),
    axis.text.y  = element_text(size = ts, angle = 0L,  hjust = .0, vjust = .5,
                                face = "plain", colour = myDarkBlue),
    axis.title.x = element_text(size = ts, angle = 0L,  hjust = .5, vjust = .0,
                                face = "plain", colour = myDarkBlue),
    axis.title.y = element_text(size = ts, angle = 90L, hjust = .5, vjust = .5,
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

  if (plotKind == "size-plot") {
    return(ggthemes::theme_tufte() +
           theme(legend.position = "none"))
                 # axis.text.x  = element_blank(),
                 # axis.ticks.x = element_blank()) )
  }

  warning("plotTheme: no match found in listed themes for: ", plotKind)
  return(basicTheme)
}

#' getColorsVector
#'
#' @description This function returns a list of colors based on the
#'   [[brewer.pal()]] function
#'
#' @details The colors are taken from the [[brewer.pal.info()]] sets with
#'   `Set1`, `Set2`, `Set3` placed first.
#'
#' @param numNeededColors The number of returned colors
#'
#' @returns an array of `RGB` colors of the wanted size
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @export
#'
#' @examples
#' colorsVector <- getColorsVector(17)
#'
#' @rdname getColorsVector
#'
getColorsVector <- function(numNeededColors) {
  qualColPalets <- brewer.pal.info[brewer.pal.info[["category"]] == "qual", ]
  numColPalets <- nrow(qualColPalets)

  qualColPalets <- rbind(qualColPalets[(numColPalets - 2L):numColPalets, ],
                         qualColPalets[1L:(numColPalets - 3L), ])

  colVector <- unlist(mapply(brewer.pal, qualColPalets[["maxcolors"]],
                             rownames(qualColPalets)))

  rm(qualColPalets, numColPalets)

  assert_that(numNeededColors <= length(colVector),
              msg = paste("Needed more colors than the number",
                          "of possible supported colors:", length(colVector)))

  return(colVector[1L:numNeededColors])
}

#----------------- legacy functions --------------------

#' Handle symmetric matrix <-> vector conversions
#'
#' @description Converts a symmetric matrix into a compacted symmetric matrix
#'   and vice-versa.
#'
#' @details This is a legacy function related to old `scCOTAN` objects. Use the
#'   more appropriate `Matrix::dspMatrix` type for similar functionality.
#'
#'   `mat2vec_rfast`will forcibly make its argument symmetric.
#'
#' @param x a `list` formed by two arrays: `genes` with the unique gene names
#'   and `values` with all the values.
#' @param genes an array with all wanted genes or the string `"all"`. When equal
#'   to `"all"` (the default), it recreates the entire matrix.
#' @param mat a square (possibly symmetric) matrix with all genes as row and
#'   column names.
#'
#' @returns `vec2mat_rfast` returns the reconstructed symmetric matrix
#'
#'   `mat2vec_rfast` a `list` formed by two arrays:
#'   * `genes` with the unique gene names,
#'   * `values` with all the values.
#'
#' @examples
#' v <- list("genes" = paste0("gene_", c(1:9)), "values" = c(1:45))
#'
#' M <- vec2mat_rfast(v)
#' all.equal(rownames(M), v[["genes"]])
#' all.equal(colnames(M), v[["genes"]])
#'
#' genes <- paste0("gene_", sample.int(ncol(M), 3))
#'
#' m <- vec2mat_rfast(v, genes)
#' all.equal(rownames(m), v[["genes"]])
#' all.equal(colnames(m), genes)
#'
#' v2 <- mat2vec_rfast(M)
#' all.equal(v, v2)
#'

#' @importFrom Rfast lower_tri.assign
#' @importFrom Rfast upper_tri.assign
#' @importFrom Rfast upper_tri
#' @importFrom Rfast transpose
#'
#' @export
#'
#' @rdname LegacyFastSymmMatrix
#'
vec2mat_rfast <- function(x, genes = "all") {
  if (!isa(x, "list") || !identical(names(x), c("genes", "values"))) {
    stop("Passed 'x' argument is not a list ",
         "with 2 elements 'genes' and 'values'")
  }

  numGenes <- length(x[["genes"]])
  if (genes[[1L]] == "all") {
    m <- matrix(0.0, numGenes, numGenes)

    m <- lower_tri.assign(m, diag = TRUE, v = x[["values"]])
    m <- upper_tri.assign(m, v = upper_tri(transpose(m)))
    rownames(m) <- x[["genes"]]
    colnames(m) <- x[["genes"]]
  } else {
    m <- matrix(0.0, nrow = numGenes, ncol = length(genes))
    rownames(m) <- x[["genes"]]
    colnames(m) <- genes

    for (posGene in match(genes, x[["genes"]])) {
      tempArray <- x[["values"]][posGene]
      p <- 1L
      l <- numGenes
      s <- posGene
      while (p < (posGene)) {
        l <- l - 1L
        s <- s + l
        tempArray <- c(tempArray, x[["values"]][s])
        p <- p + 1L
      }
      if (posGene < numGenes) {
        #linear part
        startReadingPos <- 1L
        i <- 1L
        while (i < (posGene)) {
          startReadingPos <- startReadingPos + (numGenes - (i - 1L))
          i <- i + 1L
        }
        startReadingPos <- startReadingPos + 1L

        endReadingPos <- 0L
        for (i in c(0L:(posGene - 1L))) {
          endReadingPos <- endReadingPos + (numGenes - i)
        }

        tempArray <-
          c(tempArray, x[["values"]][startReadingPos:endReadingPos])
      }
      m[, x[["genes"]][posGene]] <- tempArray
    }
  }
  return(m)
}


#'
#' @importFrom Rfast lower_tri
#'
#' @export
#'
#' @rdname LegacyFastSymmMatrix
#'
mat2vec_rfast <- function(mat) {
  mat <- as.matrix(mat)

  if (!dim(mat)[[1L]] == dim(mat)[[2L]]) {
    stop("The matrix is not square!")
  }

  v <- lower_tri(mat, diag = TRUE)
  return(list("genes" = rownames(mat), "values" = v))
}
