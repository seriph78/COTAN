#' @aliases checkObjIsUniform
#'
#' @details `checkObjIsUniform()` performs the check whether the given object is
#'   uniform according to the given checker
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @rdname UT_Check
#'
setMethod(
  "checkObjIsUniform",
  signature(objCOTAN = "COTAN",
            currentC = "SimpleGDIUniformityCheck",
            previousC = "ANY"),
  function(objCOTAN, currentC, previousC = NULL) {
    invisible(validObject(currentC))

    previousCheckIsSufficient <- FALSE
    currentC@clusterSize <- getNumCells(objCOTAN)

    if(!is.null(previousC)) {
      invisible(validObject(previousC))
      # in this case we assume previousC is the result of a previous check,
      # maybe with different thresholds: if possible we try to avoid using the
      # GDI to determine the check result

      assert_that(currentC@clusterSize == previousC@clusterSize)

      if (is.finite(previousC@fractionAboveThreshold) &&
          (previousC@GDIThreshold == currentC@GDIThreshold)) {
        currentC@fractionAboveThreshold <- previousC@fractionAboveThreshold
        previousCheckIsSufficient <- TRUE
      }

      if (is.finite(previousC@quantileAtRatio) &&
          (previousC@ratioAboveThreshold == currentC@ratioAboveThreshold)) {

        currentC@quantileAtRatio <- previousC@quantileAtRatio
        previousCheckIsSufficient <- TRUE
      }

      # if neither threshold match previous check we cannot avoid re-doing the
      # check via the GDI, so fall-back from here
    }

    # run the check via pre-calculated GDI
    if (!previousCheckIsSufficient) {

      gdi <- getGDI(objCOTAN)
      if (is_empty(gdi)) {
        # if GDI was not stored recalculate it now
        gdi <- getColumnFromDF(calculateGDI(objCOTAN), "GDI")
      }

      assert_that(!is_empty(gdi))

      currentC@quantileAtRatio <-
        quantile(gdi, probs = 1.0 - currentC@ratioAboveThreshold, names = FALSE)

      currentC@fractionAboveThreshold <-
        sum(gdi >= currentC@GDIThreshold) / length(gdi)
    }

    if (is.finite(currentC@fractionAboveThreshold)) {
      currentC@isUniform <-
        (currentC@fractionAboveThreshold <= currentC@ratioAboveThreshold)
    } else {
      assert_that(is.finite(currentC@quantileAtRatio))
      currentC@isUniform <-
        (currentC@quatileAtRatio <= currentC@GDIThreshold)
    }

    return(currentC)
  }
)


#' @details `checkersToDF()` converts a `list` of checkers (i.e. objects that
#'   derive from `BaseUniformityCheck`) into a `data.frame` with the values of
#'   the members
#'
#' @returns a `data.frame` with col-names being the member names and row-names
#'   the names attached to each checker
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_list
#'
#' @export
#'
#' @rdname UT_Check
#'
checkersToDF <- function(checkers) {
  if (!is_list(checkers)) {
    checkers <- list(checkers)
  }
  df <- data.frame()
  for (checker in checkers) {
    # check argument is a checker
    assert_that(is(checker, "BaseUniformityCheck"))

    memberNames <- slotNames(checker)
    checkerAsList <- lapply(memberNames, function(name) slot(checker, name))

    df <- rbind(df, checkerAsList)

    if (nrow(df) == 1L) {
      colnames(df) <- memberNames
    } else {
      assert_that(identical(colnames(df), memberNames))
    }
  }

  rownames(df) <- names(checkers)

  return(df)
}


#' @details `dfToCheckers()` converts a `data.frame` of checkers values into an
#'   array of checkers ensuring given `data.frame` is compatible with member
#'   types
#'
#' @returns `dfToCheckers()` returns a `list` of checkers of the requested type,
#'   each created from one of `data.frame` rows
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_character
#'
#' @export
#'
#' @rdname UT_Check
#'
dfToCheckers <- function(df, checkerClass) {
  checkers <- list()

  if (is_empty(df)) {
    return(checkers)
  }

  assert_that(is_character(checkerClass) && !isEmptyName(checkerClass),
              msg = "A valid checker type must be given")

  checker <- new(checkerClass)
  assert_that(is(checker, "BaseUniformityCheck"))

  memberNames <- slotNames(checker)
  assert_that(all(memberNames %in% colnames(df)),
              msg = paste0("Given checkerClass [", class(checker), "] is ",
                           "inconsistent with the data.frame column names"))

  for (r in seq_len(nrow(df))) {
    checker <- new(checkerClass)

    # Assign values from the data.frame row to the corresponding slots
    for (name in memberNames) {
      slot(checker, name) <- df[r, name]
    }

    checkers <- append(checkers, checker)
  }

  names(checkers) <- rownames(df)

  return(checkers)
}


#' @description `retrieveMainGDIThreshold()` extracts the main GDI threshold
#'   from the given checker class
#'
#' @param checker An checker object that defines how to check for *uniform
#'   transcript*. It is derived from [BaseUniformityCheck-class]
#'
#' @returns `retrieveMainGDIThreshold()` returns the appropriate member of the
#'   checker class or falls-back to 1.3
#'
# #' @export
#'
#' @rdname UT_Check

retrieveMainGDIThreshold <- function(checker) {
  if (is(checker, "SimpleGDIUniformityCheck")) {
    return(checker@GDIThreshold)
  } else if (is(checker, "AdvancedGDIUniformityCheck")) {
    return(checker@lowCheckThreshold)
  } else {
    return(1.4)
  }
}

#
# #'
# #' @details `isClusterUniform()` takes in the current thresholds and used them
# #'   to check whether the calculated cluster parameters are sufficient to
# #'   determine whether the cluster is **uniform** and in the positive scenario
# #'   the corresponding answer
# #'
# #' @param GDIThreshold the threshold level that discriminates uniform
# #'   *clusters*. It defaults to \eqn{1.43}
# #' @param ratioAboveThreshold the fraction of genes allowed to be above the
# #'   `GDIThreshold`. It defaults to \eqn{1\%}
# #' @param ratioQuantile the `GDI` quantile corresponding to the `usedRatioAbove`
# #' @param fractionAbove the fraction of genes above the `usedGDIThreshold`
# #' @param usedGDIThreshold the threshold level actually used to calculate fourth
# #'   argument
# #' @param usedRatioAbove the fraction of genes actually used to calculate the
# #'   third argument
# #'
# #' @returns a single `Boolean` value when it is possible to decide the answer
# #'   with the given information and `NA` otherwise
# #'
# #' @importFrom assertthat assert_that
# #'
# #' @rdname UniformClusters
# #'
#
# isClusterUniform <- function(GDIThreshold, ratioAboveThreshold,
#                              ratioQuantile, fractionAbove,
#                              usedGDIThreshold, usedRatioAbove) {
#   assert_that(!is.na(GDIThreshold), !is.na(ratioAboveThreshold),
#               !is.na(usedGDIThreshold), !is.na(usedRatioAbove),
#               GDIThreshold >= 0.0, ratioAboveThreshold >= 0.0,
#               ratioAboveThreshold <= 1.0, msg = "wrong thresholds passed in")
#
#   if (!is.na(fractionAbove) && GDIThreshold == usedGDIThreshold) {
#     return(fractionAbove <= ratioAboveThreshold)
#   } else if (!is.na(ratioQuantile) && ratioAboveThreshold == usedRatioAbove) {
#     return(ratioQuantile < GDIThreshold)
#   } else {
#     return(NA)
#   }
# }
#


#' @details `checkClusterUniformity()` takes a `COTAN` object and a cells'
#'   *cluster* and checks whether the latter is **uniform** by `GDI`. The
#'   function runs `COTAN` to check whether the `GDI` is lower than the given
#'   `GDIThreshold` (1.43) for all but at the most `ratioAboveThreshold`
#'   (\eqn{1\%}) genes. If the `GDI` results to be too high for too many genes,
#'   the *cluster* is deemed
#'   **non-uniform**.
#'
#' @param objCOTAN a `COTAN` object
#' @param clusterName the tag of the *cluster*
#' @param cells the cells belonging to the *cluster*
#' @param checker the object that defines the method and the threshold to
#'   discriminate whether a *cluster* is *uniform transcript*. See [UT_Check]
#'   for more details
#' @param cores number of cores to use. Default is 1.
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
#' @param saveObj Boolean flag; when `TRUE` saves intermediate analyses and
#'   plots to file(s)
#' @param outDir an existing directory for the analysis output. The effective
#'   output will be paced in a sub-folder.
#'
#' @returns `checkClusterUniformity` returns a checker object of the same type
#'   as the input one, that contains both threshold and results of the check:
#'   see [UT_Check] for more details
#'
#' @importFrom utils head
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @importFrom withr local_options
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname UniformClusters
#'

checkClusterUniformity <- function(
    objCOTAN,
    clusterName,
    cells,
    checker,
    cores = 1L,
    optimizeForSpeed = TRUE,
    deviceStr = "cuda",
    saveObj = TRUE,
    outDir = ".") {
  # handle legacy usage
  cellsToDrop <- getCells(objCOTAN)[!getCells(objCOTAN) %in% cells]

  objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)

  objCOTAN <- proceedToCoex(objCOTAN, cores = cores,
                            optimizeForSpeed = optimizeForSpeed,
                            deviceStr = deviceStr, saveObj = FALSE)
  gc()

  checker@clusterSize <- getNumCells(objCOTAN)

  logThis(paste0("Checking uniformity for the cluster '", clusterName,
                 "' with ", checker@clusterSize, " cells"), logLevel = 2L)

  GDIData <- calculateGDI(objCOTAN)
  objCOTAN <- storeGDI(objCOTAN, GDIData)

  # Plots
  if (isTRUE(saveObj) && !dir.exists(outDir)) {
    saveObj <- FALSE
    warning(paste("Asked to save check results,",
                  "but given output folder does not exist"))
  }

  if (isTRUE(saveObj)) tryCatch({
      # this will be restored to previous value on function exit
      local_options(list(ggrepel.max.overlaps = Inf))

      pdf(file.path(outDir, paste0("cluster_", clusterName, "_plots.pdf")))

      c(..., nuPlot, zoomedNuPlot) %<-%
        cleanPlots(objCOTAN, includePCA = FALSE)

      genesToLabel <-
        head(rownames(GDIData[order(GDIData[["GDI"]],
                                    decreasing = TRUE), ]), n = 10L)
      gdiPlot <- GDIPlot(objCOTAN, GDIIn = GDIData,
                         GDIThreshold = retrieveMainGDIThreshold(checker),
                         genes = list("top 10 GDI genes" = genesToLabel))

      plot(nuPlot)
      plot(zoomedNuPlot)
      plot(gdiPlot)

      rm(nuPlot, zoomedNuPlot, gdiPlot)
      dev.off()
    }, error = function(err) {
      logThis(paste("While saving cluster plots", err), logLevel = 0L)
      dev.off()
    }
  )

  checker <- checkObjIsUniform(objCOTAN, currentC = checker, previousC = NULL)
  rm(objCOTAN)
  gc()

  logThis(paste0(
    "Cluster ", clusterName, ", with size ", checker@clusterSize, ", is ",
    ifelse(checker@isUniform, "", "not "), "uniform"), logLevel = 1)

  if (TRUE) {
    dumpDF <- checkersToDF(checker)
    logThis(paste0(colnames(dumpDF), " = ", dumpDF, collapse = ", "),
            logLevel = 3L)
    rm(dumpDF)
  }

  if (isTRUE(saveObj)) tryCatch({
      outFile <- file.path(outDir,
                           paste0(ifelse(checker@isUniform, "", "non-"),
                                  "uniform_cluster_", clusterName, ".csv"))
      write.csv(cells, file = outFile)
    },
    error = function(err) {
      logThis(paste("While saving current clusterization", err),
              logLevel = 0L)
    }
  )

  return(checker)
}
