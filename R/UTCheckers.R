# ---- checkObjIsUniform ----

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
#' @rdname UniformTranscriptCheckers
#'
setMethod(
  "checkObjIsUniform",
  signature(currentC  = "SimpleGDIUniformityCheck",
            previousC = "ANY",
            objCOTAN  = "ANY"),
  function(currentC, previousC = NULL, objCOTAN = NULL) {
    invisible(validObject(currentC))
    assert_that(!currentC@check@isCheckAbove)
    assert_that(currentC@check@maxRankBeyond == 0L)
    assert_that(is.null(previousC) ||
                  class(previousC) == "SimpleGDIUniformityCheck")
    assert_that(is.null(objCOTAN) || class(objCOTAN) == "COTAN")

    previousCheckIsSufficient <- FALSE
    currentC@clusterSize <- ifelse(is.null(objCOTAN), 0L, getNumCells(objCOTAN))

    cCheck <- currentC@check

    if(!is.null(previousC)) {
      invisible(validObject(previousC))
      assert_that(!previousC@check@isCheckAbove)
      assert_that(previousC@check@maxRankBeyond == 0L)

      # in this case we assume previousC is the result of a previous check,
      # maybe with different thresholds: if possible we try to avoid using the
      # GDI to determine the check result

      if (currentC@clusterSize == 0L) {
        currentC@clusterSize <- previousC@clusterSize
      } else {
        assert_that(currentC@clusterSize == previousC@clusterSize)
      }

      pCheck <- previousC@check

      if (is.finite(pCheck@fractionBeyond) &&
          (pCheck@GDIThreshold == cCheck@GDIThreshold)) {
        cCheck@fractionBeyond <- pCheck@fractionBeyond
        previousCheckIsSufficient <- TRUE
      }

      if (is.finite(pCheck@quantileAtRatio) &&
          (pCheck@maxRatioBeyond == cCheck@maxRatioBeyond)) {
        cCheck@quantileAtRatio <- pCheck@quantileAtRatio
        previousCheckIsSufficient <- TRUE
      }

      # if neither threshold match previous check we cannot avoid re-doing the
      # check via the GDI, so fall-back from here
    }

    # run the check via pre-calculated GDI if possible
    if (!previousCheckIsSufficient && !is.null(objCOTAN)
        && getNumCells(objCOTAN) == currentC@clusterSize
        && isCoexAvailable(objCOTAN)) {
      gdi <- getGDI(objCOTAN)
      if (is_empty(gdi)) {
        # if GDI was not stored recalculate it now
        gdi <- getColumnFromDF(calculateGDI(objCOTAN), colName = "GDI")
      }

      assert_that(!is_empty(gdi))

      cCheck@quantileAtRatio <-
        quantile(gdi, probs = 1.0 - cCheck@maxRatioBeyond, names = FALSE)

      cCheck@fractionBeyond <-
        sum(gdi >= cCheck@GDIThreshold) / length(gdi)
    }

    if (is.finite(cCheck@fractionBeyond)) {
      currentC@isUniform <-
        (cCheck@fractionBeyond <= cCheck@maxRatioBeyond)
    } else if (is.finite(cCheck@quantileAtRatio)) {
      currentC@isUniform <-
        (cCheck@quantileAtRatio <= cCheck@GDIThreshold)
    } else {
      # signal no check result can be established
      currentC@clusterSize <- 0L
    }

    currentC@check <- cCheck

    return(currentC)
  }
)

# ------  serialization ------

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
#' @rdname UniformTranscriptCheckers
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
    checkerAsList <- list()
    for (name in memberNames) {
      member <- slot(checker, name)
      membersList <- list()
      if (class(member) == "GDICheck") {
        membersList <- lapply(slotNames(member),
                              function(subName) slot(member, subName))
        names(membersList) <- paste0(name, ".", slotNames(member))
      } else {
        membersList <- list(member)
        names(membersList)[[1]] <- name
      }
      checkerAsList <- append(checkerAsList, membersList)
    }

    df <- rbind(df, checkerAsList)

    if (nrow(df) == 1L) {
      colnames(df) <- names(checkerAsList)
    } else {
      assert_that(identical(colnames(df), names(checkerAsList)))
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
#' @rdname UniformTranscriptCheckers
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


# ------  getCheckerThreshold ------

#' @aliases getCheckerThreshold
#'
#' @description `getCheckerThreshold()` extracts the main `GDI` threshold
#'   from the given checker object
#'
#' @param checker An checker object that defines how to check for *uniform
#'   transcript*. It is derived from [BaseUniformityCheck-class]
#'
#' @returns `getCheckerThreshold()` returns the appropriate member of the
#'   checker object representing the main `GDI` threshold
#'
#' @export
#'
#' @rdname UniformTranscriptCheckers

setMethod(
  "getCheckerThreshold",
  signature(checker = "SimpleGDIUniformityCheck"),
  function(checker) {
    invisible(validObject(checker))
    return(checker@check@GDIThreshold)
  })

setMethod(
  "getCheckerThreshold",
  signature(checker = "AdvancedGDIUniformityCheck"),
  function(checker) {
    invisible(validObject(checker))
    return(checker@firstCheck@GDIThreshold)
  })


# ------  calculateThresholdShiftToUniformity ------

#' @aliases calculateThresholdShiftToUniformity
#'
#' @description `calculateThresholdShiftToUniformity()` calculates by how much
#'   the `GDI` thresholds in the given checker must be increased in order to
#'   have that the relevant cluster is deemed **uniform transcript**
#'
#' @param checker An checker object that defines how to check for *uniform
#'   transcript*. It is derived from [BaseUniformityCheck-class]
#'
#' @returns `calculateThresholdShiftToUniformity()` returns the positive shift
#'   that would make the `@isUniform` slot `TRUE` in the checker. It returns
#'   zero if the result is already `TRUE` and `NaN` in case no such shift can
#'   exist (e.g. the check have been not done yet)
#'
#' @export
#'
#' @rdname UniformTranscriptCheckers

setMethod(
  "calculateThresholdShiftToUniformity",
  signature(checker = "SimpleGDIUniformityCheck"),
  function(checker) {
    invisible(validObject(checker))
    if (checker@clusterSize == 0L ||
        !is.finite(checker@check@quantileAtRatio)) {
      return(NaN)
    }
    return(checker@check@quantileAtRatio - checker@check@GDIThreshold)
  })

# setMethod(
#   "calculateThresholdShiftToUniformity",
#   signature(checker = "AdvancedGDIUniformityCheck"),
#   function(checker) {
#     invisible(validObject(checker))
#     if (checker@clusterSize == 0L || TRUE) {
#       return(NaN)
#     }
#     # delta1 <- GDI_q95 - 1.297
#     # if (GDI_q98 > 1.307 + max(0.0, delta1)) {
#     #   return(delta1)
#     # }
#     # delta2 <- max(delta1, GDI_2nd - 1.4)
#     # return(delta2)
#   })

# ------  shiftCheckerThresholds ------


# equivFractionAbove <- function(GDIThreshold, ratioAboveThreshold,
#                                ratioQuantile, fractionAbove,
#                                usedGDIThreshold, usedRatioAbove) {
#   assert_that(!is.na(GDIThreshold), !is.na(ratioAboveThreshold),
#               !is.na(usedGDIThreshold), !is.na(usedRatioAbove),
#               GDIThreshold >= 0.0, ratioAboveThreshold >= 0.0,
#               ratioAboveThreshold <= 1.0, msg = "wrong thresholds passed in")
#   if (GDIThreshold == usedGDIThreshold) {
#     return(fractionAbove)
#   } else if (ratioAboveThreshold == usedRatioAbove) {
#     # here we assume exponential taper
#     fractionAbove <- max(fractionAbove, 1.0e-4)
#     if (abs(usedGDIThreshold - ratioQuantile) <= 1e-4) {
#       return(NA)
#     }
#     exponent <- (GDIThreshold - usedGDIThreshold) /
#                   (usedGDIThreshold - ratioQuantile)
#     return(fractionAbove * (fractionAbove/usedRatioAbove)^exponent)
#   } else {
#     return(NA)
#   }
# }
