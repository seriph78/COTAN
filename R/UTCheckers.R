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
  })

#' @export
#'
#' @rdname UniformTranscriptCheckers

setMethod(
  "checkObjIsUniform",
  signature(currentC  = "AdvancedGDIUniformityCheck",
            previousC = "ANY",
            objCOTAN  = "ANY"),
  function(currentC, previousC = NULL, objCOTAN = NULL) {
    invisible(validObject(currentC))
    assert_that(!currentC@check@isCheckAbove)
    assert_that(currentC@check@maxRankBeyond == 0L)
    assert_that(is.null(previousC) ||
                  class(previousC) == "AdvancedGDIUniformityCheck")
    assert_that(is.null(objCOTAN) || class(objCOTAN) == "COTAN")

    previousCheckIsSufficient <- FALSE
    currentC@clusterSize <- ifelse(is.null(objCOTAN), 0L, getNumCells(objCOTAN))

    cCheck1 <- currentC@firstCheck
    cCheck2 <- currentC@secondCheck
    cCheck3 <- currentC@thirdCheck

    if(!is.null(previousC)) {
      invisible(validObject(previousC))
      assert_that(!previousC@firstCheck@isCheckAbove  &&
                   previousC@secondCheck@isCheckAbove &&
                  !previousC@thirdCheck@isCheckAbove)
      assert_that(previousC@firstCheck@maxRankBeyond  == 0L &&
                  previousC@secondCheck@maxRankBeyond == 0L &&
                  !is.finite(previousC@thirdCheck@maxRatioBeyond))

      # in this case we assume previousC is the result of a previous check,
      # maybe with different thresholds: if possible we try to avoid using the
      # GDI to determine the check result

      if (currentC@clusterSize == 0L) {
        currentC@clusterSize <- previousC@clusterSize
      } else {
        assert_that(currentC@clusterSize == previousC@clusterSize)
      }

      pCheck1 <- previousC@firstCheck
      pCheck2 <- previousC@secondCheck
      pCheck3 <- previousC@thirdCheck

      if (is.finite(pCheck1@fractionBeyond) &&
          is.finite(pCheck2@fractionBeyond) &&
          pCheck3@thresholdRank > 0L &&
          pCheck1@GDIThreshold == cCheck1@GDIThreshold &&
          pCheck2@GDIThreshold == cCheck2@GDIThreshold &&
          pCheck3@GDIThreshold == cCheck3@GDIThreshold) {
        cCheck1@fractionBeyond <- pCheck1@fractionBeyond
        cCheck2@fractionBeyond <- pCheck2@fractionBeyond
        cCheck3@thresholdRank  <- pCheck3@thresholdRank

        previousCheckIsSufficient <- TRUE
      }

      if (is.finite(pCheck1@quantileAtRatio) &&
          is.finite(pCheck2@quantileAtRatio) &&
          is.finite(pCheck3@quantileAtRank)  &&
          pCheck1@maxRatioBeyond == cCheck1@maxRatioBeyond &&
          pCheck2@maxRatioBeyond == cCheck2@maxRatioBeyond &&
          pCheck3@maxRankBeyond  == cCheck1@maxRankBeyond) {
        cCheck1@quantileAtRatio <- pCheck1@quantileAtRatio
        cCheck2@quantileAtRatio <- pCheck2@quantileAtRatio
        cCheck3@quantileAtRank  <- pCheck3@quantileAtRank

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

      cCheck1@quantileAtRatio <-
        quantile(gdi, probs = 1.0 - cCheck1@maxRatioBeyond, names = FALSE)

      cCheck1@fractionBeyond <-
        sum(gdi >= cCheck1@GDIThreshold) / length(gdi)

      cCheck2@quantileAtRatio <-
        quantile(gdi, probs = 1.0 - cCheck2@maxRatioBeyond, names = FALSE)

      cCheck2@fractionBeyond <-
        sum(gdi >= cCheck2@GDIThreshold) / length(gdi)

      cCheck3@quantileAtRank <- gdi[order(gdi)][cCheck2@maxRankBeyond]

      cCheck3@thresholdRank <- sum(gdi >= cCheck3@GDIThreshold)
    }

    firstCheckOK <- FALSE
    if (is.finite(cCheck1@fractionBeyond)) {
      firstCheckOK <-
        (cCheck1@fractionBeyond <= cCheck1@maxRatioBeyond)
    } else if (is.finite(cCheck1@quantileAtRatio)) {
      firstCheckOK <-
        (cCheck1@quantileAtRatio <= cCheck1@GDIThreshold)
    } else {
      # signal no check result can be established
      currentC@clusterSize <- 0L
    }

    secondCheckOK <- FALSE
    if (is.finite(cCheck2@fractionBeyond)) {
      secondCheckOK <-
        (cCheck2@fractionBeyond >= cCheck2@maxRatioBeyond)
    } else if (is.finite(cCheck2@quantileAtRatio)) {
      secondCheckOK <-
        (cCheck2@quantileAtRatio >= cCheck2@GDIThreshold)
    } else {
      # signal no check result can be established
      currentC@clusterSize <- 0L
    }

    thirdCheckOK <- FALSE
    if (is.finite(cCheck3@thresholdRank)) {
      thirdCheckOK <-
        (cCheck3@thresholdRank <= cCheck3@maxRankBeyond)
    } else if (is.finite(cCheck3@quantileAtRank)) {
      thirdCheckOK <-
        (cCheck3@quantileAtRank <= cCheck3@GDIThreshold)
    } else {
      # signal no check result can be established
      currentC@clusterSize <- 0L
    }

    currentC@isUniform <- firstCheckOK && secondCheckOK && thirdCheckOK

    currentC@firstCheck  <- cCheck1
    currentC@secondCheck <- cCheck2
    currentC@thirdCheck  <- cCheck3

    return(currentC)
  })


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
  dfRootNames <- vapply(str_split(colnames(df), "[.]"),
                        function(a) { return(a[[1]]) },
                        character(1L))
  assert_that(all(memberNames %in% dfRootNames),
              msg = paste0("Given checkerClass [", class(checker), "] is ",
                           "inconsistent with the data.frame column names"))

  for (r in seq_len(nrow(df))) {
    checker <- new(checkerClass)

    # Assign values from the data.frame row to the corresponding slots
    for (name in memberNames) {
      member <- slot(checker, name)
      if (class(member) == "GDICheck") {
        subNames <- slotNames(member)
        if (TRUE) {
          relCols <- dfRootNames == name
          dfSubNames <- vapply(str_split(colnames(df)[relCols], "[.]"),
                               function(a) { return(a[[2]]) },
                               character(1L))
          assert_that(all(subNames %in% dfSubNames),
                      msg = paste0("Given checkerClass [", class(checker),
                                   "] is inconsistent with the data.frame",
                                   "column names"))
          rm(relCols, dfSubNames)
        }
        for (subName in subNames) {
          slot(member, subName) <- df[r, paste0(name, ".", subName)]
        }
      } else {
        member <- df[r, name]
      }

      slot(checker, name) <- member
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

#' @export
#'
#' @rdname UniformTranscriptCheckers

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

# ------  shiftCheckerThresholds ------

#' @aliases shiftCheckerThresholds
#'
#' @description `shiftCheckerThresholds()` returns a new checker object where
#'   the `GDI` thresholds where increased in order to *relax* the conditions to
#'   achieve **uniform transcript**
#'
#' @param checker An checker object that defines how to check for *uniform
#'   transcript*. It is derived from [BaseUniformityCheck-class]
#'
#' @param shift The amount by which to shift the `GDI` thresholds in the checker
#'
#' @returns `shiftCheckerThresholds()` returns a copy of the checker object
#'   where all `GDI` thresholds have been shifted by the same given `shift`
#'   amount
#'
#' @export
#'
#' @rdname UniformTranscriptCheckers

setMethod(
  "shiftCheckerThresholds",
  signature(checker = "SimpleGDIUniformityCheck",
            shift = "numeric"),
  function(checker, shift) {
    invisible(validObject(checker))
    assert_that(is.finite(shift) && shift >= 0.0,
                msg = "Given shift must be a non negative number")

    return(new("SimpleGDIUniformityCheck",
               check = new("GDICheck",
                           isCheckAbove   = checker@check@isCheckAbove,
                           GDIThreshold   = checker@check@GDIThreshold + shift,
                           maxRatioBeyond = checker@check@maxRatioBeyond,
                           maxRankBeyond  = checker@check@maxRankBeyond)))
  })

#' @export
#'
#' @rdname UniformTranscriptCheckers

setMethod(
  "shiftCheckerThresholds",
  signature(checker = "AdvancedGDIUniformityCheck",
            shift = "numeric"),
  function(checker, shift) {
    invisible(validObject(checker))
    assert_that(is.finite(shift) && shift >= 0.0,
                msg = "Given shift must be a non negative number")

    return(new("AdvancedGDIUniformityCheck",
               firstCheck =
                 new("GDICheck",
                     isCheckAbove   = checker@firstCheck@isCheckAbove,
                     GDIThreshold   = checker@firstCheck@GDIThreshold + shift,
                     maxRatioBeyond = checker@firstCheck@maxRatioBeyond,
                     maxRankBeyond  = checker@firstCheck@maxRankBeyond),
               secondCheck =
                 new("GDICheck",
                     isCheckAbove   = checker@secondCheck@isCheckAbove,
                     GDIThreshold   = checker@secondCheck@GDIThreshold + shift,
                     maxRatioBeyond = checker@secondCheck@maxRatioBeyond,
                     maxRankBeyond  = checker@secondCheck@maxRankBeyond),
               thirdCheck =
                 new("GDICheck",
                     isCheckAbove   = checker@thirdCheckk@isCheckAbove,
                     GDIThreshold   = checker@thirdCheckk@GDIThreshold + shift,
                     maxRatioBeyond = checker@thirdCheckk@maxRatioBeyond,
                     maxRankBeyond  = checker@thirdCheckk@maxRankBeyond)))
  })


