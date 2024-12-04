# ---- checkObjIsUniform ----

#' @aliases checkObjIsUniform
#'
#' @details `checkObjIsUniform()` performs the check whether the given object is
#'   uniform according to the given checker
#'
#' @param currentC the object that defines the method and the threshold to
#'   discriminate whether a *cluster* is *uniform transcript*.
#' @param previousC the optional result object of an already done check
#' @param objCOTAN an optional `COTAN` object
#'
#' @returns a copy of `currentC` with the results of the check. Note that the
#'   slot `clusterSize` will be set to zero if it is not possible to get the
#'   result of the check
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom rlang is_empty
#'
#' @importFrom methods is
#' @importFrom methods validObject
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
    checkerIsOK <- function(checker) {
      invisible(validObject(checker))
      return(!checker@check@isCheckAbove &&
              checker@check@maxRankBeyond  == 0L)
    }

    assert_that(checkerIsOK(currentC))
    assert_that(is.null(previousC) ||
                  is(previousC, "SimpleGDIUniformityCheck"))
    assert_that(is.null(objCOTAN) || is(objCOTAN, "COTAN"))

    previousCheckIsSufficient <- FALSE
    currentC@clusterSize <- ifelse(is.null(objCOTAN), 0L, getNumCells(objCOTAN))

    cCheck <- currentC@check

    if (!is.null(previousC)) {
      assert_that(checkerIsOK(previousC))

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

#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom rlang is_empty
#'
#' @importFrom methods is
#' @importFrom methods validObject
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @rdname UniformTranscriptCheckers

setMethod(
  "checkObjIsUniform",
  signature(currentC  = "AdvancedGDIUniformityCheck",
            previousC = "ANY",
            objCOTAN  = "ANY"),
  function(currentC, previousC = NULL, objCOTAN = NULL) {
    checkerIsOK <- function(checker) {
      invisible(validObject(checker))
      return(!checker@firstCheck@isCheckAbove  &&
              checker@secondCheck@isCheckAbove &&
             !checker@thirdCheck@isCheckAbove  &&
             !checker@fourthCheck@isCheckAbove &&
             checker@firstCheck@maxRankBeyond  == 0L &&
             checker@secondCheck@maxRankBeyond == 0L &&
             checker@thirdCheck@maxRankBeyond  == 0L &&
             !is.finite(checker@fourthCheck@maxRatioBeyond))
    }
    assert_that(checkerIsOK(currentC))
    assert_that(is.null(previousC) ||
                  is(previousC, "AdvancedGDIUniformityCheck"))
    assert_that(is.null(objCOTAN) || is(objCOTAN, "COTAN"))

    previousCheckIsSufficient <- FALSE
    currentC@clusterSize <- ifelse(is.null(objCOTAN), 0L, getNumCells(objCOTAN))

    cCheck1 <- currentC@firstCheck
    cCheck2 <- currentC@secondCheck
    cCheck3 <- currentC@thirdCheck
    cCheck4 <- currentC@fourthCheck

    if (!is.null(previousC)) {
      assert_that(checkerIsOK(previousC))

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
      pCheck4 <- previousC@fourthCheck

      if (is.finite(pCheck1@fractionBeyond) &&
          is.finite(pCheck2@fractionBeyond) &&
          is.finite(pCheck3@fractionBeyond) &&
          pCheck4@thresholdRank > 0L &&
          pCheck1@GDIThreshold == cCheck1@GDIThreshold &&
          pCheck2@GDIThreshold == cCheck2@GDIThreshold &&
          pCheck3@GDIThreshold == cCheck3@GDIThreshold &&
          pCheck4@GDIThreshold == cCheck4@GDIThreshold) {
        cCheck1@fractionBeyond <- pCheck1@fractionBeyond
        cCheck2@fractionBeyond <- pCheck2@fractionBeyond
        cCheck3@fractionBeyond <- pCheck3@fractionBeyond
        cCheck4@thresholdRank  <- pCheck4@thresholdRank

        previousCheckIsSufficient <- TRUE
      }

      if (is.finite(pCheck1@quantileAtRatio) &&
          is.finite(pCheck2@quantileAtRatio) &&
          is.finite(pCheck3@quantileAtRatio) &&
          is.finite(pCheck4@quantileAtRank)  &&
          pCheck1@maxRatioBeyond == cCheck1@maxRatioBeyond &&
          pCheck2@maxRatioBeyond == cCheck2@maxRatioBeyond &&
          pCheck3@maxRatioBeyond == cCheck3@maxRatioBeyond &&
          pCheck4@maxRankBeyond  == cCheck4@maxRankBeyond) {
        cCheck1@quantileAtRatio <- pCheck1@quantileAtRatio
        cCheck2@quantileAtRatio <- pCheck2@quantileAtRatio
        cCheck3@quantileAtRatio <- pCheck3@quantileAtRatio
        cCheck4@quantileAtRank  <- pCheck4@quantileAtRank

        previousCheckIsSufficient <- TRUE
      }

      # if neither threshold match previous check we cannot avoid re-doing the
      # check via the GDI, so fall-back from here
    }
    # run the check via pre-calculated GDI if possible
    if (!previousCheckIsSufficient && !is.null(objCOTAN) &&
        getNumCells(objCOTAN) == currentC@clusterSize &&
        isCoexAvailable(objCOTAN)) {
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

      cCheck3@quantileAtRatio <-
        quantile(gdi, probs = 1.0 - cCheck3@maxRatioBeyond, names = FALSE)

      cCheck3@fractionBeyond <-
        sum(gdi >= cCheck3@GDIThreshold) / length(gdi)

      cCheck4@quantileAtRank <-
        sort(gdi, decreasing = TRUE)[[cCheck4@maxRankBeyond]]

      cCheck4@thresholdRank <- sum(gdi >= cCheck4@GDIThreshold)
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
    if (is.finite(cCheck3@fractionBeyond)) {
      thirdCheckOK <-
        (cCheck3@fractionBeyond <= cCheck3@maxRatioBeyond)
    } else if (is.finite(cCheck3@quantileAtRatio)) {
      thirdCheckOK <-
        (cCheck3@quantileAtRatio <= cCheck3@GDIThreshold)
    } else {
      # signal no check result can be established
      currentC@clusterSize <- 0L
    }

    fourthCheckOK <- FALSE
    if (cCheck4@thresholdRank > 0L) {
      fourthCheckOK <-
        (cCheck4@thresholdRank < cCheck4@maxRankBeyond)
    } else if (is.finite(cCheck4@quantileAtRank)) {
      fourthCheckOK <-
        (cCheck4@quantileAtRank <= cCheck4@GDIThreshold)
    } else {
      # signal no check result can be established
      currentC@clusterSize <- 0L
    }

    currentC@isUniform <- firstCheckOK &&
      ((secondCheckOK && thirdCheckOK) || fourthCheckOK)

    currentC@firstCheck  <- cCheck1
    currentC@secondCheck <- cCheck2
    currentC@thirdCheck  <- cCheck3
    currentC@fourthCheck <- cCheck4

    return(currentC)
  })


# ------  serialization ------

#' @details `checkersToDF()` converts a `list` of checkers (i.e. objects that
#'   derive from `BaseUniformityCheck`) into a `data.frame` with the values of
#'   the members
#'
#' @param checkers a `list` of objects that defines the method, the thresholds
#'   and the results of the checks to discriminate whether a *cluster* is deemed
#'   *uniform transcript*.
#'
#' @returns a `data.frame` with col-names being the member names and row-names
#'   the names attached to each checker
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_list
#'
#' @importFrom methods is
#' @importFrom methods slot
#' @importFrom methods slot<-
#' @importFrom methods slotNames
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
      if (is(member, "GDICheck")) {
        membersList <- lapply(slotNames(member),
                              function(subName) slot(member, subName))
        names(membersList) <- paste0(name, ".", slotNames(member))
      } else {
        membersList <- list(member)
        names(membersList)[[1L]] <- name
      }
      checkerAsList <- append(checkerAsList, membersList)
    }

    # Convert checkerAsList to a data frame row
    checkerDF <- as.data.frame(checkerAsList)

    # Bind the data frame row to the main data frame
    if (is_empty(df)) {
      df <- checkerDF
    } else {
      assert_that(identical(colnames(df), colnames(checkerDF)))
      df <- rbind(df, checkerDF)
    }
  }

  rownames(df) <- names(checkers)

  return(df)
}


#' @details `dfToCheckers()` converts a `data.frame` of checkers values into an
#'   array of checkers ensuring given `data.frame` is compatible with member
#'   types
#'
#' @param df a `data.frame` with col-names being the member names and row-names
#'   the names attached to each checker
#' @param checkerClass the type of the checker to be reconstructed from the
#'   given `data.frame`
#'
#' @returns `dfToCheckers()` returns a `list` of checkers of the requested type,
#'   each created from one of `data.frame` rows
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_empty
#' @importFrom rlang is_character
#'
#' @importFrom stringr str_split_i
#' @importFrom stringr fixed
#'
#' @importFrom methods new
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

  assert_that(is_character(checkerClass), !isEmptyName(checkerClass),
              msg = "A valid checker type must be given")

  checker <- new(checkerClass)
  assert_that(is(checker, "BaseUniformityCheck"))

  memberNames <- slotNames(checker)
  dfRootNames <- str_split_i(colnames(df), fixed("."), 1L)
  assert_that(all(memberNames %in% dfRootNames),
              msg = paste0("Given checkerClass [", class(checker), "] is ",
                           "inconsistent with the data.frame column names"))

  for (r in seq_len(nrow(df))) {
    checker <- new(checkerClass)

    # Assign values from the data.frame row to the corresponding slots
    for (name in memberNames) {
      member <- slot(checker, name)
      if (is(member, "GDICheck")) {
        subNames <- slotNames(member)
        if (TRUE) {
          relCols <- dfRootNames == name
          dfSubNames <- str_split_i(colnames(df)[relCols], fixed("."), 2L)
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
#' @importFrom methods validObject
#'
#' @rdname UniformTranscriptCheckers
#'
setMethod(
  "getCheckerThreshold",
  signature(checker = "SimpleGDIUniformityCheck"),
  function(checker) {
    invisible(validObject(checker))
    return(checker@check@GDIThreshold)
  })

#' @export
#'
#' @importFrom methods validObject
#'
#' @rdname UniformTranscriptCheckers
#'
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
#' @importFrom methods validObject
#'
#' @rdname UniformTranscriptCheckers
#'
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

#' @export
#'
#' @importFrom methods validObject
#'
#' @rdname UniformTranscriptCheckers
#'
setMethod(
  "calculateThresholdShiftToUniformity",
  signature(checker = "AdvancedGDIUniformityCheck"),
  function(checker) {
    invisible(validObject(checker))
    if (checker@clusterSize == 0L ||
        !is.finite(checker@firstCheck@quantileAtRatio)  ||
        !is.finite(checker@secondCheck@quantileAtRatio) ||
        !is.finite(checker@thirdCheck@quantileAtRatio) ||
        !is.finite(checker@fourthCheck@quantileAtRank)) {
      return(NaN)
    }
    delta1 <- checker@firstCheck@quantileAtRatio -
                checker@firstCheck@GDIThreshold
    delta2 <- checker@secondCheck@quantileAtRatio -
                checker@secondCheck@GDIThreshold
    delta3 <- checker@thirdCheck@quantileAtRatio -
                checker@thirdCheck@GDIThreshold
    delta4 <- checker@fourthCheck@quantileAtRank -
                checker@fourthCheck@GDIThreshold

    # delta14 is always sufficient to make the cluster UT,
    # if delta13 is lower and delta2 collaborates then delta13 can be used.
    # We use delta13^+ against delta2 as we never lower
    # the secondCheck@GDIThreshold!

    delta13 <- max(delta1, delta3)
    delta14 <- max(delta1, delta4)

    if (delta13 < delta14 && delta2 > max(0.0, delta13)) {
      return(delta13)
    } else {
      return(delta14)
    }
  })

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
#' @importFrom methods new
#' @importFrom methods validObject
#'
#' @rdname UniformTranscriptCheckers
#'
setMethod(
  "shiftCheckerThresholds",
  signature(checker = "SimpleGDIUniformityCheck",
            shift = "numeric"),
  function(checker, shift) {
    invisible(validObject(checker))
    assert_that(is.finite(shift), shift >= 0.0,
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
#' @importFrom methods new
#' @importFrom methods validObject
#'
#' @rdname UniformTranscriptCheckers
#'
setMethod(
  "shiftCheckerThresholds",
  signature(checker = "AdvancedGDIUniformityCheck",
            shift = "numeric"),
  function(checker, shift) {
    invisible(validObject(checker))
    assert_that(is.finite(shift), shift >= 0.0,
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
                     isCheckAbove   = checker@thirdCheck@isCheckAbove,
                     GDIThreshold   = checker@thirdCheck@GDIThreshold + shift,
                     maxRatioBeyond = checker@thirdCheck@maxRatioBeyond,
                     maxRankBeyond  = checker@thirdCheck@maxRankBeyond),
               fourthCheck =
                 new("GDICheck",
                     isCheckAbove   = checker@fourthCheck@isCheckAbove,
                     GDIThreshold   = checker@fourthCheck@GDIThreshold + shift,
                     maxRatioBeyond = checker@fourthCheck@maxRatioBeyond,
                     maxRankBeyond  = checker@fourthCheck@maxRankBeyond)))
  })
