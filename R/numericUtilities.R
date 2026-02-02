#'
#' @title Numeric Utilities
#'
#' @description A set of function helper related to the statistical model
#'   underlying the `COTAN` package
#'
#' @name NumericUtilities
NULL

#------------------- negative binomial ----------

#'
#' @details `funProbZero` is a private function that gives the probability that
#'   a sample gene's reads are zero, given the `dispersion` and `mu` parameters.
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
#' @param dispersion the estimated `dispersion` (a \eqn{n}-sized vector)
#' @param mu the `lambda` times `nu` values (a \eqn{n \times m} matrix)
#'
#' @returns the probability `matrix` that a *read count* is identically zero
#'
#' @rdname NumericUtilities
#'
funProbZero <- function(dispersion, mu) {
  disp <- as.numeric(dispersion)
  numGenes <- length(disp)

  assert_that(numGenes != 0L, msg = "empty input given")

  returnArray <- is.null(dim(mu))

  if (returnArray) {
    mu <- matrix(mu, nrow = numGenes)
  }

  # if mu is a matrix
  assert_that(nrow(mu) == numGenes,
              msg = "misaligned dimensions between given `mu` and `dispersion`")

  eps <- 1.0e-5

  mif_idx <- which(disp == -Inf)
  neg_idx <- which(disp <  eps & disp != -Inf)
  pos_idx <- which(disp >= eps & disp != +Inf)
  pif_idx <- which(disp == +Inf)

  out <- disp * mu

  if (length(mif_idx)) {
    out[mif_idx, ] <- -Inf
  }
  if (length(neg_idx)) {
    # out = dispersion * mu - mu
    out[neg_idx, ] <-
      out[neg_idx, , drop = FALSE] - mu[neg_idx, , drop = FALSE]
  }
  if (length(pos_idx)) {
    # out = -(1/disp)*log1p(disp*mu)
    out[pos_idx, ] <-
      -(1.0 / disp[pos_idx]) * log1p(out[pos_idx, , drop = FALSE])
  }
  if (length(pif_idx)) {
    out[pif_idx, ] <- 0.0
  }

  out <- exp(out)

  if (returnArray) {
    out <- as.vector(out)
  }

  return(out)
}


#--------------------- dispersion solvers ----------------

#' @details `dispersionBisection` is a private function for the estimation of
#'   `dispersion` slot of a `COTAN` object via a bisection solver
#'
#' @details The goal is to find a `dispersion` value that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param sumZeros the number of cells that didn't express the gene
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the `dispersion` value
#'
#' @rdname NumericUtilities
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
    } else if (sumZeros == length(nu)) {
      # in case of zero lambda dispersion is irrelevant. We return 1.0
      return(1.0)
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
    logThis(paste("Dispersion bisection: straddling used",
                  iter, "iterations"), logLevel = 3L)

    # once we have found the two bounds to the dispersion value we use bisection
    iter <- 1L
    repeat {
      disp <- (disp1 + disp2) / 2.0

      diff <- sum(funProbZero(disp, mu)) - sumZeros

      if (abs(diff) <= threshold) {
        logThis(paste("Dispersion bisection: used",
                      iter, "iterations"), logLevel = 3L)
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


#' @details `parallelDispersionBisection` is a private function that was invoked
#'   by [estimateDispersionViaSolver()] for the estimation of the `dispersion`
#'   slot of a `COTAN` object via a parallel bisection solver. It is now
#'   deprecated
#'
#' @details The goal is to find a `dispersion` `array` that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param genes names of the relevant genes
#' @param sumZeros the number of cells that didn't express the relevant gene (a
#'   \eqn{n}-sized vector)
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion values
#'
#' @importFrom Rfast rowsums
#'
#' @rdname NumericUtilities
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

    # cannot match exactly zero prob of zeros with finite values
    # so we ignore the rows with no zeros from the solver and return -Inf
    output <- rep(-Inf, length(sumZeros))

    # in case of zero lambda dispersion is irrelevant. We return 1.0
    goodPos <- sumZeros != length(nu)
    output[!goodPos] <- 1.0

    goodPos <- goodPos & sumZeros != 0L

    if (!any(goodPos)) {
      return(output)
    }

    sumZeros <- sumZeros[goodPos]
    lambda <- lambda[goodPos]

    mu <- lambda %o% nu

    # we look for two dispersion values where the first leads to a
    # diffZeros negative and the second positive
    disps1 <- rep(0.0, length(sumZeros))
    diffs1 <- rowsums(funProbZero(disps1, mu)) - sumZeros
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
      diffs2[runPos] <- (rowsums(funProbZero(disps2[runPos],
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
    logThis(paste("parallel dispersion bisection: straddling used up to",
                  iter, "iterations"), logLevel = 3L)

    # once we have found the two bounds to the dispersion value we use bisection
    runNum <- length(diffs1)
    runPos <- rep(TRUE, runNum)
    disps <- disps1
    diffs <- diffs1
    iter <- 1L
    repeat {
      disps[runPos] <- (disps1[runPos] + disps2[runPos]) / 2.0

      diffs[runPos] <- (rowsums(funProbZero(disps[runPos],
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
    logThis(paste("Parallel dispersion bisection: used up to",
                  iter, "iterations"), logLevel = 3L)

    rm(mu)
    gc()

    output[goodPos] <- disps
    return(output)
  }


#' @details `dispersionNewton` is a private function for the estimation of
#'   `dispersion` slot of a `COTAN` object via a Newton-Raphson solver
#'
#' @details The goal is to find a `dispersion` value that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param sumZeros the number of cells that didn't express the gene
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the `dispersion` value
#'
#' @rdname NumericUtilities
#'
dispersionNewton <-
  function(sumZeros,
           lambda,
           nu,
           threshold = 0.001,
           maxIterations = 100L) {
    if (sumZeros == 0L) {
      # cannot match exactly zero prob of zeros with finite values
      return(-Inf)
    } else if (sumZeros == length(nu)) {
      # in case of zero lambda dispersion is irrelevant. We return 1.0
      return(1.0)
    }
    mu <- lambda * nu

    # we look for two dispersion values where the first leads to a
    # diffZeros negative and the second positive
    disp <- 0.0
    diff <- sum(funProbZero(disp, mu)) - sumZeros
    if (abs(diff) <= threshold) {
      return(disp)
    }

    # we assume error is an increasing function of disp
    dispIsNeg <- sign(diff) >= 0.0
    disp <- ifelse(dispIsNeg, -1.0, 1.0)
    iter <- 1L
    repeat {
      dispMu <- as.matrix(disp * mu)

      if (dispIsNeg) {
        logProbZero <- dispMu - mu
      } else {
        logProbZero <- -1.0 * log1p(dispMu) / disp
      }

      probZero <- exp(logProbZero)

      diff <- sum(probZero) - sumZeros

      if (abs(diff) <= threshold) {
        logThis(paste("Dispersion Newton-Raphson: used up to",
                      iter, "iterations"), logLevel = 3L)
        return(disp)
      }

      if (iter >= maxIterations) {
        stop("Max number of iterations reached ",
             "while finding the solution straddling intervals")
      }
      iter <- iter + 1L

      denom <- ifelse(dispIsNeg,
                      sum(dispMu * probZero),
                      sum((mu / (dispMu + 1.0) + logProbZero) * probZero))

      factor <-
        max((denom + ifelse(dispIsNeg, -1.0, 1.0) * diff) / denom, 1.0e-8)

      disp <- disp * factor
    }
  }




#' @details `parallelDispersionNewton` is a private function invoked by
#'   [parallelDispersionNewton()] for the estimation of the `dispersion` slot
#'   of a `COTAN` object via a parallel Newton-Raphson solver
#'
#' @details The goal is to find a `dispersion` `array` that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param genes names of the relevant genes
#' @param sumZeros the number of cells that didn't express the relevant gene (a
#'   \eqn{n}-sized vector)
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion values
#'
#' @importFrom Rfast rowsums
#'
#' @rdname NumericUtilities
#'
parallelDispersionNewton <-
  function(genes,
           sumZeros,
           lambda,
           nu,
           threshold = 0.001,
           maxIterations = 100L) {
    sumZeros <- sumZeros[genes]
    lambda <- lambda[genes]

    # cannot match exactly zero prob of zeros with finite values
    # so we ignore the rows with no zeros from the solver and return -Inf
    output <- rep(-Inf, length(sumZeros))

    # in case of zero lambda dispersion is irrelevant. We return 1.0
    goodPos <- sumZeros != length(nu)
    output[!goodPos] <- 1.0

    goodPos <- goodPos & sumZeros != 0L

    if (!any(goodPos)) {
      return(output)
    }

    sumZeros <- sumZeros[goodPos]
    lambda <- lambda[goodPos]

    mu <- lambda %o% nu

    # we look for two dispersion values where the first leads to a
    # diffZeros negative and the second positive
    disps <- rep(0.0, length(sumZeros))
    diffs <- rowsums(funProbZero(disps, mu)) - sumZeros
    if (all(abs(diffs) <= threshold)) {
      output[goodPos] <- disps
      return(output)
    }

    # we assume error is an increasing function of the dispersion
    dispsIsNeg <- sign(diffs) >= 0.0
    disps <- ifelse(dispsIsNeg, -1.0, 1.0)
    runPos <- rep(TRUE, length(disps))
    iter <- 1L
    repeat {
      dispsMu <- disps[runPos] * mu[runPos, , drop = FALSE]

      isNeg <- dispsIsNeg[runPos]

      logProbsZero <- matrix(NaN, nrow(dispsMu), ncol(dispsMu))
      logProbsZero[ isNeg, ] <-
        (dispsMu[isNeg, , drop = FALSE] -
           mu[runPos & dispsIsNeg, , drop = FALSE])
      logProbsZero[!isNeg, ] <-
        ((-1.0 / disps[runPos & !dispsIsNeg]) *
           log1p(dispsMu[!isNeg, , drop = FALSE]))

      probsZero <- exp(logProbsZero)

      diffs[runPos] <- rowsums(probsZero) - sumZeros[runPos]

      denoms <- rep(1.0, length(disps))
      denoms[runPos &  dispsIsNeg] <-
        rowsums(dispsMu[isNeg, , drop = FALSE] *
                  probsZero[isNeg, , drop = FALSE])
      denoms[runPos & !dispsIsNeg] <-
        rowsums(((mu[runPos & !dispsIsNeg, , drop = FALSE] /
                    (dispsMu[!isNeg, , drop = FALSE] + 1.0)) +
                   logProbsZero[!isNeg, , drop = FALSE]) *
                  probsZero[!isNeg, , drop = FALSE])

      rm(dispsMu, logProbsZero, probsZero)

      diffs[runPos & dispsIsNeg] <- -1.0 * diffs[runPos & dispsIsNeg]

      runPos <- abs(diffs) > threshold
      if (!any(runPos)) {
        break
      }

      if (iter >= maxIterations) {
        stop("Max number of iterations reached while finding the solutions")
      }
      iter <- iter + 1L

      factors <- pmax((denoms + diffs) / denoms, 1e-8)

      disps[runPos] <- disps[runPos] * factors[runPos]
    }
    logThis(paste("Parallel dispersion Newton-Raphson: used up to",
                  iter, "iterations"), logLevel = 3L)
    rm(mu)
    gc()

    output[goodPos] <- disps
    return(output)
  }


#------------------------- nu solvers ---------------------

#' @details `nuBisection` is a private function for the estimation of `nu` slot
#'   of a `COTAN` object via a bisection solver
#'
#' @details The goal is to find a `nu` value that reduces to zero the difference
#'   between the number of estimated and counted zeros
#'
#' @param sumZeros the number non expressed genes in the cell
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param dispersion the estimated `dispersion` (a \eqn{n}-sized vector)
#' @param initialGuess the initial guess for `nu`
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the nu value
#'
#' @rdname NumericUtilities
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

    # once we have found the two bounds to the dispersion value we use bisection
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


#' @details `parallelNuBisection` is a private function invoked by
#'   [estimateNuBisection()] for the estimation of `nu` slot of a `COTAN` object
#'   via a parallel bisection solver
#'
#' @details The goal is to find a `nu` `array` that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param cells names of the relevant cells
#' @param sumZeros the number of genes not expressed in the relevant cell (a
#'   \eqn{m}-sized vector)
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param dispersion the estimated `dispersion` (a \eqn{n}-sized vector)
#' @param initialGuess the initial guess for `nu`  (a \eqn{m}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion values
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom Rfast colsums
#'
#' @rdname NumericUtilities
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

    diffs1 <- colsums(funProbZero(dispersion, lambda %o% nus1)) - sumZeros
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
      diffs2[runPos] <- (colsums(funProbZero(dispersion,
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


      # nus2 are closer to producing 0
      nus1[runPos] <- nus2[runPos]

      # we double (or half) at each step
      nus2[runPos] <- nus2[runPos] * factors[runPos]
    }

    # once we have found the two bounds to the dispersion value we use bisection
    runNum <- length(diffs1)
    runPos <- rep(TRUE, runNum)
    nus <- nus1
    diffs <- diffs1
    iter <- 1L
    repeat {
      nus[runPos] <- (nus1[runPos] + nus2[runPos]) / 2.0

      diffs[runPos] <- (colsums(funProbZero(dispersion,
                                            lambda %o% nus[runPos])) -
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


#----------------- distance functions --------------------


#' @details `calcDist` is a wrapper function that invokes
#'   [parallelDist::parDist()]: the main goal is to recover and finish the
#'   calculations via a fallback when there is a problem with the main algorithm
#'
#' @param data a matrix or a data.frame of which we want to calculate the
#'   distance between columns
#' @param method type of distance to use. Can be chosen among those supported by
#'   [parallelDist::parDist()]
#' @param diag logical value indicating whether the diagonal of the distance
#'   matrix should be printed by print.dist.upper
#' @param upper logical value indicating whether the upper triangle of the
#'   distance matrix should be printed by print.dist
#'
#' @returns a `dist` object with all distances
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom rlang is_null
#'
#' @importFrom parallelDist parDist
#'
#' @importFrom proxy dist
#'
#' @rdname NumericUtilities
#'

calcDist <- function(data, method, diag = FALSE, upper = FALSE) {
  data <- as.matrix(data)

  ret <- NULL

  sI <- Sys.info()
  if (!grepl("Linux", sI[["sysname"]], ignore.case = TRUE) ||
      !grepl("aarch64", sI[["machine"]], ignore.case = TRUE)) {
    ret <- tryCatch(parDist(data, method = method,
                            diag = diag, upper = upper),
                    error = function(err) {
                      logThis(paste("While calculating distance", err),
                              logLevel = 1L)
                      logThis("Falling back to single threaded algo",
                              logLevel = 1L)
                      return(NULL)
                    })
  }

  if (is_null(ret)) {
    ret <- dist(data, method = method, diag = diag, upper = upper)
  }

  return(ret)
}

#----------------- Matrix utilities --------------------

coerceToDgCMatrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # DelayedArray / HDF5-backed etc.
  if (inherits(x, "DelayedMatrix")) {
    # prefer sparse realization if available
    x <- tryCatch(as(x, "dgCMatrix"),
                  error = function(e) as.matrix(x))
  }

  if (!inherits(x, "Matrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }

  return(as(x, "dgCMatrix"))
}

validateRawCounts <- function(rawData) {
  if (!inherits(rawData, "dgCMatrix")) {
    stop("Input 'raw' data must be of type dgCMatrix for this check")
  }
  x <- rawData@x
  if (anyNA(x)) {
    stop("Input 'raw' data contains NA!")
  }
  if (any(x < 0.0)) {
    stop("Input 'raw' data must contain only non negative integers.")
  }
  # integer check on nonzeros only
  if (any(x != round(x))) {
    stop("Input 'raw' data contains non integer numbers.")
  }
  return(TRUE)
}


#----------------- depreacetd functions --------------------

#' @details This is a deprecated function related to old `scCOTAN` objects. Use
#'   the more appropriate `Matrix::dspMatrix` type for similar functionality.
#'
#'   `mat2vec_rfast` converts a compacted symmetric matrix (that is an array)
#'   into a symmetric matrix.
#'
#' @param mat a square (possibly symmetric) matrix with all genes as row and
#'   column names.
#'
#' @returns `mat2vec_rfast` returns a `list` formed by two arrays:
#'   * `"genes"` with the unique gene names,
#'   * `"values"` with all the values.
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
#' @rdname COTAN_Legacy
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
        for (i in (0L:(posGene - 1L))) {
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

#' @details This is a deprecated function related to old `scCOTAN` objects. Use
#'   the more appropriate `Matrix::dspMatrix` type for similar functionality.
#'
#'   `vec2mat_rfast` converts a symmetric matrix into a compacted symmetric
#'   matrix. It will forcibly make its argument symmetric.
#'
#' @param x a `list` formed by two arrays: `genes` with the unique gene names
#'   and `values` with all the values.
#' @param genes an array with all wanted genes or the string `"all"`. When equal
#'   to `"all"` (the default), it recreates the entire matrix.
#'
#' @returns `vec2mat_rfast` returns the reconstructed symmetric matrix
#'
#' @importFrom Rfast lower_tri
#'
#' @export
#'
#' @rdname COTAN_Legacy
#'
mat2vec_rfast <- function(mat) {
  mat <- as.matrix(mat)

  if (dim(mat)[[1L]] != dim(mat)[[2L]]) {
    stop("The matrix is not square!")
  }

  v <- lower_tri(mat, diag = TRUE)
  return(list("genes" = rownames(mat), "values" = v))
}
