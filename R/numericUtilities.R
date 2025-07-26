#'
#' @title Numeric Utilities
#'
#' @description A set of function helper related to the statistical model
#'   underlying the `COTAN` package
#'
#' @name NumericUtilities
NULL

#------------------- UMI count models ----------

#'
#' @details `funProbZeroNegBin` is a private function that gives the probability
#'   that a sample gene's reads are zero, given the `dispersion` and `mu`
#'   parameters of a *Negative Binomial* model
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
#'   It returns \eqn{1.0} when \eqn{\mu = 0}{`mu = 0`}.
#'
#' @param dispersion the estimated `dispersion` (a \eqn{n}-sized vector)
#' @param mu the `lambda` times `nu` values (a \eqn{n \times m} matrix)
#'
#' @returns `funProbZeroNegBin()` returns the probability `matrix` that a *read
#'   count* is identically zero
#'
#' @rdname NumericUtilities
#'
funProbZeroNegBin <- function(dispersion, mu) {
  neg <- (dispersion <= 0.0)
  ad <- abs(dispersion)
  return(( neg) * (exp(-(1.0 + ad) * mu)) +
         (!neg) * (1.0 + ad * mu)^(-1.0 / ad))
}


#'
#' @details `funProbZeroMixPoi` is a private function that gives the probability
#'   that a sample gene's reads are zero, given the `pi` and `mu` parameters of
#'   a *Mixture with Poisson* model
#'
#' @details Using \eqn{\pi}{`pi`} for `pi` and \eqn{\mu}{`mu`} for `mu`, it
#'   returns: \eqn{\pi + (1 - \pi) e^{-\mu}}
#'   It returns \eqn{1.0} when \eqn{\mu = 0}{`mu = 0`}.
#' @param pi the estimated *probability* that a cell does not express a given
#'   gene (a \eqn{n}-sized vector)
#' @param mu the `lambda` times `nu` values (a \eqn{n \times m} matrix)
#'
#' @returns `funProbZeroMixPoi()` returns the probability `matrix` that a *read
#'   count* is identically zero
#'
#' @rdname NumericUtilities
#'
funProbZeroMixPoi <- function(pi, mu) {
  return(pi + (1 - pi) * exp(-mu))
}


#---------------- dispersion solvers negative binomial model ------------

#' @details `dispersionBisection()` is a private function for the estimation of
#'   `dispersion` slot of a `COTAN` object via a bisection solver
#'
#' @details The goal of `dispersionBisection()` is to find a `dispersion` value
#'   that reduces to zero the difference between the number of estimated and
#'   counted zeros
#'
#' @param sumZeros the number of cells that didn't express the gene
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns `dispersionBisection()` returns the `dispersion` value
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
    diff1 <- sum(funProbZeroNegBin(disp1, mu)) - sumZeros
    if (abs(diff1) <= threshold) {
      return(disp1)
    }

    # we assume error is an increasing function of disp
    disp2 <- -1.0 * sign(diff1)
    iter <- 1L
    repeat {
      diff2 <- sum(funProbZeroNegBin(disp2, mu)) - sumZeros

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

    # once we have found the two bounds to the dispersion value we use bisection
    iter <- 1L
    repeat {
      disp <- (disp1 + disp2) / 2.0

      diff <- sum(funProbZeroNegBin(disp, mu)) - sumZeros

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


#' @details `parallelDispersionBisection()` is a private function invoked by
#'   [estimateDispersionBisection()] for the estimation of the `dispersion` slot
#'   of a `COTAN` object via a parallel bisection solver
#'
#' @details The goal of `parallelDispersionBisection()` is to find a
#'   `dispersion` `array` that reduces to zero the difference between the number
#'   of estimated and counted zeros
#'
#' @param genes names of the relevant genes
#' @param sumZeros the number of cells that didn't express the relevant gene (a
#'   \eqn{n}-sized vector)
#' @param lambda the estimated `lambda` (a \eqn{n}-sized vector)
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns `parallelDispersionBisection()` returns the dispersion values
#'
#' @importFrom Rfast rowsums
#'
#' @importFrom rlang rep_along
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
    output <- rep_along(sumZeros, -Inf)

    # in case of zero lambda dispersion is irrelevant. We return 1.0
    goodPos <- sumZeros != length(nu)
    output[!goodPos] <- 1.0

    goodPos <- goodPos & sumZeros != 0L

    if (sum(goodPos) == 0L) {
      return(output)
    }

    sumZeros <- sumZeros[goodPos]
    lambda <- lambda[goodPos]

    mu <- lambda %o% nu

    # we look for two dispersion values where the first leads to a
    # diffZeros negative and the second positive
    disps1 <- rep_along(sumZeros, 0.0)
    diffs1 <- rowsums(funProbZeroNegBin(disps1, mu)) - sumZeros
    if (all(abs(diffs1) <= threshold)) {
      output[goodPos] <- disps1
      return(output)
    }

    # we assume error is an increasing function of the dispersion
    disps2 <- -1.0 * sign(diffs1)
    diffs2 <- diffs1
    runPos <- rep_along(diffs1, TRUE)
    iter <- 1L
    repeat {
      diffs2[runPos] <-
        (rowsums(funProbZeroNegBin(disps2[runPos],
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

    # once we have found the two bounds to the dispersion value we use bisection
    runPos <- rep_along(diffs1, TRUE)
    disps <- disps1
    diffs <- diffs1
    iter <- 1L
    repeat {
      disps[runPos] <- (disps1[runPos] + disps2[runPos]) / 2.0

      diffs[runPos] <-
        (rowsums(funProbZeroNegBin(disps[runPos],
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

    rm(mu)
    gc()

    output[goodPos] <- disps
    return(output)
  }



#--------------------- lambda solvers mixture model ----------------

#' @details `lambdaNewton()` is a private function for the estimation of
#'   `lambda` slot of a `COTAN` object via a Newton-Raphson solver
#'
#' @details The goal of `lambdaNewton()` is to find a `lambda` value by which
#'   the equation \eqn{\frac{1 - e^{-\lambda}}{\lambda} = \frac{1 - z}{m}},
#'   coming from the mixture model, so to match contemporaneously both average
#'   and number of zero in the model. The formula is actually \eqn{\nu} weighted
#'
#' @param avgNumNonZeros the average number of cells that express the gene
#' @param avgCounts the average expression of the gene
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param zeroPi whether this genes is a special case where \eqn{p \equiv 0} and
#'   only the zero probability marginals are calibrated
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns `lambdaNewton()` returns the lambda value
#'
#' @importFrom assertthat assert_that
#'
#' @rdname NumericUtilities
#'
lambdaNewton <-
  function(avgNumNonZeros,
           avgCounts,
           nu,
           zeroPi,
           threshold = 0.0001,
           maxIterations = 100L) {
    assert_that(avgNumNonZeros != 1.0 || zeroPi,
                avgNumNonZeros != avgCounts || zeroPi,
                msg = "Inconsistent `zeroPi` given")
    # handle boundary cases
    if (zeroPi) {
      if (avgNumNonZeros == 1.0) {
        return(+Inf)
      }
      if (avgNumNonZeros == 0.0) {
        assert_that(avgCounts == 0.0,
                    msg = "Inconsistent `avgCounts` given")
        return(0.0)
      }
    }

    #initial guess
    ratio <- avgNumNonZeros / avgCounts
    lambda <- 1.0 / ratio

    # Newton-Raphson loop
    iter <- 1L
    repeat {
      mu <- lambda * nu
      expMu <- exp(-mu)

      avgExpMu <- mean(expMu)

      if (zeroPi) {
        tmp1 <- (1.0 - avgExpMu - avgNumNonZeros)
      } else {
        tmp1 <- (1.0 - avgExpMu - lambda * ratio)
      }

      relDiff <- (tmp1 / lambda) / lambda

      if (abs(relDiff) <= threshold) {
        return(lambda)
      }

      if (iter >= maxIterations) {
        stop("Max number of iterations reached ",
             "while finding the solution straddling intervals")
      }

      avgMuExpMu <- mean(mu * expMu)

      #update lambda
      if (zeroPi) {
        fact <- 1.0 - tmp1 / avgMuExpMu
      } else {
        fact <- 1.0 + tmp1 / (1.0 - avgExpMu - avgMuExpMu)
      }

      # print(paste0(iter, ": factor ", fact, ", lambda ", lambda))

      iter <- iter + 1L
      lambda <- lambda * fact
    }

    return(lambda)
  }



#'
#' @returns `lambdaNewton()` returns the lambda value
#'

#' @details `parallelLambdaNewton()` is a private function invoked by
#'   [estimateLambdaPiNewton()] for the estimation of the `lambda` slot of a
#'   `COTAN` object via a Newton-Raphson solver
#'
#' @details The goal of `parallelLambdaNewton()` is to find a `lambda array`
#'   that solves the equation \eqn{\frac{1 - e^{-\lambda}}{\lambda} = \frac{1 -
#'   z}{m}}, coming from the mixture model, so to match contemporaneously both
#'   average and number of zero in the model. The formula is actually \eqn{\nu}
#'   weighted
#'
#' @param genes names of the relevant genes (a \eqn{n}-sized vector)
#' @param avgNumNonZeros the average number of cells that express the gene  (a
#'   \eqn{n}-sized vector)
#' @param avgCounts the average expression of the gene  (a \eqn{n}-sized vector)
#' @param nu the estimated `nu` (a \eqn{m}-sized vector)
#' @param zeroPi whether this genes is a special case where \eqn{p \equiv 0} and
#'   only the zero probability marginals are calibrated (a \eqn{n}-sized vector)
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns `parallelLambdaNewton()` returns the `lambda` values
#'
#' @importFrom Rfast rowmeans
#'
#' @importFrom rlang rep_along
#'
#' @importFrom assertthat assert_that
#'
#' @rdname NumericUtilities
#'
parallelLambdaNewton <-
  function(genes,
           avgNumNonZeros,
           avgCounts,
           nu,
           zeroPi,
           threshold = 0.0001,
           maxIterations = 100L) {
    avgNumNonZeros <- avgNumNonZeros[genes]
    avgCounts <- avgCounts[genes]
    zeroPi <- zeroPi[genes]

    assert_that(all(avgNumNonZeros != 1.0 | zeroPi),
                all(avgNumNonZeros != avgCounts | zeroPi),
                msg = "Inconsistent `zeroPi` given")

    # handle boundary cases
    output <- rep_along(avgCounts, NaN)

    noZeros <- avgNumNonZeros == 1.0
    output[noZeros] <- +Inf

    assert_that(identical(avgNumNonZeros == 0.0, avgCounts == 0.0),
                msg = "Inconsistent `avgCounts` given")

    onlyZeros <- avgNumNonZeros == 0.0
    output[onlyZeros] <- 0.0

    goodPos <- !(noZeros | onlyZeros)

    if (all(!goodPos)) {
      return(output)
    }

    #initial guess
    ratio <- rep_along(avgCounts, 1.0)
    ratio[goodPos] <- avgNumNonZeros[goodPos] / avgCounts[goodPos]
    lambda <- 1.0 / ratio

    # Newton-Raphson loop
    runPos <- goodPos
    relDiff <- rep_along(runPos, 0.0)
    tmp2 <- rep_along(runPos, NaN)
    iter <- 1L
    repeat {
      lambda1 <- lambda[runPos]
      zeroPi1 <- zeroPi[runPos]

      mu <- lambda1 %o% nu
      expMu <- exp(-mu)
      avgExpMu <- rowmeans(expMu)

      mu <- mu * expMu
      avgMuExpMu <- rowmeans(mu)

      tmp1 <- rep_along(lambda1, NaN)
      tmp1[ zeroPi1] <- (1.0 - avgExpMu[ zeroPi1] -
                           avgNumNonZeros[runPos][zeroPi1])
      tmp1[!zeroPi1] <- (1.0 - avgExpMu[!zeroPi1] -
                           lambda1[!zeroPi1] * ratio[runPos][!zeroPi1])

      relDiff[runPos] <- (tmp1 / lambda1) / lambda1

      runPos1 <- abs(relDiff) > threshold
      if (all(!runPos1)) {
        break
      }

      if (iter >= maxIterations) {
        stop("Max number of iterations reached while finding the solutions")
      }
      iter <- iter + 1L

      # update lambda
      tmp2[runPos][ zeroPi1] <-
        lambda1[ zeroPi1] * (1.0 - tmp1[ zeroPi1] / avgMuExpMu[zeroPi1])
      tmp2[runPos][!zeroPi1] <-
        lambda1[!zeroPi1] * (1.0 + tmp1[!zeroPi1] / (1.0 - avgExpMu[!zeroPi1] -
                                                      avgMuExpMu[!zeroPi1]))

      # update only the lambda not already OK
      runPos <- runPos1
      lambda[runPos] <- tmp2[runPos]
    }

    output[goodPos] <- lambda[goodPos]
    return(output)
  }


#------------------ nu solvers negative binomial model -------------------

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
    diff1 <- sum(funProbZeroNegBin(dispersion, nu1 * lambda)) - sumZeros
    if (abs(diff1) <= threshold) {
      return(nu1)
    }

    factor <- 2.0 ^ sign(diff1)
    nu2 <- nu1 * factor # we assume error is an decreasing function of nu
    iter <- 1L
    repeat {
      diff2 <- sum(funProbZeroNegBin(dispersion, nu2 * lambda)) - sumZeros

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

      diff <- sum(funProbZeroNegBin(dispersion, nu * lambda)) - sumZeros

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
#' @importFrom rlang rep_along
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
    output <- rep_along(initialGuess, Inf)

    if (sum(goodPos) == 0L) {
      return(output)
    }

    sumZeros <- sumZeros[goodPos]

    # we look for two dispersion values where the first leads to a
    # diffZeros negative and the second positive
    nus1 <- initialGuess[goodPos]

    diffs1 <- colsums(funProbZeroNegBin(dispersion, lambda %o% nus1)) - sumZeros
    if (all(abs(diffs1) <= threshold)) {
      output[goodPos] <- nus1
      return(output)
    }

    # we assume error is an increasing function of the dispersion
    factors <- 2.0 ^ sign(diffs1)
    nus2 <- nus1 * factors
    diffs2 <- diffs1
    runPos <- rep_along(diffs1, TRUE)
    iter <- 1L
    repeat {
      diffs2[runPos] <-
        (colsums(funProbZeroNegBin(dispersion,
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
    runPos <- rep_along(diffs1, TRUE)
    nus <- nus1
    diffs <- diffs1
    iter <- 1L
    repeat {
      nus[runPos] <- (nus1[runPos] + nus2[runPos]) / 2.0

      diffs[runPos] <-
        (colsums(funProbZeroNegBin(dispersion,
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

  if (!dim(mat)[[1L]] == dim(mat)[[2L]]) {
    stop("The matrix is not square!")
  }

  v <- lower_tri(mat, diag = TRUE)
  return(list("genes" = rownames(mat), "values" = v))
}
