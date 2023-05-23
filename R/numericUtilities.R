#------------------- negative binomial ----------

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
#' @param disp the estimated dispersion (can be a \eqn{n}-sized vector)
#' @param mu the lambda times nu value  (can be a \eqn{n \times m} matrix)
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


#--------------------- dispersion solvers ----------------

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

    # once we have found the two bounds to the dispersion value we use bisection
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

    # once we have found the two bounds to the dispersion value we use bisection
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


#------------------------- nu solvers ---------------------

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

      diffs[runPos] <- (colSums(funProbZero(dispersion,
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


#----------------- legacy functions --------------------

#' Handle symmetric matrix <-> vector conversions
#'
#' @description Converts a symmetric matrix into a compacted symmetric matrix
#'   and vice-versa.
#'
#' @name LegacyFastSymmMatrix
NULL

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
