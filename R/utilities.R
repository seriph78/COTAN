#' setLoggingLevel
#'
#' @description Set the `COTAN` logging level
#'
#' @details
#' * 0 - Always on log messages.
#' * 1 - Major log messages.
#' * 2 - Minor log messages.
#' * 3 - All log messages.
#'
#' @seealso [logThis()]
#'
#' @param newLevel the new default logging level. It defaults to 1
#'
#' @export
#'
#' @examples
#' setLoggingLevel(3) # for debugging purposes only
#'
#' @rdname setLoggingLevel
#'
setLoggingLevel <- function(newLevel = 1) {
  message(paste0("Setting new log level to ", newLevel))
  options(COTAN.LogLevel = newLevel)
}

#' logThis
#'
#' @description Print the given message string if the current log level is
#'   greater or equal to the given log level.
#'
#' @details It uses [message()] to actually print the messages on the [strerr()]
#'   connection, so it is subject to [suppressMessages()]
#'
#' @seealso [setLoggingLevel()] for details on the levels
#'
#' @param msg the message to print
#' @param logLevel the logging level of the current message. It defaults to 2
#' @param appendLF whther to add a new-line character at the end of the message
#'
#' @returns whether the message has been printed
#'
#' @export
#'
#' @examples
#' logThis("LogLevel 0 messages will always show, ", logLevel = 0, appendLF = FALSE)
#' suppressMessages(logThis("unless all messages are suppressed", logLevel = 0))
#'
#' @rdname logThis
#'
logThis <- function(msg, logLevel = 2, appendLF = TRUE) {
  # set the logging level global variable
  if (is.null(getOption("COTAN.LogLevel"))) {
    setLoggingLevel() # to default
  }
  currentLevel <- getOption("COTAN.LogLevel")
  if (currentLevel >= logLevel) {
    message(msg, appendLF= appendLF)
    return(TRUE)
  }
  return(FALSE)
}


#' setColumnInDF
#'
#' @description Private function that append, if missing, or resets, if present,
#'   a column into a `data.frame`, whether the `data.frame` is empty or not.
#'
#' @param df the `data.frame`
#' @param colToSet the the column to add
#' @param colName the name of the new or existing column in the `data.frame`
#' @param rowNames when not empty, if the input `data.frame` has no real row
#'   names, the new row names of the resulting `data.frame`
#'
#' @returns the updated/created `data.frame`
#'
#' @importFrom rlang is_empty
#'
#' @noRd
#'
setColumnInDF <- function(df, colToSet, colName, rowNames = c()) {
  out <- df
  if(is_empty(out)) {
    out <- data.frame(colToSet)
    colnames(out) <- colName
  }
  else {
    if (colName %in% colnames(out)) {
      out[[colName]] <- colToSet
    }
    else {
      out <- cbind(out, colToSet)
      colnames(out)[ncol(out)] <- colName
    }
  }

  # default assigned rownames are numbers...
  if (!is_empty(rowNames) && is.integer(attr(out, "row.names"))) {
    rownames(out) <- rowNames
  }

  return(out)
}


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
#'   The function is continuous in \eqn{d = 0}{`disp = 0`}, increasing in \eqn{d}{`disp`}
#'   and decreasing in \eqn{\mu}{`mu`}. It returns 0 when
#'   \eqn{d = -\infty}{`disp = -Inf`} or \eqn{\mu = \infty}{`mu = Inf`}.
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
  neg <- disp <= 0
  ad <- abs(disp)
  return( ( neg) * (exp(-(1 + ad) * mu)) +
          (!neg) * (1 + ad * mu)^(-1 / ad) )
}


#' dispersionBisection
#'
#' @description Private function invoked by [estimateDispersion()] for the
#'   estimation of 'dispersion' slot of a `COTAN` object via a bisection solver
#'
#' @details The goal is to find a dispersion value that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param gene name of the relevant gene
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param muEstimator a matrix estimator of vector mu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion value
#'
#' @noRd
#'
dispersionBisection <-
  function(gene,
           zeroOneMatrix,
           muEstimator,
           threshold = 0.001,
           maxIterations = 1000) {
  sumZeros <- ncol(zeroOneMatrix) - sum(zeroOneMatrix[gene, ])
  if (sumZeros == 0){
    # cannot match exactly zero prob of zeros with finite values
    return(-Inf)
  }
  muEstimator <- muEstimator[gene, ]

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  disp1 <- 0
  diff1 <- sum(funProbZero(disp1, muEstimator)) - sumZeros
  if (abs(diff1) <= threshold) {
    return(disp1)
  }

  disp2 <- -1 * sign(diff1) # we assume error is an increasing function of disp
  iter <- 1
  repeat {
    diff2 <- sum(funProbZero(disp2, muEstimator)) - sumZeros

    if (diff2 * diff1 < 0) {
      break
    }

    if (iter >= maxIterations) {
      stop(paste("Max number of iterations reached",
                 "while finding the solution straddling intervals"))
    }
    iter <- iter + 1

    disp1 <- disp2 # disp2 is closer to producing 0

    disp2 <- 2 * disp2 # we double at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  iter <- 1
  repeat {
    disp <- (disp1 + disp2) / 2

    diff <- sum(funProbZero(disp, muEstimator)) - sumZeros

    if (abs(diff) <= threshold) {
      return(disp)
    }

    if (iter >= maxIterations) {
      stop(paste("Max number of iterations reached",
                 "while finding the solution straddling intervals"))
    }
    iter <- iter + 1

    # drop same sign diff point
    if (diff * diff2 > 0) {
      disp2 <- disp
    } else {
      disp1 <- disp
    }
  }
}


#' parallelDispersionBisection
#'
#' @description Private function invoked by [estimateDispersion()] for the
#'   estimation of the `dispersion` slot of a `COTAN` object via a parallel
#'   bisection solver
#'
#' @details The goal is to find a dispersion array that reduces to zero the
#'   difference between the number of estimated and counted zeros
#'
#' @param genes names of the relevant genes
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param muEstimator a matrix estimator of vector mu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion values
#'
#' @noRd
#'
parallelDispersionBisection <-
  function(genes,
           zeroOneMatrix,
           muEstimator,
           threshold = 0.001,
           maxIterations = 1000) {
  sumZeros <- ncol(zeroOneMatrix) - rowSums(zeroOneMatrix[genes, , drop = FALSE])
  muEstimator <- muEstimator[genes, , drop = FALSE]

  goodPos <- sumZeros != 0

  # cannot match exactly zero prob of zeros with finite values
  # so we ignore the rows with no zeros from the solver and return -Inf
  output <- rep(-Inf, length(sumZeros))

  if (sum(goodPos) == 0) {
    return(output)
  }

  sumZeros <- sumZeros[goodPos]
  muEstimator <- muEstimator[goodPos, , drop = FALSE]

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  disps1 <- rep(0, length(sumZeros))
  diffs1 <- rowSums(funProbZero(disps1, muEstimator)) - sumZeros
  if (all(abs(diffs1) <= threshold)) {
    output[goodPos] <- disp1
    return(output)
  }

  disps2 <- -1 * sign(diffs1) # we assume error is an increasing function of the dispersion
  diffs2 <- diffs1
  runPos <- rep(TRUE, length(diffs1))
  iter <- 1
  repeat {
    diffs2[runPos] <- (rowSums(funProbZero(disps2[runPos],
                                          muEstimator[runPos, , drop = FALSE])) -
                       sumZeros[runPos])

    runPos <- (diffs2 * diffs1 >= 0)

    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solution straddling intervals")
    }
    iter <- iter + 1

    disps1[runPos] <- disps2[runPos] # disps2 are closer to producing 0

    disps2[runPos] <- 2 * disps2[runPos] # we double at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  runNum <- length(diffs1)
  runPos <- rep(TRUE, runNum)
  disps <- disps1
  diffs <- diffs1
  iter <- 1
  repeat {
    disps[runPos] <- (disps1[runPos] + disps2[runPos]) / 2

    diffs[runPos] <- (rowSums(funProbZero(disps[runPos],
                                         muEstimator[runPos, , drop = FALSE])) -
                      sumZeros[runPos])

    runPos <- abs(diffs) > threshold
    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solutions")
    }
    iter <- iter + 1

    # drop same sign diff point
    pPos <- runPos & (diffs * diffs2 > 0)
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
#' @param cell name of the relevant cell
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
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
  function(cell,
           zeroOneMatrix,
           lambda,
           dispersion,
           initialGuess,
           threshold = 0.001,
           maxIterations = 1000) {
  sumZeros <- nrow(zeroOneMatrix) - sum(zeroOneMatrix[, cell])

  if (sumZeros == 0) {
    # cannot match exactly zero prob of zeros with finite values
    return(Inf)
  }

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  nu1 <- initialGuess[cell]
  diff1 <- sum(funProbZero(dispersion, nu1 * lambda)) - sumZeros
  if (abs(diff1) <= threshold) {
    return(nu1)
  }

  factor <- 2 ^ sign(diff1)
  nu2 <- nu1 * factor # we assume error is an decreasing function of nu
  iter <- 1
  repeat {
    diff2 <- sum(funProbZero(dispersion, nu2 * lambda)) - sumZeros

    if (diff2 * diff1 < 0) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solution straddling intervals")
    }
    iter <- iter + 1

    nu1 <- nu2 # nu2 is closer to producing 0

    nu2 <- nu2 * factor # we double/half at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  iter <- 1
  repeat {
    nu <- (nu1 + nu2) / 2

    diff <- sum(funProbZero(dispersion, nu * lambda)) - sumZeros

    if (abs(diff) <= threshold) {
      return(nu)
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solution straddling intervals")
    }
    iter <- iter + 1

    # drop same sign diff point
    if (diff * diff2 > 0) {
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
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param lambda the lambda vector
#' @param dispersion the dispersion vector
#' @param initialGuess the initial guess for nu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @returns the dispersion values
#'
#' @noRd
#'
parallelNuBisection <-
  function(cells,
           zeroOneMatrix,
           lambda,
           dispersion,
           initialGuess,
           threshold = 0.001,
           maxIterations = 1000) {
  sumZeros <- nrow(zeroOneMatrix) - colSums(zeroOneMatrix[, cells, drop = FALSE])
  initialGuess <- initialGuess[cells]

  stopifnot("initialGuess must hold only positive values" = all(initialGuess > 0))

  goodPos <- sumZeros != 0

  # cannot match exactly zero prob of zeros with finite values
  output <- rep(Inf, length(initialGuess))

  if (sum(goodPos) == 0) {
    return(output)
  }

  sumZeros <- sumZeros[goodPos]

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  nus1 <- initialGuess[goodPos]

  diffs1 <- colSums(funProbZero(dispersion, lambda %*% t(nus1))) - sumZeros
  if (all(abs(diffs1) <= threshold)) {
    output[goodPos] <- nus1
    return(output)
  }

  factors <- 2 ^ sign(diffs1)
  nus2 <- nus1 * factors # we assume error is an increasing function of the dispersion
  diffs2 <- diffs1
  runPos <- rep(TRUE, length(diffs1))
  iter <- 1
  repeat {
    diffs2[runPos] <- (colSums(funProbZero(dispersion, lambda %*% t(nus2[runPos]))) -
                       sumZeros[runPos])

    runPos <- (diffs2 * diffs1 >= 0)

    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solution straddling intervals")
    }
    iter <- iter + 1

    nus1[runPos] <- nus2[runPos] # nus2 are closer to producing 0

    nus2[runPos] <- nus2[runPos] * factors[runPos] # we double/half at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  runNum <- length(diffs1)
  runPos <- rep(TRUE, runNum)
  nus <- nus1
  diffs <- diffs1
  iter <- 1
  repeat {
    nus[runPos] <- (nus1[runPos] + nus2[runPos]) / 2

    diffs[runPos] <- (colSums(funProbZero(dispersion, lambda %*% t(nus[runPos]))) -
                      sumZeros[runPos])

    runPos <- (abs(diffs) > threshold)
    if (!any(runPos)) {
      break
    }

    if (iter >= maxIterations) {
      stop("Max number of iterations reached while finding the solutions")
    }
    iter <- iter + 1

    # drop same sign diff point
    pPos <- runPos & (diffs * diffs2 > 0)
    nus2[pPos] <- nus[pPos]

    nPos <- runPos & !pPos
    nus1[nPos] <- nus[nPos]
  }

  output[goodPos] <- nus
  return(output)
}


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
#' @importFrom ggplot2 unit
#'
#' @importFrom ggthemes theme_tufte
#'
#' @export
#'
#' @examples
#' theme <- plotTheme("pca")
#' theme
#'
#' @rdname plotTheme
#'
plotTheme <- function(plotKind = "common", textSize = 14) {
  myDarkBlue <- "#3C5488FF"
  ts <- textSize

  basic_theme <- theme(
    axis.text.x  = element_text(size = ts, angle = 0,  hjust = .5, vjust = .5, face = "plain", colour = myDarkBlue),
    axis.text.y  = element_text(size = ts, angle = 0,  hjust = 0,  vjust = .5, face = "plain", colour = myDarkBlue),
    axis.title.x = element_text(size = ts, angle = 0,  hjust = .5, vjust = 0,  face = "plain", colour = myDarkBlue),
    axis.title.y = element_text(size = ts, angle = 90, hjust = .5, vjust = .5, face = "plain", colour = myDarkBlue) )

  if (plotKind == "common") {
    return(basic_theme)
  }

  if (plotKind == "pca") {
    return( basic_theme +
            theme(legend.title = element_blank(),
                  legend.text = element_text(size = 12, face = "italic", color = myDarkBlue),
                  legend.position = "bottom") )
  }

  if (plotKind == "genes") {
    return( basic_theme +
            theme(plot.title = element_text(size = 20, hjust = 0.02, vjust = -10,
                                            face = "italic", color = myDarkBlue),
                  plot.subtitle = element_text(vjust = -15, hjust = 0.01,
                                               color = "darkred")) )
  }

  if (plotKind == "UDE") {
    return( basic_theme +
            theme(plot.title   = element_text(size = 20, color = myDarkBlue),
                  legend.title = element_text(size = 14, color = myDarkBlue, face = "italic"),
                  legend.text  = element_text(size = 11, color = myDarkBlue),
                  legend.key.width = unit(2, "mm"),
                  legend.position  ="right") )
  }

  if (plotKind == "heatmap") {
    return( basic_theme +
            theme( axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   panel.spacing = unit(0, "lines"),
                   strip.background = element_rect(fill = "#8491B44C"),
                   strip.text.y = element_text(size = ts, colour = myDarkBlue),
                   strip.text.x = element_text(size = ts, angle = 90, colour = myDarkBlue),
                   legend.text = element_text(color = myDarkBlue, face = "italic"),
                   legend.position = "bottom",
                   legend.title = element_blank(),
                   legend.key.height = unit(2, "mm") ) )
  }

  if (plotKind == "GDI") {
    return( basic_theme +
            theme( legend.title = element_blank(),
                   plot.title = element_text(size = ts + 2, face = "bold.italic", color = myDarkBlue),
                   legend.text = element_text(color = myDarkBlue, face = "italic"),
                   legend.position = "bottom") )
  }

  if (plotKind == "size-plot") {
    return( ggthemes::theme_tufte() +
            theme(legend.position = "none") )
                  # axis.text.x  = element_blank(),
                  # axis.ticks.x = element_blank()) )
  }

  warning(paste("plotTheme: no match found in listed themes for:", plotKind))
  return(basic_theme)
}


#----------------- legacy functions


#' vec2mat_rfast
#'
#' @description Converts a compacted symmetric matrix back into a proper matrix
#'
#' @details This is a lagacy function related to old `scCOTAN` objects. Use the
#'   more appropriate `Matrix::dspMatrix` type for similar functionality
#'
#' @seealso [calculateCoex()] or [expectedContingencyTables()] for actual
#'   examples
#'
#' @param x a list formed by two arrays: `genes` with the unique gene names and
#'   `values` with all the values.
#' @param genes a vector with all wanted genes or the string `"all"`. When equal
#'   to `"all"` (the default), it recreates the entire matrix.
#'
#' @returns a matrix
#'
#' @importFrom Rfast lower_tri.assign
#' @importFrom Rfast upper_tri.assign
#' @importFrom Rfast transpose
#'
#' @export
#'
#' @examples
#' v <- list("genes" = paste0("gene.", c(1:10)), "values" = c(1:55))
#' genes <- c("gene.3", "gene.4", "gene.7")
#' M <- vec2mat_rfast(v)
#' m <- vec2mat_rfast(v, genes)
#'
#' @rdname vec2mat_rfast
#'
setMethod(
  "vec2mat_rfast",
  "list",
  function(x, genes = "all") {
    numGenes <- length(x$genes)
    if (genes[1] == "all") {
      m <- matrix(0, numGenes, numGenes)

      m <- Rfast::lower_tri.assign(m, diag = TRUE, v = x$values)
      m <- Rfast::upper_tri.assign(m, v = Rfast::upper_tri(Rfast::transpose(m)))
      rownames(m) <- x$genes
      colnames(m) <- x$genes
    } else {
      m <- matrix(0, nrow = numGenes, ncol = length(genes))
      rownames(m) <- x$genes
      colnames(m) <- genes

      for (pos.gene in match(genes, x$genes)) {
        temp.array <- x$values[pos.gene]
        p = 1
        l <- numGenes
        s <- pos.gene
        while (p < (pos.gene)) {
          l <- l - 1
          s <- s + l
          temp.array <- c(temp.array, x$values[s])
          p <- p + 1
        }
        if(pos.gene < numGenes) {
          #linear part
          start.reading.position <- 1
          i <- 1
          while (i < (pos.gene)) {
            start.reading.position <- start.reading.position + (numGenes - (i - 1))
            i <- i + 1
          }
          start.reading.position = start.reading.position+1

          end.reading.position <- 0
          for (i in c(0:(pos.gene-1))) {
            end.reading.position <- end.reading.position + (numGenes-i)
          }

          temp.array <- c(temp.array, x$values[start.reading.position:end.reading.position])
        }
        m[,x$genes[pos.gene]] <- temp.array
      }
    }
    return(m)
  }
)


#' mat2vec_rfast
#'
#' @description Converts a square matrix into a compact form, by forcibly make
#'   it symmetric
#'
#' @details This is a lagacy function related to old `scCOTAN` objects. Use the
#'   more appropriate `Matrix::dspMatrix` type for similar functionality
#'
#' @seealso [calculateCoex()] or [expectedContingencyTables()] for actual
#'   examples
#'
#' @param mat a possibly symmetric matrix with all genes as row and column names
#'
#' @returns a `list` formed by two arrays: `genes` with the unique gene names
#'   and `values` with all the values.
#'
#' @importFrom Rfast lower_tri.assign
#' @importFrom Rfast upper_tri.assign
#' @importFrom Rfast transpose
#'
#' @export
#'
#' @examples
#' mat <- matrix(0, nrow = 10, ncol = 10)
#' mat <- Rfast::lower_tri.assign(mat, c(1:55), diag = TRUE)
#' mat <- Rfast::upper_tri.assign(mat, v = Rfast::upper_tri(Rfast::transpose(mat)))
#' v <- mat2vec_rfast(mat)
#'
#' @rdname mat2vec_rfast
#'
setMethod(
  "mat2vec_rfast",
  "matrix",
  function(mat) {
    if (!dim(mat)[1] == dim(mat)[2]) {
      stop("The matrix is not square!")
    }

    v <- Rfast::lower_tri(mat, diag = TRUE)
    return(list("genes" = rownames(mat), "values" = v))
  }
)

