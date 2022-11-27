#' setColumnInDF
#'
#' private function that append if missing or resets if present a column
#' into a data.frame, whether the data.frame is empty or not.
#'
#' @param df the data.frame
#' @param colToSet the the column to add
#' @param colName the data.frame
#' @param rowNames when not empty the new row names of the result
#'
#' @return the updated/created data.frame
#'
#' @importFrom rlang is_empty
#'
#' @rdname setColumnInDF
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

  if (!is_empty(rowNames) && is_empty(rownames(out))) {
    rownames(out) <- rowNames
  }

  return(out)
}


#' funProbZero
#'
#' private function that gives the probability of a sample gene count
#' being zero given the given the dispersion and mu
#'
#' @param disp the estimated dispersion
#' @param mu the lambda times nu value
#'
#' @return the probability that a count is identically zero
#'
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
#' private function invoked by 'estimateDispersion' for the estimation
#' of 'dispersion' slot of a COTAN object via a bisection solver
#'
#' the goal is to find a dispersion value that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param gene name of the relevant gene
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param muEstimator a matrix estimator of vector mu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @return the dispersion value
#'
#' @rdname dispersionBisection
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
      stop("Max number of iterations reached while finding the solution straddling intervals")
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
      stop("Max number of iterations reached while finding the solution straddling intervals")
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
#' private function invoked by 'estimateDispersion' for the estimation
#' of 'dispersion' slot of a COTAN object via a parallel bisection solver
#'
#' the goal is to find a dispersion array that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param genes names of the relevant genes
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param muEstimator a matrix estimator of vector mu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @return the dispersion values
#'
#' @rdname parallelDispersionBisection
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
#' private function invoked by 'estimateNuBisection' for the estimation
#' of 'nu' slot of a COTAN object via a bisection solver
#'
#' the goal is to find a nu value that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param cell name of the relevant cell
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param lambda the lambda vector
#' @param dispersion the dispersion vector
#' @param initialGuess the initial guess for nu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @return the nu value
#'
#' @rdname nuBisection
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
    return(initialGuess)
  }

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  nu1 <- initialGuess
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
#' private function invoked by 'estimateNuBisection' for the estimation
#' of 'nu' slot of a COTAN object via a parallel bisection solver
#'
#' the goal is to find a nu array that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param cells names of the relevant cells
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param lambda the lambda vector
#' @param dispersion the dispersion vector
#' @param initialGuess the initial guess for nu
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#'
#' @return the dispersion values
#'
#' @rdname parallelNuBisection
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
  output <- rep(initialGuess, length(sumZeros))

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
#' This function returns the appropriate theme for the selected plot kind
#'
#' @param plotKind a string indicating the plot kind; supported kinds are:
#' "common", "pca", "genes", "UDE", "heatmap", "GDI", "size-plot"
#' @param textSize axes and strip text size (default=14)
#' @return a ggplot2 theme object
#
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 unit
#'
#' @importFrom ggthemes theme_tufte
#'
#' @export
#' @rdname plotTheme
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

