funProbZero <- function(disp, mu) {
  (disp <= 0) * (exp(-(1 + abs(disp)) * mu)) +
    (disp > 0) * (1 + abs(disp) * mu)^(-1 / abs(disp))
}

#' dispersionBisection
#'
#' private function invoked by 'estimateDispersionBisection' for the estimation
#' of 'dispersion' field of a COTAN object with bisection
#'
#' the goal is to find dispersion value that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param genes name of the genes
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param muEstimator estimator of vector mu
#' @param threshold minimal solution precision
#' @return r, data.frame(a, u)
dispersionBisection <- function(genes,
                                zeroOneMatrix,
                                muEstimator,
                                threshold = 0.001) {
  sumZeros <- sum(zeroOneMatrix[genes, ] == 0)
  muEstimator <- muEstimator[genes, ]

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  dispersion1 <- 0
  dispersion2 <- 0
  negativeDiff <- sum(funProbZero(dispersion1, muEstimator)) - sumZeros
  positiveDiff <- negativeDiff

  if (negativeDiff > 0) {
    dispersion1 <- -1
    negativeDiff <- sum(funProbZero(dispersion1, muEstimator)) - sumZeros
    while (negativeDiff > 0) {
      dispersion2 <- dispersion1 # dispersion1 is closer to producing 0
      positiveDiff <- negativeDiff
      dispersion1 <- 2 * dispersion1 # we double at each step
      negativeDiff <- sum(funProbZero(dispersion1, muEstimator)) - sumZeros
    }
  } else {
    dispersion2 <- 1
    positiveDiff <- sum(funProbZero(dispersion2, muEstimator)) - sumZeros
    while (positiveDiff < 0) {
      dispersion1 <- dispersion2 # dispersion2 is closer to producing 0
      negativeDiff <- positiveDiff
      dispersion2 <- 2 * dispersion2 # we double at each step
      positiveDiff <- sum(funProbZero(dispersion2, muEstimator)) - sumZeros
    }
  }

  # once we have found the two dispersion values, we use bisection
  dispersion <- (dispersion1 + dispersion2) / 2
  diff <- sum(funProbZero(dispersion, muEstimator)) - sumZeros
  while (abs(diff) > threshold) {
    if (diff > 0) {
      dispersion2 <- dispersion
      positiveDiff <- diff
    } else {
      dispersion1 <- dispersion
      negativeDiff <- diff
    }
    dispersion <- (dispersion1 + dispersion2) / 2
    diff <- sum(funProbZero(dispersion, muEstimator)) - sumZeros
  }

  r <- data.frame(dispersion, diff)
  rownames(r) <- genes
  return(r)
}
