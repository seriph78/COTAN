funProbZero <- function(disp, mu) {
  (disp <= 0) * (exp(-(1 + abs(disp)) * mu)) +
    (disp > 0) * (1 + abs(disp) * mu)^(-1 / abs(disp))
}

#' dispertionBisection
#'
#' private function invoked by 'estimateDispertionBisection' for the estimation
#' of 'dispertion' field of a COTAN object with bisection
#'
#' the goal is to find dispertion value that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param genes name of the genes
#' @param zeroOneMatrix
#' @param muEstimator estimator of vector mu
#' @param threshold
#' @return r, data.frame(a, u)
dispertionBisection <- function(genes,
                                zeroOneMatrix,
                                muEstimator,
                                threshold = 0.001) {
  sumZeros <- sum(zeroOneMatrix[genes, ] == 0)
  muEstimator <- muEstimator[genes, ]

  # we look for two dispertion values where the first leads to a
  # diffZeros negative and the second positive
  dispertion1 <- 0
  dispertion2 <- 0
  negativeDiff <- sum(funProbZero(dispertion1, muEstimator)) - sumZeros
  positiveDiff <- negativeDiff

  if (negativeDiff > 0) {
    dispertion1 <- -1
    negativeDiff <- sum(funProbZero(dispertion1, muEstimator)) - sumZeros
    while (negativeDiff > 0) {
      dispertion2 <- dispertion1 # dispertion1 is closer to producing 0
      positiveDiff <- negativeDiff
      dispertion1 <- 2 * dispertion1 # we double at each step
      negativeDiff <- sum(funProbZero(dispertion1, muEstimator)) - sumZeros
    }
  } else {
    dispertion2 <- 1
    positiveDiff <- sum(funProbZero(dispertion2, muEstimator)) - sumZeros
    while (positiveDiff < 0) {
      dispertion1 <- dispertion2 # dispertion2 is closer to producing 0
      negativeDiff <- positiveDiff
      dispertion2 <- 2 * dispertion2 # we double at each step
      positiveDiff <- sum(funProbZero(dispertion2, muEstimator)) - sumZeros
    }
  }

  # once we have found the two dispersion values, we use bisection
  dispertion <- (dispertion1 + dispertion2) / 2
  diff <- sum(funProbZero(dispertion, muEstimator)) - sumZeros
  while (abs(diff) > threshold) {
    if (diff > 0) {
      dispertion2 <- dispertion
      positiveDiff <- diff
    } else {
      dispertion1 <- dispertion
      negativeDiff <- diff
    }
    dispertion <- (dispertion1 + dispertion2) / 2
    diff <- sum(funProbZero(dispertion, muEstimator)) - sumZeros
  }

  r <- data.frame(dispertion, diff)
  rownames(r) <- genes
  return(r)
}
