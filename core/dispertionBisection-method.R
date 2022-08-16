#' dispertionBisection
#'
#' private method invoked by 'estimateDispertionBisection' for the estimation 
#' of 'dispertion' field of a COTAN object with bisection
#' @param genes name of the genes
#' @param zeroOneMatrix 
#' @param muEstimator estimator of vector mu
#' @param threshold
#' @return r, data.frame(a, u)
setMethod(
  "dispertionBisection",
  "character",
  function(genes, zeroOneMatrix, muEstimator, threshold) {
    sumZeros <- sum(zeroOneMatrix[genes, ] == 0)
    muEstimatorgenes <- muEstimator[genes, ]

    a1 <- 0
    u1 <- diffZeros(a1, sumZeros, muEstimatorgenes)
    a2 <- a1
    u2 <- u1
    if (u1 > 0) {
      a1 <- a1 - 1
      u1 <- diffZeros(a1, sumZeros, muEstimatorgenes)
      while (u1 > 0) {
        a2 <- a1
        u2 <- u1
        a1 <- 2 * a1
        u1 <- diffZeros(a1, sumZeros, muEstimatorgenes)
      }
    } else {
      a2 <- 1
      u2 <- diffZeros(a2, sumZeros, muEstimatorgenes)
      while (u2 < 0) {
        a1 <- a2
        u1 <- u2
        a2 <- 2 * a2
        u2 <- diffZeros(a2, sumZeros, muEstimatorgenes)
      }
    }
    a <- (a1 + a2) / 2
    u <- diffZeros(a, sumZeros, muEstimatorgenes)
    while (abs(u) > threshold) {
      if (u > 0) {
        a2 <- a
        u2 <- u
      } else {
        a1 <- a
        u1 <- u
      }
      a <- (a1 + a2) / 2
      u <- diffZeros(a, sumZeros, muEstimatorgenes)
    }
    
    r <- data.frame(a, u)
    rownames(r) <- genes
    return(r)
  }
)
