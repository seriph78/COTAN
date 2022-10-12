# method for estimating the expected values of contingency tables
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells, 
#' otherwise for the genes 
#' @export
setMethod(
  "expectedContingencyTables",
  "COTAN",
  function(objCOTAN, cells) {
    zeroOne <- as.matrix(objCOTAN@raw)
    # zeroOne matrix : formed by row data matrix changed to 0-1 matrix
    zeroOne <- sign(objCOTAN@raw)

    mu <- estimateMu(objCOTAN)
    mu <- mu[!rownames(mu) %in% objCOTAN@hkGenes, ]
    
    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(objCOTAN@dispersion, mu[, colnames(zeroOne)])

    message <- "Error: some Na in matrix of probability of zero UMI counts. "
    stopifnot(message = !any(is.na(probZero)))

    if(cells) {
      # dimension m x m (m number of cells)
      expectedNN <- t(probZero) %*% probZero
      expectedNY <- t(1 - probZero) %*% probZero
      expectedYN <- t(expectedNY)
      expectedYY <- t(1 - probZero) %*% (1 - probZero) 
    } else {
      # dimension n x n (n number of genes)
      expectedNN <- probZero %*% t(probZero)
      expectedNY <- probZero %*% t(1 - probZero)
      expectedYN <- t(expectedNY)
      expectedYY <- (1 - probZero) %*% t(1 - probZero)
    }

    out <- list(
      "expectedNN" = expectedNN,
      "expectedNY" = expectedNY,
      "expectedYN" = expectedYN,
      "expectedYY" = expectedYY
    )
    
    return(out)
  }
)