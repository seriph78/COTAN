# method for estimating the expected values of contingency tables
#' @export
setMethod(
  "expectedContingencyTables",
  "COTAN",
  function(objCOTAN) {
    cells <- as.matrix(objCOTAN@raw)
    # Cells matrix : formed by row data matrix changed to 0-1 matrix
    cells[cells > 0] <- 1
    cells[cells <= 0] <- 0

    mu <- estimateMu(objCOTAN)
    mu <- mu[!rownames(mu) %in% objCOTAN@hkGenes, ]
    
    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(objCOTAN@dispertion, mu[, colnames(cells)])

    message <- "Error: some Na in matrix of probability of zero UMI counts. "
    stopifnot(message = !any(is.na(probZero)))

    expectedNN <- probZero %*% t(probZero)
    expectedNY <- probZero %*% t(1 - probZero)
    expectedYN <- t(expectedNY)
    expectedYY <- (1 - probZero) %*% t(1 - probZero)

    out <- list(
      "expectedNN" = expectedNN,
      "expectedNY" = expectedNY,
      "expectedYN" = expectedYN,
      "expectedYY" = expectedYY
    )
    
    return(out)
  }
)