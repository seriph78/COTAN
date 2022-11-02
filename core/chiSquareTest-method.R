# Pearson's chi-squared test statistics
setMethod(
  "chiSquaredTest", "COTAN",
  function(objCOTAN) {
    test <- (objCOTAN@coex$values)^2 * objCOTAN@nCells
    test <- list("genes" = objCOTAN@coex$genes, "values" = test)
    return(test)
  }
)