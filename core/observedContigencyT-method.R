setGeneric("observedContigencyT", function(objCOTAN) standardGeneric("observedContigencyT"))
setMethod(
  "observedContigencyT", "COTAN",
  function(objCOTAN) {
    zeroOne <- sign(objCOTAN@raw)
    
    sum <- rowSums(zeroOne)
    sum <- as.matrix(somma)
    yesAny <- do.call("cbind", 
                      replicate(length(rownames(sum)), 
                      sum, 
                      simplify = FALSE))
    colnames(yesAny) <- rownames(yesAny)
    
    yesYes <- as.matrix(observedContingencyYY(objCOTAN))
    yesNo  <- yesAny - yesYes
    noYes  <- t(yesAny) - yesYes
    noNo   <- length(colnames(zeroOne)) - (yesYes + noYes + yesNo)
    
    out <- list("yesYes" = yesYes,
                "noYes" = noYes, 
                "yesNo" = yesNo, 
                "noNo" = noNo)
    return(out)
  }
)
