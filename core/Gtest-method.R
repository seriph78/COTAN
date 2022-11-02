# G-test statistics
setMethod(
  "Gtest", "COTAN",
  function(objCOTAN) {
    hk <- objCOTAN@hkGenes
    
    # OBSERVED CONTINGENCY TABLE
    table <- observedContingencyT(objCOTAN) 
    
    # exclude housekeeping genes
    table$yesYes <- table$yesYes[!rownames(table$yesYes) %in% hk, 
                                 !colnames(table$yesYes) %in% hk]
    table$yesNo  <- table$yesNo[!rownames(table$yesNo) %in% hk, 
                                !colnames(table$yesNo) %in% hk]
    table$noYes  <- table$noYes[!rownames(table$noYes) %in% hk, 
                                !colnames(table$noYes) %in% hk]
    table$noNo   <- table$noNo[!rownames(table$noNo) %in% hk, 
                               !colnames(table$noNo) %in% hk]
    
    # EXPECTED CONTINGENCY TABLE
    expectedCT <- expectedContingencyTables(object)
    
    for (i in expectedCT) {
      stopifnot("Some expected values are 0!" = !any(i == 0))
    }

    # estimation of G
    t1 <- as.matrix(table$yesYes) * log(as.matrix(table$yesYes) /
                                   as.matrix(expectedCT$expectedYY))
    t1[which(as.matrix(table$yesYes) == 0)] <- 0
    
    t2 <- as.matrix(table$yesNo) * log(as.matrix(table$yesNo) /
                                          as.matrix(expectedCT$expectedYN))
    t2[which(as.matrix(table$yesNo) == 0)] <- 0
    
    t3 <- as.matrix(table$noYes) * log(as.matrix(table$noYes) /
                                         as.matrix(expectedCT$expectedNY))
    t3[which(as.matrix(table$noYes) == 0)] <- 0
    
    t4 <- as.matrix(table$noNo) * log(as.matrix(table$noNo) /
                                         as.matrix(expectedCT$expectedNN))
    t4[which(as.matrix(table$noNo) == 0)] <- 0
    
    G <- 2 * (t1 + t2 + t3 + t4)
    
    return(G)
  }
)
