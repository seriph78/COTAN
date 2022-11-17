#' observedContingencyYY
#'
#' calculate observed yes/yes field of contingency table
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells,
#' otherwise for the genes
#' @return the YesYes observed contengency table
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix symmetricMatrix
#' @export
#'
#' @rdname observedContingencyYY
setMethod(
  "observedContingencyYY",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE) {
    zeroOne <- getZeroOneProj(objCOTAN)

    if(isTRUE(actOnCells)){
      # for cells
      YY <- t(zeroOne) %*% zeroOne
    } else{
      # for genes
      YY <- zeroOne %*% t(zeroOne)
    }

    #return(as(YY, "symmetricMatrix"))
    return(YY)
  }
)


#' expectedContingencyTables
#'
#' method for estimating the expected values of contingency tables
#' @param objCOTAN A COTAN object
#' @param actOnCells Boolean, if TRUE, the function works for the cells,
#' otherwise for the genes
#'
#' @return a list with the expected contengency tables
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix symmetricMatrix
#'
#' @export
#'
#' @rdname expectedContingencyTables
setMethod(
  "expectedContingencyTables",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE) {
    mu <- estimateMu(objCOTAN)[flagNotHousekeepingGenes(objCOTAN), ]

    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(getDispersion(objCOTAN), mu)

    errMgs <- "Error: some NA in matrix of probability of zero UMI counts. "
    stopifnot(errMgs = !anyNA(probZero))

    cat("calculating NN..")
    if(actOnCells) {
      # dimension m x m (m number of cells)
      expectedNN <- t(probZero) %*% probZero
      cat("NY..")
      expectedNY <- t(1 - probZero) %*% probZero
      cat("YN..")
      expectedYN <- t(expectedNY)
      cat("YY..")
      expectedYY <- t(1 - probZero) %*% (1 - probZero)
    } else {
      # dimension n x n (n number of genes)
      expectedNN <- probZero %*% t(probZero)
      cat("NY..")
      expectedNY <- probZero %*% t(1 - probZero)
      cat("YN..")
      expectedYN <- t(expectedNY)
      cat("YY..")
      expectedYY <- (1 - probZero) %*% t(1 - probZero)
    }
    cat(" done\n")

    out <- list( "expectedNN" = as.matrix(expectedNN),
                 "expectedNY" = as.matrix(expectedNY),
                 "expectedYN" = as.matrix(expectedYN),
                 "expectedYY" = as.matrix(expectedYY) )

    return(out)
  }
)


#' calculateCoex
#'
#' This function estimates and stores the coex matrix in the 'cellCoex' or
#' coex field depending on given actOnCells flag
#'
#' @param objCOTAN A COTAN object
#' @param actOnCells Boolean, if TRUE, the function works for the cells,
#' otherwise for the genes
#'
#' @return the updated COTAN object
#'
#' @export
#'
#' @rdname calculateCoex
setMethod(
  "calculateCoex",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE) {
    if(actOnCells) { kind <- "cells" } else { kind <- "genes" }
    print(paste("Calculate", kind, "coex: START"))

    if (is_empty(getHousekeepingGenes(objCOTAN))) {
      objCOTAN <- findHousekeepingGenes(objCOTAN)
    }

    print("Retrieving observed yes/yes contingency table")

    observedYY <- observedContingencyYY(objCOTAN, actOnCells)

    if(!actOnCells){
      # remove the housekeeping genes
      noHKFlags <- flagNotHousekeepingGenes(objCOTAN)
      observedYY <- observedYY[noHKFlags, noHKFlags]
    }

    observedYY <- mat2vec_rfast(as.matrix(observedYY))

    gc()

    print("Retrieving expected contingency table")

    # four estimators: expectedNN, expectedNY, expectedYN, expectedYY
    list[expectedNN, expectedNY, expectedYN, expectedYY] <-
      expectedContingencyTables(objCOTAN, actOnCells)

    expectedNN <- mat2vec_rfast(expectedNN)
    expectedNY <- mat2vec_rfast(expectedNY)
    expectedYN <- mat2vec_rfast(expectedYN)
    expectedYY <- mat2vec_rfast(expectedYY)

    gc()

    print("Estimating coex")

    # sum for division
    sumForDiv <- ( 1 / pmax(1, expectedYY$values) +
                   1 / pmax(1, expectedYN$values) +
                   1 / pmax(1, expectedNY$values) +
                   1 / pmax(1, expectedNN$values) )

    if (actOnCells) {
      normFact <- 1 / sqrt(getNumGenes(objCOTAN)) # divided by sqrt(n)
    }
    else {
      normFact <- 1 / sqrt(getNumCells(objCOTAN)) # divided by sqrt(m)
    }

    rm(expectedYN); rm(expectedNY); rm(expectedNN);
    gc()

    # coex estimation
    coex <- (observedYY$values - expectedYY$values)
    coex <- coex * sqrt(sumForDiv) * normFact

    coex <- list("genes" = observedYY$genes, "values" = coex)

    if(actOnCells) {
      objCOTAN@cellsCoex <- coex
    }
    else {
      objCOTAN@coex <- coex
    }

    rm(expectedYY)
    rm(coex)
    gc()

    print(paste("Calculate", kind, "coex: DONE"))
    return(objCOTAN)
  }
)
