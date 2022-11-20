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

    cat("calculating YY..")
    if (isTRUE(actOnCells)) {
      # for cells
      YY <- t(zeroOne) %*% zeroOne
    }
    else{
      # for genes
      YY <- zeroOne %*% t(zeroOne)
    }
    rm(zeroOne)
    gc()
    cat(" done\n")

    #return(as(YY, "symmetricMatrix"))
    return(as.matrix(YY))
  }
)


#' observedContingency
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
#' @rdname observedContingency
setMethod(
  "observedContingency",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE) {
    zeroOne <- getZeroOneProj(object)

    numGenes <- getNumGenes(objCOTAN)
    numCells <- getNumCells(objCOTAN)

    cat("calculating YY..")
    if (isTRUE(actOnCells)) {
      # dimension m x m (m number of cells)
      observedYY <- t(zeroOne) %*% zeroOne

      cat("NY..")
      # Any/Yes vector [cycled] = Yes/Yes + No/Yes
      observedAY <- colSums(probZero)

      observedNY <- observedAY - observedYY

      cat("YN..")
      observedYN <- t(observedNY)

      cat("NN..")
      observedNN <- numGenes - observedAY - observedYN
    }
    else {
      # dimension n x n (n number of genes)
      observedYY <- zeroOne %*% t(zeroOne)

      cat("YN..")
      # Yes/Any vector [cycled] = Yes/Yes + Yes/No
      observedYA <- rowSums(zeroOne)

      observedYN <- observedYA - observedYY

      cat("NY..")
      observedNY <- t(observedYN)

      cat("NN..")
      observedNN <- numCells - observedYA - observedNY
    }
    rm(probZero)
    gc()
    cat(" done\n")

    return( list("observedNN" = as.matrix(observedNN),
                 "observedNY" = as.matrix(observedNY),
                 "observedYN" = as.matrix(observedYN),
                 "observedYY" = as.matrix(observedYY)) )
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
#' @importFrom Rfast mat.mult
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
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

    numGenes <- getNumGenes(objCOTAN)
    numCells <- getNumCells(objCOTAN)

    cat("calculating NN..")
    if (isTRUE(actOnCells)) {
      # dimension m x m (m number of cells)
      expectedNN <- mat.mult(t(probZero) %*% probZero)

      cat("YN..")
      # Any/No vector [cycled] = No/No + Yes/No
      expectedAN <- colsums(probZero)

      expectedYN <- expectedNA - expectedNN

      cat("NY..")
      expectedNY <- t(expectedYN)

      cat("YY..")
      expectedYY <- numGenes - expectedNA - expectedNY
    }
    else {
      # dimension n x n (n number of genes)
      expectedNN <- mat.mult(probZero, t(probZero))

      cat("NY..")
      # No/Any vector [cycled] = No/No + No/Yes
      expectedNA <- rowsums(probZero)

      expectedNY <- expectedNA - expectedNN

      cat("YN..")
      expectedYN <- t(expectedNY)

      cat("YY..")
      expectedYY <- numCells - expectedNA - expectedYN
    }
    rm(probZero)
    gc()
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

    if (isFALSE(actOnCells)) {
      # remove the housekeeping genes
      noHKFlags <- flagNotHousekeepingGenes(objCOTAN)
      observedYY <- observedYY[noHKFlags, noHKFlags]
    }

    observedYY <- mat2vec_rfast(observedYY)

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

    if (isTRUE(actOnCells)) {
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

    # TODO: the "genes" name is misleading in case of actOnCells
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


#' calculateS
#'
#' calculate the statistics S for genes contingency tables
#'
#' @param objCOTAN A COTAN object
#'
#' @return the S matrix
#'
#' @export
#'
#' @rdname calculateS
#'
setMethod(
  "calculateS",
  "COTAN",
  function(objCOTAN) {
    print("Calculating S: START")

    objCOTAN <- standardizeCoex(objCOTAN)
    coex <- getCoex(objCOTAN, asMatrix = FALSE)

    stopifnot("Coex is missing" = !is_empty(coex))

    S <- (coex$values)^2 * getNumCells(objCOTAN)

    print("Calculating S: DONE")
    return( vec2mat_rfast( list("genes" = coex$genes, "values" = S) ) )
  }
)


#' calculateG
#'
#' calculate the statistics G-test for genes contingency tables
#' It is proportional to the genes' precence mutual information.
#'
#' @param objCOTAN A COTAN object
#'
#' @return the G matrix
#'
#' @export
#'
#' @rdname calculateG
#'
setMethod(
  "calculateG",
  "COTAN",
  function(objCOTAN) {
    print("Calculating G: START")

    list[observedNN, observedNY, observedYN, observedYY] <-
      observedContingency(objCOTAN, actOnCells = FALSE)

    noHKFlags <- flagNotHousekeepingGenes(objCOTAN)
    observedNN <- observedNN[noHKFlags, noHKFlags]
    observedNY <- observedNY[noHKFlags, noHKFlags]
    observedYN <- observedYN[noHKFlags, noHKFlags]
    observedYY <- observedYY[noHKFlags, noHKFlags]

    list[expectedNN, expectedNY, expectedYN, expectedYY] <-
      expectedContingencyTables(objCOTAN, actOnCells = FALSE)

    for (i in c(expectedNN, expectedNY, expectedYN, expectedYY)) {
      stopifnot("Some expected values are 0!" = !any(i == 0))
    }

    print("Estimating G")

    tNN <- observedNN * log(observedNN / expectedNN)
    tNN[which(observedNN == 0)] <- 0
    rm(observedNN); rm(expectedNN)

    tNY <- observedNY * log(observedNY / expectedNY)
    tNY[which(observedNY == 0)] <- 0
    rm(observedNY); rm(expectedNY)

    tYN <- observedYN * log(observedYN / expectedYN)
    tYN[which(observedYN == 0)] <- 0
    rm(observedYN); rm(expectedYN)

    tYY <- observedYY * log(observedYY / expectedYY)
    tYY[which(observedYY == 0)] <- 0
    rm(observedYY); rm(expectedYY)

    G <- 2 * (tNN + tNY + tYN + tYY)

    print("Calculating G: DONE")
    return( as.matrix(G) )
  }
)
