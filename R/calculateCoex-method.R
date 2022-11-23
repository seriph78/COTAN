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
    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(getDispersion(objCOTAN), estimateMu(objCOTAN))

    errMgs <- "Error: some NA in matrix of probability of zero UMI counts. "
    stopifnot(errMgs = !anyNA(probZero))

    numGenes <- getNumGenes(objCOTAN)
    numCells <- getNumCells(objCOTAN)

    cat("calculating NN..")
    if (isTRUE(actOnCells)) {
      # dimension m x m (m number of cells)
      expectedNN <- mat.mult(t(probZero), probZero)

      cat("YN..")
      # Any/No vector [cycled] = No/No + Yes/No
      expectedAN <- colsums(probZero)

      expectedYN <- expectedAN - expectedNN

      cat("NY..")
      expectedNY <- t(expectedYN)

      cat("YY..")
      expectedYY <- numGenes - expectedAN - expectedNY
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
#' This function estimates and stores the coex matrix in the 'genesCoex' or
#' 'cellCoex' field depending on given 'actOnCells' flag
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
    if(actOnCells) { kind <- "cells'" } else { kind <- "genes'" }
    print(paste("Calculate", kind, "coex: START"))

    print("Retrieving observed yes/yes contingency table")

    observedYY <- observedContingencyYY(objCOTAN, actOnCells)

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

    print(paste("Estimating", kind, "coex"))

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
    coex <- list("genes" = getGenes(objCOTAN), "values" = coex)

    if(actOnCells) {
      objCOTAN@cellsCoex <- coex
    }
    else {
      objCOTAN@genesCoex <- coex
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
    coex <- getGenesCoex(objCOTAN, asMatrix = FALSE)

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


#' calculateGDI
#'
#' This function produce a dataframe with the GDI for each genes.
#'
#' @param object A COTAN object
#' @param type Type of statistic to be used. Default is "S":
#' Pearson's chi-squared test statistics. "G" is G-test statistics
#'
#' @return A dataframe
#'
#' @export
#'
#' @importFrom stats pchisq
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#'
#' @examples
#' data("ERCC.cotan")
#' quant.p <- calculateGDI(ERCC.cotan)
#'
#' @rdname calculateGDI
setMethod(
  "calculateGDI",
  "COTAN",
  function(objCOTAN, type = "S") {
    print("Calculate GDI dataframe: START")

    if (type == "S") {
      print("Using S")
      S <- calculateS(objCOTAN)
    } else if (type == "G") {
      print("Using G")
      S <- calculateG(objCOTAN)
    }

    diag(S) <- 0
    CDSorted <- apply(S, 2, sort, decreasing = TRUE)
    rg <- round(nrow(CDSorted) / 20, digits = 0)
    CDSorted <- CDSorted[seq_len(rg), ]
    CDSorted <- pchisq(as.matrix(CDSorted), df = 1, lower.tail = FALSE)

    GDI <- log(-log(colMeans(CDSorted)))
    GDI <- as.data.frame(GDI)
    colnames(GDI) <- "GDI"
    rm(CDSorted)
    gc()

    sum.raw.norm <- log(rowSums(getNormalizedData(objCOTAN)))
    GDI <- merge(GDI, as.data.frame(sum.raw.norm), by = "row.names", all.x = TRUE)
    GDI <- column_to_rownames(GDI, var = "Row.names")
    rm(sum.raw.norm)
    gc()

    exp.cells <- (rowSums(getZeroOneProj(objCOTAN)) / getNumCells(objCOTAN)) * 100
    GDI <- merge(GDI, as.data.frame(exp.cells), by = "row.names", all.x = TRUE)
    GDI <- column_to_rownames(GDI, var = "Row.names")
    rm(exp.cells)
    gc()

    GDI <- GDI[, c("sum.raw.norm", "GDI", "exp.cells")]

    print("Calculate GDI dataframe: DONE")
    return(GDI)
  }
)
