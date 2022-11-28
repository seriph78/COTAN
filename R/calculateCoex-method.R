#' calculateMu
#'
#' estimate vector mu
#' @param objCOTAN a COTAN object
#'
#' @return the Mu matrix
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix dgeMatrix
#'
#' @export
#'
#' @rdname calculateMu
setMethod(
  "calculateMu",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getLambda(objCOTAN))) {
      stop("lambda must not be empty, estimate it")
    }

    if (is_empty(getNu(objCOTAN))) {
      stop("nu must not be empty, estimate it")
    }

    muEstimator <- getLambda(objCOTAN) %*% t(getNu(objCOTAN))

    rownames(muEstimator) <- getGenes(objCOTAN)

    return(as.matrix(muEstimator))
  }
)


#' observedContingencyYY
#'
#' calculate observed yes/yes field of contingency table
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells,
#' otherwise for the genes
#' @return the YesYes observed contengency table
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix dgeMatrix
#' @importClassesFrom Matrix dsyMatrix
#'
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

    YY <- forceSymmetric(as(YY, "generalMatrix"))

    cat(" done\n")

    return(YY)
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

    return( list("observedNN" = as(observedNN, "generalMatrix"),
                 "observedNY" = as(observedNY, "generalMatrix"),
                 "observedYN" = as(observedYN, "generalMatrix"),
                 "observedYY" = as(observedYY, "generalMatrix")) )
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
    probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

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

    out <- list( "expectedNN" = as(expectedNN, "generalMatrix"),
                 "expectedNY" = as(expectedNY, "generalMatrix"),
                 "expectedYN" = as(expectedYN, "generalMatrix"),
                 "expectedYY" = as(expectedYY, "generalMatrix") )

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

    print("Retrieving expected contingency table")

    # four estimators: expectedNN, expectedNY, expectedYN, expectedYY
    list[expectedNN, expectedNY, expectedYN, expectedYY] <-
      expectedContingencyTables(objCOTAN, actOnCells)

    print(paste("Estimating", kind, "coex"))

    # sum for division
    sumForDiv <- ( 1 / pmax(1, expectedYY) +
                   1 / pmax(1, expectedYN) +
                   1 / pmax(1, expectedNY) +
                   1 / pmax(1, expectedNN) )

    if (isTRUE(actOnCells)) {
      normFact <- 1 / sqrt(getNumGenes(objCOTAN)) # divided by sqrt(n)
    }
    else {
      normFact <- 1 / sqrt(getNumCells(objCOTAN)) # divided by sqrt(m)
    }

    rm(expectedYN); rm(expectedNY); rm(expectedNN);

    # coex estimation
    coex <- (observedYY - expectedYY)
    coex <- coex * sqrt(sumForDiv) * normFact

    rm(sumForDiv)
    rm(observedYY)
    rm(expectedYY)
    gc()

    coex <- pack(forceSymmetric(coex))

    if(actOnCells) {
      objCOTAN@cellsCoex <- coex
    }
    else {
      objCOTAN@genesCoex <- coex
    }

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

    coex <- getGenesCoex(objCOTAN)

    stopifnot("Coex is missing" = !is_empty(coex))

    S <- (coex)^2 * getNumCells(objCOTAN)

    print("Calculating S: DONE")

    return(S)
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
    gc()

    G <- 2 * (tNN + tNY + tYN + tYY)

    rm(tNN); rm(tNY); rm(tYN); rm(tYY)
    gc()

    print("Calculating G: DONE")
    return(G)
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


#' calculatePValue
#'
#' This function computes the p-values for genes in the COTAN object.
#' It can be used genome-wide or by setting some specific genes of interest.
#' By default it computes the p-values using the S statistics (\eqn{\chi^{2}})
#'
#' @param objCOTAN a COTAN object
#' @param statType Whic statistics to use to compute the p-values.
#' By default it will use the S (\eqn{\chi^{2}})
#' @param geneSubsetCol an array of genes. It will be put in columns.
#' If left empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows.
#' If left empty the function will do it genome-wide.
#'
#' @return a p-value matrix
#'
#' @export
#'
#' @importFrom Matrix forceSymmetric
#'
#' @examples
#'
#' data("ERCC.cotan")
#' ERCC.cotan <- calculatePValue(ERCC.cotan, statType = "S")
#'
#' @rdname calculatePValue
setMethod(
  "calculatePValue",
  "COTAN",
  function(objCOTAN, statType = "S", geneSubsetCol = c(), geneSubsetRow = c()) {
    if(is_empty(geneSubsetCol) && !is_empty(geneSubsetRow)) {
      stop(paste0("can't have genome wide on columns and not rows!",
                  " Use a subset on geneSubsetCol, not on rows."))
    }

    if (statType == "S") {
      print("Using function S")
      S <- calculateS(objCOTAN)
    }
    else if (statType == "G") {
      print("Using function G")
      S <- calculateG(objCOTAN)
    }
    else {
      stop("Unrecognised stat type: must be either 'S' or 'G'")
    }

    print("calculating PValues: START")

    genes <- getGenes(objCOTAN)
    if (is_empty(geneSubsetCol)) {
      geneSubsetCol <- genes
    }
    else {
      geneSubsetCol = genes[genes %in% geneSubsetCol]
    }

    if (is_empty(geneSubsetRow)) {
      geneSubsetRow <- geneSubsetCol
    }
    else {
      geneSubsetRow <- geneSubsetCol[geneSubsetCol %in% geneSubsetRow]
    }

    strCol <- (if (all(genes %in% geneSubsetCol)) "genome wide" else "on a set of genes")
    strRow <- (if (all(genes %in% geneSubsetRow)) "genome wide" else "on a set of genes")
    print(paste("Get p-values", strCol, "on columns and", strRow, "on rows"))

    S <- S[geneSubsetCol, geneSubsetRow, drop = FALSE]

    p_value <- pchisq(as.matrix(S), df = 1, lower.tail = FALSE)
    print("calculating PValues: DONE")

    return(p_value)
  }
)
