#' calculateMu
#'
#' calculate vector mu = lambda × nu^t
#'
#' @param objCOTAN a COTAN object
#'
#' @return the Mu matrix
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix dgeMatrix
#'
#' @export
#'
#' @examples
#' mu <- calculateMu(objCOTAN)
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


#' observedContingencyTablesYY
#'
#' calculate observed Yes/Yes field of contingency table
#'
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells,
#' otherwise for the genes
#' @param asDspMatrices Boolean, if TRUE the function will return only packed
#' dense symmetric matrices
#'
#' @return a list with the 'Yes/Yes' observed contingency table
#' and the 'Yes' observed vector
#'
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
#' @importFrom Matrix tcrossprod
#'
#' @importClassesFrom Matrix dgeMatrix
#' @importClassesFrom Matrix dsyMatrix
#'
#' @export
#'
#' @examples
#' obsYY <- observedContingencyTablesYY(objCOTAN, asDspMatrices = TRUE)
#'
#' @rdname observedContingencyTablesYY
setMethod(
  "observedContingencyTablesYY",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE) {
    zeroOne <- getZeroOneProj(objCOTAN)

    cat("calculating YY..")
    if (isTRUE(actOnCells)) {
      # for cells
      observedYY <- crossprod(zeroOne)
      observedY  <- colSums(zeroOne)
    } else {
      # for genes
      observedYY <- tcrossprod(zeroOne)
      observedY  <- rowSums(zeroOne)
    }
    rm(zeroOne)

    observedYY <- forceSymmetric(as(observedYY, "denseMatrix"))

    if (isTRUE(asDspMatrices)) {
      observedYY <- pack(observedYY)
    }

    cat(" done\n")

    return(list("observedYY" = observedYY, "observedY" = observedY))
  }
)


#' observedContingencyTables
#'
#' calculate observed yes/yes field of contingency table
#'
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells,
#' otherwise for the genes
#' #' @param asDspMatrices Boolean, if TRUE the function will return only packed
#' dense symmetric matrices
#'
#' @return the observed contingency tables
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix symmetricMatrix
#' @export
#'
#' @examples
#' obs <- observedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' @rdname observedContingencyTables
setMethod(
  "observedContingencyTables",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE) {
    zeroOne <- getZeroOneProj(object)

    numGenes <- getNumGenes(objCOTAN)
    numCells <- getNumCells(objCOTAN)

    list[observedYY, observedY] <-
      observedContingencyTablesYY(objCOTAN,
                                  actOnCells = actOnCells,
                                  asDspMatrices = FALSE)
    gc()

    if (isTRUE(actOnCells)) {
      # dimension m x m (m number of cells)
      cat("calculating NY..")

      # Any/Yes vector [cycled] = Yes/Yes + No/Yes
      #observedNY <- observedY - observedYY

      cat("YN..")
      observedYN <- t(observedY - observedYY)

      cat("NN..")
      observedNN <- numGenes - observedY - observedYN

      observedYN <- as(observedYN, "denseMatrix")
    } else {
      # dimension n x n (n number of genes)
      cat("calculating YN..")

      # Yes/Any vector [cycled] = Yes/Yes + Yes/No
      #observedYN <- observedY - observedYY

      cat("NY..")
      observedNY <- t(observedY - observedYY)

      cat("NN..")
      observedNN <- numCells - observedY - observedNY

      observedNY <- as(observedNY, "denseMatrix")
    }
    rm(probZero)
    gc()
    cat(" done\n")

    observedNN <- forceSymmetric(as(observedNN, "denseMatrix"))

    if (isTRUE(asDspMatrices)) {
      observedNN <- pack(observedNN)
      observedYY <- pack(observedYY)
      # these operation drops the lower triangle values
      # but the other matrix contains them anyway
      if( isTRUE(actOnCells) ) {
        observedNY <- pack(forceSymmetric(t(observedYN)))
        observedYN <- pack(forceSymmetric(  observedYN) )
      } else {
        observedYN <- pack(forceSymmetric(t(observedNY)))
        observedNY <- pack(forceSymmetric(  observedNY) )
      }
    }

    return( list("observedNN" = observedNN,
                 "observedNY" = observedNY,
                 "observedYN" = observedYN,
                 "observedYY" = observedYY) )
  }
)


#' expectedContingencyTablesNN
#'
#' calculate the expected No/No filed of contingency tables
#'
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells,
#' otherwise for the genes
#' @param asDspMatrices Boolean, if TRUE the function will return only packed
#' dense symmetric matrices
#' @param optimizeForSpeed Boolean, if TRUE, the function will use Rfast
#' parallel algorithms that on the flip side use more memory
#'
#' @return a list with the 'No/No' observed contingency table
#' and the 'No' expected vector
#'
#' @importFrom Matrix t
#'
#' @importClassesFrom Matrix dgeMatrix
#' @importClassesFrom Matrix dsyMatrix
#'
#' @importFrom Rfast Crossprod
#' @importFrom Rfast Tcrossprod
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @export
#'
#' @examples
#' expNN <- expectedContingencyTablesNN(objCOTAN, asDspMatrices = TRUE)
#'
#' @rdname expectedContingencyTablesNN
setMethod(
  "expectedContingencyTablesNN",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE, optimizeForSpeed = TRUE) {
    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))
    gc()

    stopifnot("Error: some NA in matrix of probability of zero UMI counts" = !anyNA(probZero))

    numGenes <- getNumGenes(objCOTAN)
    numCells <- getNumCells(objCOTAN)

    cat("calculating NN..")

    if (isTRUE(actOnCells)) {
      # dimension m x m (m number of cells)
      if (isTRUE(optimizeForSpeed)) {
        # use Rfast functions
        expectedNN <- Crossprod(probZero, probZero)
        expectedN  <- colsums(probZero, parallel = TRUE)
      } else {
        expectedNN <- crossprod(probZero)
        expectedN  <- colSums(probZero)
      }
    } else {
      # dimension n x n (n number of genes)
      if (isTRUE(optimizeForSpeed)) {
        # use Rfast functions
        expectedNN <- Tcrossprod(probZero, probZero)
        expectedN  <- rowsums(probZero, parallel = TRUE)
      } else {
        expectedNN <- tcrossprod(probZero)
        expectedN  <- rowSums(probZero)
      }
    }
    rm(probZero)

    expectedNN <- forceSymmetric(as(expectedNN, "denseMatrix"))

    if (isTRUE(asDspMatrices)) {
      expectedNN <- pack(expectedNN)
    }

    cat(" done\n")

    return( list("expectedNN" = expectedNN, "expectedN" = expectedN) )
  }
)


#' expectedContingencyTables
#'
#' calcualte the expected values of contingency tables
#'
#' @param objCOTAN A COTAN object
#' @param actOnCells Boolean, if TRUE, the function works for the cells,
#' otherwise for the genes
#' @param asDspMatrices Boolean, if TRUE the function will return only packed
#' dense symmetric matrices
#' @param optimizeForSpeed Boolean, if TRUE, the function will use Rfast
#' parallel algorithms that on the flip side use more memory
#'
#' @return a list with the expected contingency tables
#'
#' @importFrom Rfast Crossprod
#' @importFrom Rfast Tcrossprod
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix symmetricMatrix
#'
#' @export
#'
#' @examples
#' exp <- expectedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' @rdname expectedContingencyTables
setMethod(
  "expectedContingencyTables",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE, optimizeForSpeed = TRUE) {
    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

    stopifnot("Error: some NA in matrix of probability of zero UMI counts" = !anyNA(probZero))

    numGenes <- getNumGenes(objCOTAN)
    numCells <- getNumCells(objCOTAN)

    print(paste0("ECT: dims[", numGenes, ", ", numCells, "]"))

    list[expectedNN, expectedN] <-
      expectedContingencyTablesNN(objCOTAN,
                                  actOnCells = actOnCells,
                                  asDspMatrices = FALSE,
                                  optimizeForSpeed = optimizeForSpeed)
    gc()

    if (isTRUE(actOnCells)) {
      # dimension m x m (m number of cells)
      cat("calculating YN..")

      # Any/No vector [cycled] = No/No + Yes/No
      #expectedYN <- expectedN - expectedNN

      cat("NY..")
      expectedNY <- t(expectedN - expectedNN)

      cat("YY..")
      expectedYY <- numGenes - expectedN - expectedNY

      expectedNY <- as(expectedNY, "denseMatrix")
    } else {
      # dimension n x n (n number of genes)
      cat("calculating NY..")

      # No/Any vector [cycled] = No/No + No/Yes
      #expectedNY <- expectedN - expectedNN

      cat("YN..")
      expectedYN <- t(expectedN - expectedNN)

      cat("YY..")
      expectedYY <- numCells - expectedN - expectedYN

      expectedYN <- as(expectedYN, "denseMatrix")
    }
    rm(probZero)
    gc()
    cat(" done\n")

    expectedYY <- forceSymmetric(as(expectedYY, "denseMatrix"))

    if (isTRUE(asDspMatrices)) {
      expectedNN <- pack(expectedNN)
      expectedYY <- pack(expectedYY)
      # these operation drops the lower triangle values
      # but the other matrix contains them anyway
      if( isTRUE(actOnCells) ) {
        expectedYN <- pack(forceSymmetric(t(expectedNY)))
        expectedNY <- pack(forceSymmetric(  expectedNY) )
      } else {
        expectedNY <- pack(forceSymmetric(t(expectedYN)))
        expectedYN <- pack(forceSymmetric(  expectedYN) )
      }
    }

    return( list("expectedNN" = expectedNN,
                 "expectedNY" = expectedNY,
                 "expectedYN" = expectedYN,
                 "expectedYY" = expectedYY) )
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
#' @param optimizeForSpeed Boolean, if TRUE, the function will use Rfast
#' parallel algorithms that on the flip side use more memory
#'
#' @return the updated COTAN object
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)[["objCOTAN"]]
#' objCOTAN <- estimateDispersion(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#'
#' @rdname calculateCoex
setMethod(
  "calculateCoex",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE, optimizeForSpeed = TRUE) {
    if (isTRUE(actOnCells)) {
      kind <- "cells'"
      allNames <- getCells(objCOTAN)
      normFact <- 1 / sqrt(getNumGenes(objCOTAN)) # divided by sqrt(n)
    } else {
      kind <- "genes'"
      allNames <- getGenes(objCOTAN)
      normFact <- 1 / sqrt(getNumCells(objCOTAN)) # divided by sqrt(m)
    }

    print(paste("Calculate", kind, "coex: START"))

    print("Retrieving expected contingency table")

    # four estimators: expectedNN, expectedNY, expectedYN, expectedYY
    list[expectedNN, expectedNY, expectedYN, expectedYY] <-
      expectedContingencyTables(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = TRUE,
                                optimizeForSpeed = optimizeForSpeed)

    gc()

    print(paste("Calculating", kind, "coex normalization factor"))

    coex <- normFact * sqrt(1 / pmax(1, expectedYY@x) +
                            1 / pmax(1, expectedYN@x) +
                            1 / pmax(1, expectedNY@x) +
                            1 / pmax(1, expectedNN@x) )

    rm(expectedYN); rm(expectedNY); rm(expectedNN);
    gc()

    print("Retrieving observed yes/yes contingency table")

    list[observedYY, ] <-
      observedContingencyTablesYY(objCOTAN,
                                  actOnCells = actOnCells,
                                  asDspMatrices = TRUE)

    gc()

    # coex estimation
    print(paste("Estimating", kind, "coex"))

    coex <- coex * (observedYY@x - expectedYY@x)

    stopifnot("Matrix@x has the wrong size" = length(coex) == (length(allNames)*(length(allNames)+1)/2))

    coex <- new("dspMatrix", Dim = dim(expectedYY), x = coex)
    rownames(coex) <- colnames(coex) <- allNames

    rm(observedYY)
    rm(expectedYY)
    gc()

    if(actOnCells) {
      objCOTAN@cellsCoex <- coex
    } else {
      objCOTAN@genesCoex <- coex
    }

    rm(coex)
    gc()

    print(paste("Calculate", kind, "coex: DONE"))

    return(objCOTAN)
  }
)


# Internal function to handle genes subset...
handleGenesSubsets <- function(genes, genesSubset = c()) {
  if (is_empty(genesSubset)) {
    genesSubset <- genes
  } else {
    stopifnot("Passed genes are not a subset of the full list" = all(genesSubset %in% genes))
    genesSubset = genes[genes %in% genesSubset]
  }
  return(genesSubset)
}

#' calculateS
#'
#' calculate the statistics S for genes contingency tables
#'
#' @param objCOTAN A COTAN object
#' @param geneSubsetCol an array of genes. It will be put in columns.
#' If left empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows.
#' If left empty the function will do it genome-wide.
#'
#' @return the S matrix
#'
#' @export
#'
#' @examples
#' S <- calculateS(objCOTAN)
#'
#' @rdname calculateS
#'
setMethod(
  "calculateS",
  "COTAN",
  function(objCOTAN, geneSubsetCol = c(), geneSubsetRow = c()) {
    geneSubsetCol <- handleGenesSubsets(getGenes(objCOTAN), geneSubsetCol)
    geneSubsetRow <- handleGenesSubsets(getGenes(objCOTAN), geneSubsetRow)

    print("Calculating S: START")

    coex <- getGenesCoex(objCOTAN)[geneSubsetRow, geneSubsetCol, drop = FALSE]

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
#' @param geneSubsetCol an array of genes. It will be put in columns.
#' If left empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows.
#' If left empty the function will do it genome-wide.
#'
#' @return the G matrix
#'
#' @export
#'
#' @examples
#' G <- calculateG(objCOTAN)
#'
#' @rdname calculateG
#'
setMethod(
  "calculateG",
  "COTAN",
  function(objCOTAN, geneSubsetCol = c(), geneSubsetRow = c()) {
    geneSubsetCol <- handleGenesSubsets(getGenes(objCOTAN), geneSubsetCol)
    geneSubsetRow <- handleGenesSubsets(getGenes(objCOTAN), geneSubsetRow)

    print("Calculating G: START")

    list[observedNN, observedNY, observedYN, observedYY] <-
      observedContingencyTables(objCOTAN, actOnCells = FALSE)

    list[expectedNN, expectedNY, expectedYN, expectedYY] <-
      expectedContingencyTables(objCOTAN, actOnCells = FALSE)

    for (i in c(expectedNN, expectedNY, expectedYN, expectedYY)) {
      stopifnot("Some expected values are 0!" = !any(i == 0))
    }

    print("Estimating G")

    tNN <- observedNN * log(observedNN / expectedNN)
    tNN <- tNN[geneSubsetRow, geneSubsetCol, drop = FALSE]
    tNN[which(observedNN == 0)] <- 0
    rm(observedNN); rm(expectedNN)

    tNY <- observedNY * log(observedNY / expectedNY)
    tNY <- tNY[geneSubsetRow, geneSubsetCol, drop = FALSE]
    tNY[which(observedNY == 0)] <- 0
    rm(observedNY); rm(expectedNY)

    tYN <- observedYN * log(observedYN / expectedYN)
    tYN <- tYN[geneSubsetRow, geneSubsetCol, drop = FALSE]
    tYN[which(observedYN == 0)] <- 0
    rm(observedYN); rm(expectedYN)

    tYY <- observedYY * log(observedYY / expectedYY)
    tYY <- tYY[geneSubsetRow, geneSubsetCol, drop = FALSE]
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
#'
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)[["objCOTAN"]]
#' objCOTAN <- estimateDispersion(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' GDI <- calculateGDI(objCOTAN)
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
#' @return a p-value matrix as dspMatrix
#'
#' @export
#'
#' @importFrom Matrix forceSymmetric
#' @importFrom Matrix pack
#'
#' @importClassesFrom Matrix dspMatrix
#'
#' @examples
#'
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)[["objCOTAN"]]
#' objCOTAN <- estimateDispersion(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' pValue <- calculatePValue(objCOTAN)
#'
#' @rdname calculatePValue
setMethod(
  "calculatePValue",
  "COTAN",
  function(objCOTAN, statType = "S", geneSubsetCol = c(), geneSubsetRow = c()) {
    if (statType == "S") {
      print("Using function S")
      S <- calculateS(objCOTAN,
                      geneSubsetCol = geneSubsetCol,
                      geneSubsetRow = geneSubsetRow)
    } else if (statType == "G") {
      print("Using function G")
      S <- calculateG(objCOTAN,
                      geneSubsetCol = geneSubsetCol,
                      geneSubsetRow = geneSubsetRow)
    } else {
      stop("Unrecognised stat type: must be either 'S' or 'G'")
    }

    print("calculating PValues: START")

    strCol <- (if (all(getGenes(objCOTAN) %in% geneSubsetCol)) "genome wide" else "on a set of genes")
    strRow <- (if (all(getGenes(objCOTAN) %in% geneSubsetRow)) "genome wide" else "on a set of genes")
    print(paste("Get p-values", strCol, "on columns and", strRow, "on rows"))

    pValues <- pchisq(as.matrix(S), df = 1, lower.tail = FALSE)

    print("calculating PValues: DONE")

    return(pValues)
  }
)
