#' calculateMu
#'
#' @description Calculate the vector \eqn{\mu = \lambda \times \nu^T}
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns The Mu matrix
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix t
#'
#' @importClassesFrom Matrix dgeMatrix
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
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

    return(getLambda(objCOTAN) %o% getNu(objCOTAN))
  }
)


#' observedContingencyTablesYY
#'
#' @description Calculate observed Yes/Yes field of the contingency table
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean; when `TRUE` the function will return only
#'   packed dense symmetric matrices
#'
#' @returns A `list` with the 'Yes/Yes' observed contingency table and the 'Yes'
#'   observed vector
#'
#' @importFrom rlang is_empty
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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' observedYList <- observedContingencyTablesYY(objCOTAN, asDspMatrices = TRUE)
#' observedYY <- observedYList[[1]]
#' observedY  <- observedYList[[2]]
#'
#' @rdname observedContingencyTablesYY
#'
observedContingencyTablesYY <- function(objCOTAN,
                                        actOnCells = FALSE,
                                        asDspMatrices = FALSE) {
  zeroOne <- getZeroOneProj(objCOTAN)

  logThis("calculating YY..", logLevel = 3, appendLF = FALSE)
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

  logThis(" done", logLevel = 3)

  return(list("observedYY" = observedYY, "observedY" = observedY))
}



#' observedContingencyTables
#'
#' @description Calculate the observed contingency tables
#'
#' @details When `asDspMatrices == TRUE`, the method will effectively throw away
#'   the lower half from the returned `observedYN` and `observedNY` matrices,
#'   but, since they are transpose one of another, the full information is still
#'   available.
#'
#' @note The sum of the returned matrices will have constant value. This value
#'   is the number of genes/cells depending on `actOnCells` being
#'   `TRUE`/`FALSE`.
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean; when `TRUE` the function will return only
#'   packed dense symmetric matrices
#'
#' @returns The observed contingency tables as named `list` with elements:
#'   "observedNN", "observedNY", "observedYN" "observedYY"
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix symmetricMatrix
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' obs <- observedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' @rdname observedContingencyTables
#'
observedContingencyTables <- function(objCOTAN,
                                      actOnCells = FALSE,
                                      asDspMatrices = FALSE) {
  zeroOne <- getZeroOneProj(objCOTAN)

  numGenes <- getNumGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  list[observedYY, observedY] <-
    observedContingencyTablesYY(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = FALSE)
  gc()

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    logThis("calculating NY..", logLevel = 3, appendLF = FALSE)

    # Any/Yes vector [cycled] = Yes/Yes + No/Yes
    #observedNY <- observedY - observedYY

    logThis("YN..", logLevel = 3, appendLF = FALSE)
    observedYN <- t(observedY - observedYY)

    logThis("NN..", logLevel = 3, appendLF = FALSE)
    observedNN <- numGenes - observedY - observedYN

    observedYN <- as(observedYN, "denseMatrix")
  } else {
    # dimension n x n (n number of genes)
    logThis("calculating YN..", logLevel = 3, appendLF = FALSE)

    # Yes/Any vector [cycled] = Yes/Yes + Yes/No
    #observedYN <- observedY - observedYY

    logThis("NY..", logLevel = 3, appendLF = FALSE)
    observedNY <- t(observedY - observedYY)

    logThis("NN..", logLevel = 3, appendLF = FALSE)
    observedNN <- numCells - observedY - observedNY

    observedNY <- as(observedNY, "denseMatrix")
  }
  rm(zeroOne)
  gc()
  logThis(" done", logLevel = 3)

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
  } else {
    if( isTRUE(actOnCells) ) {
      observedNY <- t(observedYN)
    } else {
      observedYN <- t(observedNY)
    }
  }

  return( list("observedNN" = observedNN,
               "observedNY" = observedNY,
               "observedYN" = observedYN,
               "observedYY" = observedYY) )
}


#' expectedContingencyTablesNN
#'
#' @description Calculate the expected No/No field of the contingency table
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean; when `TRUE` the function will return only
#'   packed dense symmetric matrices
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use [Rfast]
#'   parallel algorithms that on the flip side use more memory
#'
#' @returns a `list` with the 'No/No' expected contingency table and the 'No'
#'   expected vector
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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#' expectedNList <- expectedContingencyTablesNN(objCOTAN, asDspMatrices = TRUE)
#' expectedNN <- expectedNList[[1]]
#' expectedN  <- expectedNList[[2]]
#'
#' @rdname expectedContingencyTablesNN
#'
expectedContingencyTablesNN <- function(objCOTAN,
                                        actOnCells = FALSE,
                                        asDspMatrices = FALSE,
                                        optimizeForSpeed = TRUE) {
  # estimate Probabilities of 0 with internal function funProbZero
  probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))
  gc()

  stopifnot("Error: some NA in matrix of probability of zero UMI counts" = !anyNA(probZero))

  numGenes <- getNumGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  logThis("calculating NN..", logLevel = 3, appendLF = FALSE)

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
    rownames(expectedNN) <- colnames(expectedNN) <- getCells(objCOTAN)
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
    rownames(expectedNN) <- colnames(expectedNN) <- getGenes(objCOTAN)
  }
  rm(probZero)

  expectedNN <- forceSymmetric(as(expectedNN, "denseMatrix"))

  if (isTRUE(asDspMatrices)) {
    expectedNN <- pack(expectedNN)
  }

  logThis(" done", logLevel = 3)

  return( list("expectedNN" = expectedNN, "expectedN" = expectedN) )
}


#' expectedContingencyTables
#'
#' @description Calculate the expected values of contingency tables
#'
#' @details When `asDspMatrices == TRUE`, the method will effectively throw away
#'   the lower half from the returned `expectedYN` and `expectedNY` matrices,
#'   but, since they are transpose one of another, the full information is still
#'   available.
#'
#' @note The sum of the returned matrices will have constant value. This value
#'   is the number of genes/cells depending on `actOnCells` being
#'   `TRUE`/`FALSE`.
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean; when `TRUE` the function will return only
#'   packed dense symmetric matrices
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use `Rfast`
#'   parallel algorithms that on the flip side use more memory
#'
#' @return The expected contingency tables as named `list` with elements:
#' "expectedNN", "expectedNY", "expectedYN", "expectedYY"
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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#' exp <- expectedContingencyTables(objCOTAN, asDspMatrices = TRUE)
#'
#' @rdname expectedContingencyTables
#'
expectedContingencyTables <- function(objCOTAN,
                                      actOnCells = FALSE,
                                      asDspMatrices = FALSE,
                                      optimizeForSpeed = TRUE) {
  # estimate Probabilities of 0 with internal function funProbZero
  probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))

  stopifnot("Error: some NA in matrix of probability of zero UMI counts" = !anyNA(probZero))

  numGenes <- getNumGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  list[expectedNN, expectedN] <-
    expectedContingencyTablesNN(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = isFALSE(optimizeForSpeed),
                                optimizeForSpeed = optimizeForSpeed)
  gc()

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    logThis("calculating YN..", logLevel = 3, appendLF = FALSE)

    # Any/No vector [cycled] = No/No + Yes/No
    #expectedYN <- expectedN - expectedNN

    logThis("NY..", logLevel = 3, appendLF = FALSE)
    expectedNY <- t(expectedN - expectedNN)

    logThis("YY..", logLevel = 3, appendLF = FALSE)
    expectedYY <- numGenes - expectedN - expectedNY

    expectedNY <- as(expectedNY, "denseMatrix")
  } else {
    # dimension n x n (n number of genes)
    logThis("calculating NY..", logLevel = 3, appendLF = FALSE)

    # No/Any vector [cycled] = No/No + No/Yes
    #expectedNY <- expectedN - expectedNN

    logThis("YN..", logLevel = 3, appendLF = FALSE)
    expectedYN <- t(expectedN - expectedNN)

    logThis("YY..", logLevel = 3, appendLF = FALSE)
    expectedYY <- numCells - expectedN - expectedYN

    expectedYN <- as(expectedYN, "denseMatrix")
  }
  rm(probZero)
  gc()

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
  } else {
    if( isTRUE(actOnCells) ) {
      expectedYN <- t(expectedNY)
    } else {
      expectedNY <- t(expectedYN)
    }
  }
  logThis(" done", logLevel = 3)

  return( list("expectedNN" = expectedNN,
               "expectedNY" = expectedNY,
               "expectedYN" = expectedYN,
               "expectedYY" = expectedYY) )
}



#' contingencyTables
#'
#' @description This function returns the observed and expected contingency
#'   tables for a given pair of genes.
#'
#' @details This code mimics the code to calculate the contingency tables from
#'   [observedContingencyTables()] and [expectedContingencyTables()], but
#'   restricted to only the relevant genes and thus much faster and less memory
#'   intensive
#'
#' @param objCOTAN a `COTAN` object
#' @param g1 a gene
#' @param g2 another gene
#'
#' @return A list containing the observed and expected contingency tables
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#' g1 <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 1)]
#' g2 <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 1)]
#' ctList <- contingencyTables(objCOTAN, g1 = g1, g2 = g2)
#' ctList[["observed"]]
#'
#' @rdname contingencyTables
#'
contingencyTables <- function(objCOTAN, g1, g2) {
  genes <- getGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  stopifnot("the first gene is unknown"  = (c(g1) %in% genes))
  stopifnot("the second gene is unknown" = (c(g2) %in% genes))

  dimnames <- list(c(paste(g2, "yes", sep = "."), paste(g2, "no", sep = ".")),
                   c(paste(g1, "yes", sep = "."), paste(g1, "no", sep = ".")))

  #-------------------------------------------------

  zeroOne <- sign(getRawData(objCOTAN)[c(g1, g2), ])
  rownames(zeroOne) <- c(g1, g2)

  observedYY <- tcrossprod(zeroOne)
  observedY  <- rowSums(zeroOne)

  rm(zeroOne)
  gc()

  observedNY <- t(observedY - observedYY)
  observedNN <- numCells - observedY - observedNY

  observedCT <- matrix(c(observedYY[g1, g2], observedNY[g2, g1],
                         observedNY[g1, g2], observedNN[g1, g2]),
                       ncol = 2, nrow = 2, dimnames = dimnames)

  #-------------------------------------------------

  probZero <- funProbZero(getDispersion(objCOTAN)[c(g1, g2)],
                          getLambda(objCOTAN)[c(g1, g2)] %o% getNu(objCOTAN))
  rownames(probZero) <- c(g1, g2)

  if (anyNA(probZero)) {
    warning("Some NA in estimated probability of zero matrix")
  }

  expectedNN <- tcrossprod(probZero)
  expectedN <- rowSums(probZero)

  rm(probZero)
  gc()

  expectedYN <- t(expectedN - expectedNN)
  expectedYY <- numCells - expectedN - expectedYN

  expectedCT <- matrix(c(expectedYY[g1, g2], expectedYN[g1, g2],
                         expectedYN[g2, g1], expectedNN[g1, g2]),
                       ncol = 2, nrow = 2, dimnames = dimnames)

  #-------------------------------------------------

  return(list("observed" = observedCT, "expected" = expectedCT))
}


#' calculateCoex
#'
#' @description This function estimates and stores the *coex* matrix in the
#'   `cellCoex` or `genesCoex` field depending on given `actOnCells` flag.
#'
#' @details The function calculates the coex matrix, defined by following
#'   formula:
#'
#'   \deqn{\frac{\sum_{i,j \in \{\text{Y, N}\}}{(-1)^{\#\{i,j\}}\frac{O_{ij}-E_{ij}}{1 \vee E_{ij}}}}
#'              {\sqrt{n \sum_{i,j \in \{\text{Y, N}\}}{\frac{1}{1 \vee E_{ij}}}}}}
#'
#'   where \eqn{O} and \eqn{E} are the observed and expected contingency tables
#'   and \eqn{n} is the relevant numerosity (the number of genes/cells depending
#'   on given `actOnCells` flag).
#'
#'   The formula can be more effectively implemented as:
#'
#'   \deqn{\sqrt{\frac{1}{n}\sum_{i,j \in \{\text{Y, N}\}}{\frac{1}{1 \vee E_{ij}}}}
#'         \, \bigl(O_\text{YY}-E_\text{YY}\bigr)}
#'
#'   once one notices that \eqn{O_{ij} - E_{ij} = (-1)^{\#\{i,j\}} \, r} for all
#'   \eqn{i,j \in \{\text{Y, N}\}}.
#'
#'
#'   The latter follows from the fact that the relevant marginal sums of the
#'   the expected contingency tables were enforced to match the marginal sums
#'   of the observed ones.
#'
#'   It also calculates the percentage of *problematic* genes/cells pairs.
#'   A pair is *problematic* when one or more of the expected counts were
#'   significantly smaller than 1 (\eqn{< 0.5}). These small expected values
#'   signal that scant information is present for such a pair.
#'
#' @seealso [estimateDispersionNuBisection()] for more details.
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use [Rfast]
#'   parallel algorithms that on the flip side use more memory
#'
#' @returns The updated `COTAN` object
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = TRUE)
#'
#' @rdname calculateCoex
#'
setMethod(
  "calculateCoex",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE,
           optimizeForSpeed = TRUE) {
    if (isTRUE(actOnCells)) { kind <- "cells'" } else { kind <- "genes'" }
    logThis(paste("Calculate", kind, "coex: START"), logLevel = 1)

    logThis(paste("Retrieving expected", kind, "contingency table"), logLevel = 3)

    # four estimators: expectedNN, expectedNY, expectedYN, expectedYY
    list[expectedNN, expectedNY, expectedYN, expectedYY] <-
      expectedContingencyTables(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = TRUE,
                                optimizeForSpeed = optimizeForSpeed)

    gc()

    logThis(paste("Calculating", kind, "coex normalization factor"), logLevel = 3)

    if (isTRUE(actOnCells)) {
      allNames <- getCells(objCOTAN)
      normFact <- 1 / sqrt(getNumGenes(objCOTAN)) # divided by sqrt(n)
    } else {
      allNames <- getGenes(objCOTAN)
      normFact <- 1 / sqrt(getNumCells(objCOTAN)) # divided by sqrt(m)
    }

    problematicPairsFraction <-
      sum(pmin(expectedYY@x, expectedYN@x,
               expectedNY@x, expectedNN@x) < 0.5) / length(expectedYY@x)

    coex <- normFact * sqrt(1 / pmax(1, expectedYY@x) +
                            1 / pmax(1, expectedYN@x) +
                            1 / pmax(1, expectedNY@x) +
                            1 / pmax(1, expectedNN@x) )

    rm(expectedYN); rm(expectedNY); rm(expectedNN);
    gc()

    logThis(paste("Retrieving observed", kind, "yes/yes contingency table"), logLevel = 3)

    list[observedYY, ] <-
      observedContingencyTablesYY(objCOTAN,
                                  actOnCells = actOnCells,
                                  asDspMatrices = TRUE)

    gc()

    # coex estimation
    logThis(paste("Estimating", kind, "coex"), logLevel = 3)

    coex <- coex * (observedYY@x - expectedYY@x)

    stopifnot("Matrix@x has the wrong size" = length(coex) == (length(allNames)*(length(allNames)+1)/2))

    coex <- new("dspMatrix", Dim = dim(expectedYY), x = coex)
    rownames(coex) <- colnames(coex) <- allNames

    rm(observedYY)
    rm(expectedYY)
    gc()

    if(actOnCells) {
      objCOTAN@cellsCoex <- coex
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[6], TRUE)
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[8],
                                             problematicPairsFraction)
    } else {
      objCOTAN@genesCoex <- coex
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[5], TRUE)
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[7],
                                             problematicPairsFraction)
    }

    rm(coex)
    gc()

    logThis(paste("Calculate", kind, "coex: DONE"), logLevel = 1)

    return(objCOTAN)
  }
)


#' calculateS
#'
#' @description This function calculates the statistics S for genes contingency
#'   tables
#'
#' @details It always has the diagonal set to zero.
#'
#' @param objCOTAN a `COTAN` object
#' @param geneSubsetCol an array of genes. It will be put in columns. If left
#'   empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows. If left empty
#'   the function will do it genome-wide.
#'
#' @return the S matrix
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)
#' S <- calculateS(objCOTAN)
#'
#' @rdname calculateS
#'
calculateS <- function(objCOTAN, geneSubsetCol = c(), geneSubsetRow = c()) {
  geneSubsetCol <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetCol)
  geneSubsetRow <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetRow)

  logThis("Calculating S: START", logLevel = 2)

  coex <- getGenesCoex(objCOTAN, zeroDiagonal = TRUE)
  if(!identical(geneSubsetRow, getGenes(objCOTAN)) ||
     !identical(geneSubsetCol, getGenes(objCOTAN))) {
    coex <- coex[geneSubsetRow, geneSubsetCol, drop = FALSE]
  }

  stopifnot("Coex is missing" = !is_empty(coex))

  S <- (coex)^2 * getNumCells(objCOTAN)

  logThis("Calculating S: DONE", logLevel = 2)

  return(S)
}



#' calculateG
#'
#' @description This function calculates the statistics G-test for genes
#'   contingency tables.
#'
#' @details  It always has the diagonal set to zero. It is proportional to the
#'   genes' presence mutual information.
#'
#' @param objCOTAN a `COTAN` object
#' @param geneSubsetCol an array of genes. It will be put in columns. If left
#'   empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows. If left empty
#'   the function will do it genome-wide.
#'
#' @return the G matrix
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)
#' G <- calculateG(objCOTAN)
#'
#' @rdname calculateG
#'
calculateG <- function(objCOTAN, geneSubsetCol = c(), geneSubsetRow = c()) {
  geneSubsetCol <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetCol)
  geneSubsetRow <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetRow)

  logThis("Calculating G: START", logLevel = 2)

  list[observedNN, observedNY, observedYN, observedYY] <-
    observedContingencyTables(objCOTAN, actOnCells = FALSE)

  list[expectedNN, expectedNY, expectedYN, expectedYY] <-
    expectedContingencyTables(objCOTAN, actOnCells = FALSE)

  logThis("Estimating G", logLevel = 3)

  tNN <- observedNN * log(observedNN / pmax(1, expectedNN))
  tNN[observedNN == 0] <- 0
  diag(tNN) <- 0
  tNN <- tNN[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedNN); rm(expectedNN)

  tNY <- observedNY * log(observedNY / pmax(1, expectedNY))
  tNY[observedNY == 0] <- 0
  diag(tNY) <- 0
  tNY <- tNY[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedNY); rm(expectedNY)

  tYN <- observedYN * log(observedYN / pmax(1, expectedYN))
  tYN[observedYN == 0] <- 0
  diag(tYN) <- 0
  tYN <- tYN[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedYN); rm(expectedYN)

  tYY <- observedYY * log(observedYY / pmax(1, expectedYY))
  tYY[observedYY == 0] <- 0
  diag(tYY) <- 0
  tYY <- tYY[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedYY); rm(expectedYY)
  gc()

  G <- 2 * (tNN + tNY + tYN + tYY)

  if(identical(geneSubsetRow, getGenes(objCOTAN)) &&
     identical(geneSubsetCol, getGenes(objCOTAN))) {
    G <- pack(forceSymmetric(G))
  }

  rm(tNN); rm(tNY); rm(tYN); rm(tYY)
  gc()

  logThis("Calculating G: DONE", logLevel = 2)

  return(G)
}


#' calculateGDI
#'
#' @description This function produces a`data.frame` with the GDI for each
#'   gene.
#'
#' @param objCOTAN a `COTAN` object
#' @param statType Which statistics to use to compute the p-values. By default
#'   it will use the "S" (Pearson's \eqn{\chi^{2}} test) otherwise the "G"
#'   (G-test)
#'
#' @returns A `data.frame` with the GDI data
#'
#' @export
#'
#' @importFrom rlang set_names
#'
#' @importFrom stats pchisq
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' GDI <- calculateGDI(objCOTAN)
#'
#' @rdname calculateGDI
#'
calculateGDI <- function(objCOTAN, statType = "S") {
  if (statType == "S") {
    logThis("Using S", logLevel = 3)
    S <- calculateS(objCOTAN)
  } else if (statType == "G") {
    logThis("Using G", logLevel = 3)
    S <- calculateG(objCOTAN)
  } else {
    stop("Unrecognised stat type: must be either 'S' or 'G'")
  }

  logThis("Calculate GDI dataframe: START", logLevel = 2)

  top5pcRows <- as.integer(max(1:round(getNumGenes(objCOTAN) / 20, digits = 0)))

  pValueSorted <- apply(S, 2, sort, decreasing = TRUE)
  pValueSorted <- pValueSorted[1:top5pcRows, , drop = FALSE]
  pValueSorted <- pchisq(as.matrix(pValueSorted), df = 1, lower.tail = FALSE)

  GDI <- log(-log(colMeans(pValueSorted)))
  GDI <- set_names(as.data.frame(GDI), "GDI")
  rm(pValueSorted)
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

  logThis("Calculate GDI dataframe: DONE", logLevel = 2)

  return(GDI)
}


#' calculatePValue
#'
#' @description This function computes the p-values for genes in the COTAN
#'   object. It can be used genome-wide or by setting some specific genes of
#'   interest. By default it computes the p-values using the S statistics
#'   (\eqn{\chi^{2}})
#'
#' @param objCOTAN a `COTAN` object
#' @param statType Which statistics to use to compute the p-values. By default
#'   it will use the "S" (Pearson's \eqn{\chi^{2}} test) otherwise the "G"
#'   (G-test)
#' @param geneSubsetCol an array of genes. It will be put in columns. If left
#'   empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows. If left empty
#'   the function will do it genome-wide.
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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#' objCOTAN <- calculateCoex(objCOTAN, actOnCells = FALSE)
#' pValue <- calculatePValue(objCOTAN)
#'
#' @rdname calculatePValue
#'
calculatePValue <- function(objCOTAN, statType = "S",
                            geneSubsetCol = c(), geneSubsetRow = c()) {
  geneSubsetCol <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetCol)
  geneSubsetRow <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetRow)

  if (statType == "S") {
    logThis("Using S", logLevel = 3)
    S <- calculateS(objCOTAN,
                    geneSubsetCol = geneSubsetCol,
                    geneSubsetRow = geneSubsetRow)
  } else if (statType == "G") {
    logThis("Using G", logLevel = 3)
    S <- calculateG(objCOTAN,
                    geneSubsetCol = geneSubsetCol,
                    geneSubsetRow = geneSubsetRow)
  } else {
    stop("Unrecognised stat type: must be either 'S' or 'G'")
  }

  logThis("calculating PValues: START", logLevel = 2)

  allCols <- identical(geneSubsetCol, getGenes(objCOTAN))
  allRows <- identical(geneSubsetRow, getGenes(objCOTAN))

  strCol <- if (allCols) "genome wide" else "on a set of genes"
  strRow <- if (allRows) "genome wide" else "on a set of genes"
  logThis(paste("Get p-values", strCol, "on columns and", strRow, "on rows"), logLevel = 2)

  pValues <- pchisq(as.matrix(S), df = 1, lower.tail = FALSE)

  if (allCols && allRows) {
    pValues <- pack(forceSymmetric(pValues))
  }

  logThis("calculating PValues: DONE", logLevel = 2)

  return(pValues)
}

