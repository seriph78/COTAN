# ----------- Contingency tables --------------

#' @description These are the functions and methods used to calculate the
#'   **COEX** matrices according to the `COTAN` model. From there it is possible
#'   to calculate the associated *pValue* and the *GDI* (*Global Differential
#'   Expression*)
#'
#' @description The **COEX** matrix is defined by following formula:
#'
#'   \deqn{\frac{\sum_{i,j \in \{\text{Y, N}\}}{
#'                     (-1)^{\#\{i,j\}}\frac{O_{ij}-E_{ij}}{1 \vee E_{ij}}}}
#'              {\sqrt{n \sum_{i,j \in \{\text{Y, N}\}}{
#'                             \frac{1}{1 \vee E_{ij}}}}}}
#'
#'   where \eqn{O} and \eqn{E} are the observed and expected contingency tables
#'   and \eqn{n} is the relevant numerosity (the number of genes/cells depending
#'   on given `actOnCells` flag).
#'
#'   The formula can be more effectively implemented as:
#'
#'   \deqn{\sqrt{\frac{1}{n}\sum_{i,j \in \{\text{Y, N}\}}{
#'                                \frac{1}{1 \vee E_{ij}}}}
#'         \, \bigl(O_\text{YY}-E_\text{YY}\bigr)}
#'
#'   once one notices that \eqn{O_{ij} - E_{ij} = (-1)^{\#\{i,j\}} \, r} for
#'   some constant \eqn{r} for all \eqn{i,j \in \{\text{Y, N}\}}.
#'
#'   The latter follows from the fact that the relevant marginal sums of the
#'   the expected contingency tables were enforced to match the marginal sums
#'   of the observed ones.
#'
#' @seealso [ParametersEstimations] for more details.
#'
#' @name CalculatingCOEX
#'
#' @rdname CalculatingCOEX
NULL

#' @aliases calculateMu
#'
#' @note The sum of the matrices returned by the function
#'   `observedContingencyTables()` and `expectedContingencyTables()` will have
#'   the same value on all elements. This value is the number of genes/cells
#'   depending on the parameter `actOnCells` being `TRUE/FALSE`.
#'
#' @details `calculateMu()` calculates the vector \eqn{\mu = \lambda \times
#'   \nu^T}
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `calculateMu()` returns the `mu` matrix
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix t
#'
#' @importClassesFrom Matrix dgeMatrix
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
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


#' @details `observedContingencyTablesYY()` calculates observed *Yes/Yes* field
#'   of the contingency table
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean - when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean - when `TRUE` the function will return only
#'   packed dense symmetric matrices
#'
#' @returns `observedContingencyTablesYY()` returns a `list` with:
#'   * `observedYY` the *Yes/Yes* observed contingency table as `matrix`
#'   * `observedY`  the full *Yes* observed `vector`
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
#' @importFrom Matrix tcrossprod
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#'
#' @importClassesFrom Matrix dgeMatrix
#' @importClassesFrom Matrix dsyMatrix
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
observedContingencyTablesYY <- function(objCOTAN,
                                        actOnCells = FALSE,
                                        asDspMatrices = FALSE) {
  logThis("calculating YY..", logLevel = 3L, appendLF = FALSE)

  zeroOne <- getZeroOneProj(objCOTAN)

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

  logThis(" done", logLevel = 3L)

  return(list("observedYY" = observedYY, "observedY" = observedY))
}


#' @details `observedPartialContingencyTablesYY()` calculates observed *Yes/Yes*
#'   field of the contingency table
#'
#' @param objCOTAN a `COTAN` object
#' @param columnsSubset a sub-set of the columns of the matrices that will be
#'   returned
#' @param zeroOne the raw count matrix projected to `0` or `1`. If not given the
#'   appropriate one will be calculated on the fly
#' @param actOnCells Boolean - when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean - when `TRUE` the function will return only
#'   packed dense symmetric matrices
#'
#' @returns `observedPartialContingencyTablesYY()` returns a `list` with:
#'   * `observedYY` the *Yes/Yes* observed contingency table as `matrix`,
#'     restricted to the selected columns as named `list` with elements
#'   * `observedY`  the full *Yes* observed `vector`
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
#' @importFrom Matrix tcrossprod
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
observedPartialContingencyTablesYY <-
  function(objCOTAN,
           columnsSubset,
           zeroOne = NULL,
           actOnCells = FALSE) {
  logThis("calculating partial YY..", logLevel = 3L, appendLF = FALSE)

  if (is_empty(zeroOne)) {
    # get zero/one projection with internal function get
    zeroOne <- getZeroOneProj(objCOTAN)
  }
  assert_that(identical(dim(zeroOne), dim(getRawData(objCOTAN))))
  assert_that(!anyNA(zeroOne),
              msg = paste0("Error: some NA in matrix of probability",
                           " of zero UMI counts"))
  gc()

  if (is.character(columnsSubset)) {
    if (isTRUE(actOnCells)) {
      allColumns <- getCells(objCOTAN)
    } else {
      allColumns <- getGenes(objCOTAN)
    }
    columnsSubset <- which(allColumns %in% columnsSubset)
  } else {
    columnsSubset <- sort(columnsSubset)
  }

  if (isTRUE(actOnCells)) {
    # for cells
    observedYY <- crossprod(zeroOne,
                            zeroOne[, columnsSubset, drop = FALSE])
    observedY  <- colSums(zeroOne)
  } else {
    # for genes
    observedYY <- tcrossprod(zeroOne,
                             zeroOne[columnsSubset, , drop = FALSE])
    observedY  <- rowSums(zeroOne)
  }
  rm(zeroOne)

  logThis(" done", logLevel = 3L)

  return(list("observedYY" = observedYY, "observedY" = observedY))
}


#' @details `observedContingencyTables()` calculates the observed contingency
#'   tables. When the parameter `asDspMatrices == TRUE`, the method will
#'   effectively throw away the lower half from the returned `observedYN` and
#'   `observedNY` matrices, but, since they are transpose one of another, the
#'   full information is still available.
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean; when `TRUE` the function will return only
#'   packed dense symmetric matrices
#'
#' @returns `observedContingencyTables()` returns the observed contingency
#'   tables as named `list` with elements:
#'   * `"observedNN"`
#'   * `"observedNY"`
#'   * `"observedYN"`
#'   * `"observedYY"`
#'
#' @importFrom Matrix t
#' @importClassesFrom Matrix symmetricMatrix
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
observedContingencyTables <- function(objCOTAN,
                                      actOnCells = FALSE,
                                      asDspMatrices = FALSE) {
  numGenes <- getNumGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  c(observedYY, observedY) %<-%
    observedContingencyTablesYY(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = FALSE)
  gc()

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    logThis("calculating NY..", logLevel = 3L, appendLF = FALSE)

    # Any/Yes vector [cycled] = Yes/Yes + No/Yes
    #observedNY <- observedY - observedYY

    logThis("YN..", logLevel = 3L, appendLF = FALSE)
    observedYN <- t(observedY - observedYY)

    logThis("NN..", logLevel = 3L, appendLF = FALSE)
    observedNN <- numGenes - observedY - observedYN

    observedYN <- as(observedYN, "denseMatrix")
  } else {
    # dimension n x n (n number of genes)
    logThis("calculating YN..", logLevel = 3L, appendLF = FALSE)

    # Yes/Any vector [cycled] = Yes/Yes + Yes/No
    #observedYN <- observedY - observedYY

    logThis("NY..", logLevel = 3L, appendLF = FALSE)
    observedNY <- t(observedY - observedYY)

    logThis("NN..", logLevel = 3L, appendLF = FALSE)
    observedNN <- numCells - observedY - observedNY

    observedNY <- as(observedNY, "denseMatrix")
  }
  gc()
  logThis("t()..", logLevel = 3L, appendLF = FALSE)

  observedNN <- forceSymmetric(as(observedNN, "denseMatrix"))

  if (isTRUE(asDspMatrices)) {
    observedNN <- pack(observedNN)
    observedYY <- pack(observedYY)
    # these operation drops the lower triangle values
    # but the other matrix contains them anyway
    if (isTRUE(actOnCells)) {
      observedNY <- pack(forceSymmetric(t(observedYN)))
      observedYN <- pack(forceSymmetric(  observedYN) )
    } else {
      observedYN <- pack(forceSymmetric(t(observedNY)))
      observedNY <- pack(forceSymmetric(  observedNY) )
    }
  } else {
    if (isTRUE(actOnCells)) {
      observedNY <- t(observedYN)
    } else {
      observedYN <- t(observedNY)
    }
  }
  logThis(" done", logLevel = 3L)

  return(list("observedNN" = observedNN,
              "observedNY" = observedNY,
              "observedYN" = observedYN,
              "observedYY" = observedYY))
}


#' @details `observedPartialContingencyTables()` calculates the observed
#'   contingency tables.
#'
#' @param objCOTAN a `COTAN` object
#' @param columnsSubset a sub-set of the columns of the matrices that will be
#'   returned
#' @param zeroOne the raw count matrix projected to `0` or `1`. If not given the
#'   appropriate one will be calculated on the fly
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#'
#' @returns `observedPartialContingencyTables()` returns the observed
#'   contingency tables, restricted to the selected columns, as named `list`
#'   with elements:
#'   * `"observedNN"`
#'   * `"observedNY"`
#'   * `"observedYN"`
#'   * `"observedYY"`
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
observedPartialContingencyTables <-
  function(objCOTAN,
           columnsSubset,
           zeroOne = NULL,
           actOnCells = FALSE) {
  numGenes <- getNumGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  if (is.character(columnsSubset)) {
    if (isTRUE(actOnCells)) {
      allColumns <- getCells(objCOTAN)
    } else {
      allColumns <- getGenes(objCOTAN)
    }
    columnsSubset <- which(allColumns %in% columnsSubset)
  } else {
    columnsSubset <- sort(columnsSubset)
  }

  c(observedYY, observedY) %<-%
    observedPartialContingencyTablesYY(objCOTAN,
                                       columnsSubset = columnsSubset,
                                       zeroOne = zeroOne,
                                       actOnCells = actOnCells)
  gc()

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    logThis("calculating partial NY..", logLevel = 3L, appendLF = FALSE)

    # Any/Yes vector [cycled] = Yes/Yes + No/Yes
    observedNY <- observedY - observedYY

    logThis("YN..", logLevel = 3L, appendLF = FALSE)
    observedYN <- t(observedY[columnsSubset] - t(observedYY))

    logThis("NN..", logLevel = 3L, appendLF = FALSE)
    observedNN <- numGenes - observedY - observedYN
  } else {
    # dimension n x n (n number of genes)
    logThis("calculating partial YN..", logLevel = 3L, appendLF = FALSE)

    # Yes/Any vector [cycled] = Yes/Yes + Yes/No
    observedYN <- observedY - observedYY

    logThis("NY..", logLevel = 3L, appendLF = FALSE)
    observedNY <- t(observedY[columnsSubset] - t(observedYY))

    logThis("NN..", logLevel = 3L, appendLF = FALSE)
    observedNN <- numCells - observedY - observedNY
  }
  logThis(" done", logLevel = 3L)

  gc()

  return(list("observedNN" = observedNN,
              "observedNY" = observedNY,
              "observedYN" = observedYN,
              "observedYY" = observedYY))
}


#' @details `expectedContingencyTablesNN()` calculates the expected *No/No*
#'   field of the contingency table
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean; when `TRUE` the function will return only
#'   packed dense symmetric matrices
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use `Rfast`
#'   parallel algorithms that on the flip side use more memory
#'
#' @returns `expectedContingencyTablesNN()` returns a `list` with:
#'   * `expectedNN` the *No/No* expected contingency table as `matrix`
#'   * `expectedN`  the *No* expected `vector`
#'
#' @importFrom Matrix t
#' @importFrom Matrix colSums
#'
#' @importClassesFrom Matrix dgeMatrix
#' @importClassesFrom Matrix dsyMatrix
#'
#' @importFrom Rfast Crossprod
#' @importFrom Rfast Tcrossprod
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
expectedContingencyTablesNN <- function(objCOTAN,
                                        actOnCells = FALSE,
                                        asDspMatrices = FALSE,
                                        optimizeForSpeed = TRUE) {
  logThis("calculating NN..", logLevel = 3L, appendLF = FALSE)

  # estimate Probabilities of 0 with internal function funProbZero
  probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))
  gc()

  assert_that(!anyNA(probZero),
              msg = paste0("Error: some NA in matrix of probability",
                           " of zero UMI counts"))


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

  logThis(" done", logLevel = 3L)

  return(list("expectedNN" = expectedNN, "expectedN" = expectedN))
}


#' @details `expectedPartialContingencyTablesNN()` calculates the expected
#'   *No/No* field of the contingency table
#'
#' @param objCOTAN a `COTAN` object
#' @param columnsSubset a sub-set of the columns of the matrices that will be
#'   returned
#' @param probZero is the expected **probability of zero** for each gene/cell
#'   pair. If not given the appropriate one will be calculated on the fly
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use `Rfast`
#'   parallel algorithms that on the flip side use more memory
#'
#' @returns `expectedPartialContingencyTablesNN()` returns a `list` with:
#'   * `expectedNN` the *No/No* expected contingency table as `matrix`,
#'     restricted to the selected columns, as named `list` with elements
#'   * `expectedN`  the full *No* expected `vector`
#'
#' @importFrom Rfast Crossprod
#' @importFrom Rfast Tcrossprod
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
expectedPartialContingencyTablesNN <-
  function(objCOTAN,
           columnsSubset,
           probZero = NULL,
           actOnCells = FALSE,
           optimizeForSpeed = TRUE) {
  logThis("calculating partial NN..", logLevel = 3L, appendLF = FALSE)

  if (is_empty(probZero)) {
    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(getDispersion(objCOTAN), calculateMu(objCOTAN))
  }
  assert_that(identical(dim(probZero), dim(getRawData(objCOTAN))))
  assert_that(!anyNA(probZero),
              msg = paste0("Error: some NA in matrix of probability",
                           " of zero UMI counts"))
  gc()


  if (is.character(columnsSubset)) {
    if (isTRUE(actOnCells)) {
      allColumns <- getCells(objCOTAN)
    } else {
      allColumns <- getGenes(objCOTAN)
    }
    columnsSubset <- which(allColumns %in% columnsSubset)
  } else {
    columnsSubset <- sort(columnsSubset)
  }

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    if (isTRUE(optimizeForSpeed)) {
      # use Rfast functions
      expectedNN <- Crossprod(probZero,
                              probZero[, columnsSubset, drop = FALSE])
      expectedN  <- colsums(probZero, parallel = TRUE)
    } else {
      expectedNN <- crossprod(probZero,
                              probZero[, columnsSubset, drop = FALSE])
      expectedN  <- colSums(probZero)
    }
    rownames(expectedNN) <- getCells(objCOTAN)
    colnames(expectedNN) <- getCells(objCOTAN)[columnsSubset]
  } else {
    # dimension n x n (n number of genes)
    if (isTRUE(optimizeForSpeed)) {
      # use Rfast functions
      expectedNN <- Tcrossprod(probZero,
                               probZero[columnsSubset, , drop = FALSE])
      expectedN  <- rowsums(probZero, parallel = TRUE)
    } else {
      expectedNN <- tcrossprod(probZero,
                               probZero[columnsSubset, , drop = FALSE])
      expectedN  <- rowSums(probZero)
    }
    rownames(expectedNN) <- getGenes(objCOTAN)
    colnames(expectedNN) <- getGenes(objCOTAN)[columnsSubset]
  }
  rm(probZero)

  logThis(" done", logLevel = 3L)

  return(list("expectedNN" = expectedNN, "expectedN" = expectedN))
}


#' @details `expectedContingencyTables()` calculates the expected values of
#'   contingency tables. When the parameter `asDspMatrices == TRUE`, the method
#'   will effectively throw away the lower half from the returned `expectedYN`
#'   and `expectedNY` matrices, but, since they are transpose one of another,
#'   the full information is still available.
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param asDspMatrices Boolean; when `TRUE` the function will return only
#'   packed dense symmetric matrices
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use `Rfast`
#'   parallel algorithms that on the flip side use more memory
#'
#' @return `expectedContingencyTables()` returns the expected contingency tables
#'   as named `list` with elements:
#'   * `"expectedNN"`
#'   * `"expectedNY"`
#'   * `"expectedYN"`
#'   * `"expectedYY"`
#'
#' @importFrom Rfast Crossprod
#' @importFrom Rfast Tcrossprod
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @importFrom Matrix t
#'
#' @importClassesFrom Matrix symmetricMatrix
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
expectedContingencyTables <- function(objCOTAN,
                                      actOnCells = FALSE,
                                      asDspMatrices = FALSE,
                                      optimizeForSpeed = TRUE) {
  numGenes <- getNumGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  c(expectedNN, expectedN) %<-%
    expectedContingencyTablesNN(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = isFALSE(optimizeForSpeed),
                                optimizeForSpeed = optimizeForSpeed)

  gc()

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    logThis("calculating YN..", logLevel = 3L, appendLF = FALSE)

    # Any/No vector [cycled] = No/No + Yes/No
    #expectedYN <- expectedN - expectedNN

    logThis("NY..", logLevel = 3L, appendLF = FALSE)
    expectedNY <- t(expectedN - expectedNN)

    logThis("YY..", logLevel = 3L, appendLF = FALSE)
    expectedYY <- numGenes - expectedN - expectedNY

    expectedNY <- as(expectedNY, "denseMatrix")
  } else {
    # dimension n x n (n number of genes)
    logThis("calculating NY..", logLevel = 3L, appendLF = FALSE)

    # No/Any vector [cycled] = No/No + No/Yes
    #expectedNY <- expectedN - expectedNN

    logThis("YN..", logLevel = 3L, appendLF = FALSE)
    expectedYN <- t(expectedN - expectedNN)

    logThis("YY..", logLevel = 3L, appendLF = FALSE)
    expectedYY <- numCells - expectedN - expectedYN

    expectedYN <- as(expectedYN, "denseMatrix")
  }
  gc()
  logThis("t()..", logLevel = 3L, appendLF = FALSE)

  expectedYY <- forceSymmetric(as(expectedYY, "denseMatrix"))

  if (isTRUE(asDspMatrices)) {
    expectedNN <- pack(expectedNN)
    expectedYY <- pack(expectedYY)
    # these operation drops the lower triangle values
    # but the other matrix contains them anyway
    if (isTRUE(actOnCells)) {
      expectedYN <- pack(forceSymmetric(t(expectedNY)))
      expectedNY <- pack(forceSymmetric(  expectedNY) )
    } else {
      expectedNY <- pack(forceSymmetric(t(expectedYN)))
      expectedYN <- pack(forceSymmetric(  expectedYN) )
    }
  } else {
    if (isTRUE(actOnCells)) {
      expectedYN <- t(expectedNY)
    } else {
      expectedNY <- t(expectedYN)
    }
  }
  logThis(" done", logLevel = 3L)

  return(list("expectedNN" = expectedNN,
              "expectedNY" = expectedNY,
              "expectedYN" = expectedYN,
              "expectedYY" = expectedYY))
}


#' @details `expectedPartialContingencyTables()` calculates the expected values
#'   of contingency tables, restricted to the specified column sub-set
#'
#' @param objCOTAN a `COTAN` object
#' @param columnsSubset a sub-set of the columns of the matrices that will be
#'   returned
#' @param probZero is the expected **probability of zero** for each gene/cell
#'   pair. If not given the appropriate one will be calculated on the fly
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use `Rfast`
#'   parallel algorithms that on the flip side use more memory
#'
#' @return `expectedPartialContingencyTables()` returns the expected contingency
#'   tables, restricted to the selected columns, as named `list` with elements:
#'   * `"expectedNN"`
#'   * `"expectedNY"`
#'   * `"expectedYN"`
#'   * `"expectedYY"`
#'
#' @importFrom Rfast Crossprod
#' @importFrom Rfast Tcrossprod
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
expectedPartialContingencyTables <-
  function(objCOTAN,
           columnsSubset,
           probZero = NULL,
           actOnCells = FALSE,
           optimizeForSpeed = TRUE) {
  numGenes <- getNumGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  if (is.character(columnsSubset)) {
    if (isTRUE(actOnCells)) {
      allColumns <- getCells(objCOTAN)
    } else {
      allColumns <- getGenes(objCOTAN)
    }
    columnsSubset <- which(allColumns %in% columnsSubset)
  } else {
    columnsSubset <- sort(columnsSubset)
  }

  c(expectedNN, expectedN) %<-%
    expectedPartialContingencyTablesNN(objCOTAN, columnsSubset, probZero,
                                       actOnCells = actOnCells,
                                       optimizeForSpeed = optimizeForSpeed)

  gc()

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    logThis("calculating partial YN..", logLevel = 3L, appendLF = FALSE)

    # Any/No vector [cycled] = No/No + Yes/No
    expectedYN <- expectedN - expectedNN

    logThis("NY..", logLevel = 3L, appendLF = FALSE)
    expectedNY <- t(expectedN[columnsSubset] - t(expectedNN))

    logThis("YY..", logLevel = 3L, appendLF = FALSE)
    expectedYY <- numGenes - expectedN - expectedNY
  } else {
    # dimension n x n (n number of genes)
    logThis("calculating partial NY..", logLevel = 3L, appendLF = FALSE)

    # No/Any vector [cycled] = No/No + No/Yes
    expectedNY <- expectedN - expectedNN

    logThis("YN..", logLevel = 3L, appendLF = FALSE)
    expectedYN <- t(expectedN[columnsSubset] - t(expectedNN))

    logThis("YY..", logLevel = 3L, appendLF = FALSE)
    expectedYY <- numCells - expectedN - expectedYN
  }
  gc()

  logThis(" done", logLevel = 3L)

  return(list("expectedNN" = expectedNN,
              "expectedNY" = expectedNY,
              "expectedYN" = expectedYN,
              "expectedYY" = expectedYY))
}


#' @details `contingencyTables()` returns the observed and expected contingency
#'   tables for a given pair of genes. The implementation runs the same
#'   algorithms used to calculate the full observed/expected contingency tables,
#'   but restricted to only the relevant genes and thus much faster and less
#'   memory intensive
#'
#' @param objCOTAN a `COTAN` object
#' @param g1 a gene
#' @param g2 another gene
#'
#' @return `contingencyTables()` returns a list containing the observed and
#'   expected contingency tables
#'
#' @importFrom Matrix rowSums
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
contingencyTables <- function(objCOTAN, g1, g2) {
  genes <- getGenes(objCOTAN)
  numCells <- getNumCells(objCOTAN)

  assert_that(c(g1) %in% genes, msg = "the first gene is unknown")
  assert_that(c(g2) %in% genes, msg = "the second gene is unknown")

  dimnames <- list(c(paste(g2, "yes", sep = "."), paste(g2, "no", sep = ".")),
                   c(paste(g1, "yes", sep = "."), paste(g1, "no", sep = ".")))

  # observed

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
                       ncol = 2L, nrow = 2L, dimnames = dimnames)

  # estimated

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
                       ncol = 2L, nrow = 2L, dimnames = dimnames)


  return(list("observed" = observedCT, "expected" = expectedCT))
}


# --------------- COEX and GDI -----------

#' @aliases calculateCoex
#'
#' @details `calculateCoex()` estimates and stores the `COEX` matrix in the
#'   `cellCoex` or `genesCoex` field depending on given `actOnCells` flag. It
#'   also calculates the percentage of *problematic* genes/cells pairs. A pair
#'   is *problematic* when one or more of the expected counts were significantly
#'   smaller than 1 (\eqn{< 0.5}). These small expected values signal that scant
#'   information is present for such a pair.
#'
#' @param objCOTAN a `COTAN` object
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use `Rfast`
#'   parallel algorithms that on the flip side use more memory
#'
#' @returns `calculateCoex()` returns the updated `COTAN` object
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
setMethod(
  "calculateCoex",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE,
           optimizeForSpeed = TRUE) {
    if (isTRUE(actOnCells)) {
      kind <- "cells'"
    } else {
      kind <- "genes'"
    }
    logThis(paste("Calculate", kind, "coex: START"), logLevel = 1L)

    logThis(paste("Retrieving expected", kind, "contingency table"),
            logLevel = 3L)

    # four estimators:
    c(expectedNN, expectedNY, expectedYN, expectedYY) %<-%
      expectedContingencyTables(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = TRUE,
                                optimizeForSpeed = optimizeForSpeed)

    gc()

    logThis(paste("Calculating", kind, "coex normalization factor"),
            logLevel = 3L)

    if (isTRUE(actOnCells)) {
      allNames <- getCells(objCOTAN)
      normFact <- 1.0 / sqrt(getNumGenes(objCOTAN)) # divided by sqrt(n)
    } else {
      allNames <- getGenes(objCOTAN)
      normFact <- 1.0 / sqrt(getNumCells(objCOTAN)) # divided by sqrt(m)
    }

    problematicPairsFraction <-
      sum(pmin(expectedYY@x, expectedYN@x,
               expectedNY@x, expectedNN@x) < 0.5) / length(expectedYY@x)

    logThis(paste0("Fraction of genes with very low",
                   " expected contingency tables: ",
                   problematicPairsFraction), logLevel = 3L)

    coex <- normFact * sqrt(1.0 / pmax(1.0, expectedYY@x) +
                            1.0 / pmax(1.0, expectedYN@x) +
                            1.0 / pmax(1.0, expectedNY@x) +
                            1.0 / pmax(1.0, expectedNN@x))

    rm(expectedYN, expectedNY, expectedNN)
    gc()

    logThis(paste("Retrieving observed", kind, "yes/yes contingency table"),
            logLevel = 3L)

    c(observedYY, .) %<-%
      observedContingencyTablesYY(objCOTAN,
                                  actOnCells = actOnCells,
                                  asDspMatrices = TRUE)

    gc()

    # coex estimation
    logThis(paste("Estimating", kind, "coex"), logLevel = 3L)

    coex <- coex * (observedYY@x - expectedYY@x)

    assert_that(2L * length(coex) == length(allNames) * (length(allNames) + 1L),
                msg = "Output coex@x has the wrong size")

    coex <- new("dspMatrix", Dim = dim(expectedYY), x = coex)
    rownames(coex) <- colnames(coex) <- allNames

    rm(observedYY, expectedYY)
    gc()

    if (actOnCells) {
      objCOTAN@cellsCoex <- coex
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["csync"]], TRUE)
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["cbad"]],
                                             problematicPairsFraction)
    } else {
      objCOTAN@genesCoex <- coex
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["gsync"]], TRUE)
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["gbad"]],
                                             problematicPairsFraction)
    }

    rm(coex)
    gc()

    logThis(paste("Calculate", kind, "coex: DONE"), logLevel = 1L)

    return(objCOTAN)
  }
)


#' @details `calculatePartialCoex()` estimates a sub-section of the `COEX`
#'   matrix in the `cellCoex` or `genesCoex` field depending on given
#'   `actOnCells` flag. It also calculates the percentage of *problematic*
#'   genes/cells pairs. A pair is *problematic* when one or more of the expected
#'   counts were significantly smaller than 1 (\eqn{< 0.5}). These small
#'   expected values signal that scant information is present for such a pair.
#'
#' @param objCOTAN a `COTAN` object
#' @param columnsSubset a sub-set of the columns of the matrices that will be
#'   returned
#' @param probZero is the expected **probability of zero** for each gene/cell
#'   pair. If not given the appropriate one will be calculated on the fly
#' @param actOnCells Boolean; when `TRUE` the function works for the cells,
#'   otherwise for the genes
#' @param optimizeForSpeed Boolean; when `TRUE` the function will use `Rfast`
#'   parallel algorithms that on the flip side use more memory
#'
#' @returns `calculatePartialCoex()` returns the asked section of the `COEX`
#'   matrix
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
calculatePartialCoex <- function(objCOTAN,
                                 columnsSubset,
                                 probZero = NULL,
                                 zeroOne = NULL,
                                 actOnCells = FALSE,
                                 optimizeForSpeed = TRUE) {
    if (isTRUE(actOnCells)) {
      kind <- "cells'"
    } else {
      kind <- "genes'"
    }
    logThis(paste("Calculate", kind, "partial coex: START"), logLevel = 1L)

    logThis(paste("Retrieving expected", kind, "partial contingency table"),
            logLevel = 3L)

    # four estimators:
    c(expectedNN, expectedNY, expectedYN, expectedYY) %<-%
      expectedPartialContingencyTables(objCOTAN, columnsSubset,
                                       probZero = probZero,
                                       actOnCells = actOnCells,
                                       optimizeForSpeed = optimizeForSpeed)

    gc()

    logThis(paste("Calculating", kind, "partial coex normalization factor"),
            logLevel = 3L)

    if (isTRUE(actOnCells)) {
      normFact <- 1.0 / sqrt(getNumGenes(objCOTAN)) # divided by sqrt(n)
    } else {
      normFact <- 1.0 / sqrt(getNumCells(objCOTAN)) # divided by sqrt(m)
    }

    problematicPairsFraction <-
      sum(pmin(expectedYY, expectedYN,
               expectedNY, expectedNN) < 0.5) / length(expectedYY)

    logThis(paste0("Fraction of genes with very low",
                   " expected contingency tables: ",
                   problematicPairsFraction), logLevel = 3L)

    coex <- normFact * sqrt(1.0 / pmax(1.0, expectedYY) +
                            1.0 / pmax(1.0, expectedYN) +
                            1.0 / pmax(1.0, expectedNY) +
                            1.0 / pmax(1.0, expectedNN))

    rm(expectedYN, expectedNY, expectedNN)
    gc()

    logThis(paste("Retrieving observed", kind,
                  "yes/yes partial contingency table"),
            logLevel = 3L)

    c(observedYY, .) %<-%
      observedPartialContingencyTablesYY(objCOTAN,
                                         columnsSubset = columnsSubset,
                                         zeroOne = zeroOne,
                                         actOnCells = actOnCells)

    gc()

    # coex estimation
    logThis(paste("Estimating", kind, "partial coex"), logLevel = 3L)

    coex <- coex * (observedYY - expectedYY)

    rm(observedYY, expectedYY)
    gc()

    logThis(paste("Calculate", kind, "partial coex: DONE"), logLevel = 1L)

    return(coex)
  }


#' @details `calculateS()` calculates the statistics **S** for genes contingency
#'   tables. It always has the diagonal set to zero.
#'
#' @param objCOTAN a `COTAN` object
#' @param geneSubsetCol an array of genes. It will be put in columns. If left
#'   empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows. If left empty
#'   the function will do it genome-wide.
#'
#' @returns `calculateS()` returns the `S` matrix
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @rdname CalculatingCOEX
#'
calculateS <- function(objCOTAN, geneSubsetCol = vector(mode = "character"),
                       geneSubsetRow = vector(mode = "character")) {
  geneSubsetCol <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetCol)
  geneSubsetRow <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetRow)

  logThis("Calculating S: START", logLevel = 3L)

  assert_that(isCoexAvailable(objCOTAN), msg = "Coex is missing")

  coex <- getGenesCoex(objCOTAN, zeroDiagonal = TRUE)
  if (!identical(geneSubsetRow, getGenes(objCOTAN)) ||
      !identical(geneSubsetCol, getGenes(objCOTAN))) {
    coex <- coex[geneSubsetRow, geneSubsetCol, drop = FALSE]
  }

  assert_that(!is_empty(coex), msg = "Genes subset result in empty matrix")

  S <- coex^2L * getNumCells(objCOTAN)

  logThis("Calculating S: DONE", logLevel = 3L)

  return(S)
}


#' @details `calculateG()` calculates the statistics *G-test* for genes
#'   contingency tables. It always has the diagonal set to zero. It is
#'   proportional to the genes' presence mutual information.
#'
#' @param objCOTAN a `COTAN` object
#' @param geneSubsetCol an array of genes. It will be put in columns. If left
#'   empty the function will do it genome-wide.
#' @param geneSubsetRow an array of genes. It will be put in rows. If left empty
#'   the function will do it genome-wide.
#'
#' @returns `calculateG()` returns the G matrix
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
calculateG <- function(objCOTAN, geneSubsetCol = vector(mode = "character"),
                       geneSubsetRow = vector(mode = "character")) {
  geneSubsetCol <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetCol)
  geneSubsetRow <- handleNamesSubsets(getGenes(objCOTAN), geneSubsetRow)

  logThis("Calculating G: START", logLevel = 2L)

  c(observedNN, observedNY, observedYN, observedYY) %<-%
    observedContingencyTables(objCOTAN, actOnCells = FALSE)

  c(expectedNN, expectedNY, expectedYN, expectedYY) %<-%
    expectedContingencyTables(objCOTAN, actOnCells = FALSE)

  logThis("Estimating G", logLevel = 3L)

  tNN <- observedNN * log(observedNN / pmax(1.0, expectedNN))
  tNN[observedNN == 0L] <- 0.0
  diag(tNN) <- 0.0
  tNN <- tNN[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedNN, expectedNN)

  tNY <- observedNY * log(observedNY / pmax(1.0, expectedNY))
  tNY[observedNY == 0L] <- 0.0
  diag(tNY) <- 0.0
  tNY <- tNY[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedNY, expectedNY)

  tYN <- observedYN * log(observedYN / pmax(1.0, expectedYN))
  tYN[observedYN == 0L] <- 0.0
  diag(tYN) <- 0.0
  tYN <- tYN[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedYN, expectedYN)

  tYY <- observedYY * log(observedYY / pmax(1.0, expectedYY))
  tYY[observedYY == 0L] <- 0.0
  diag(tYY) <- 0.0
  tYY <- tYY[geneSubsetRow, geneSubsetCol, drop = FALSE]
  rm(observedYY, expectedYY)
  gc()

  G <- 2.0 * (tNN + tNY + tYN + tYY)

  if (identical(geneSubsetRow, getGenes(objCOTAN)) &&
      identical(geneSubsetCol, getGenes(objCOTAN))) {
    G <- pack(forceSymmetric(G))
  }

  rm(tNN, tNY, tYN, tYY)
  gc()

  logThis("Calculating G: DONE", logLevel = 2L)

  return(G)
}
