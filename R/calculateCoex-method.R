# ----------- Contingency tables --------------

#'
#' @title Calculating the `COEX` matrix
#'
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
#'   expected contingency tables were enforced to match the marginal sums of the
#'   observed ones.
#'
#' @description The new implementation of the function relies on the `torch`
#'   package. This implies that it is potentially able to use the system `GPU`
#'   to run the heavy duty calculations required by this method. However
#'   installing the `torch` package on a system can be *finicky*, so we
#'   tentatively provide a short help page [Installing_torch] hoping that it
#'   will help...
#'
#' @seealso [ParametersEstimations] for more details.
#'
#' @seealso [Installing_torch] about the `torch` package
#'
#' @name CalculatingCOEX
#'
#' @rdname CalculatingCOEX
#'
#' @note The sum of the matrices returned by the function
#'   `observedContingencyTables()` and `expectedContingencyTables()` will have
#'   the same value on all elements. This value is the number of genes/cells
#'   depending on the parameter `actOnCells` being `TRUE/FALSE`.
#'
NULL

#' @details `calculateLikelihoodOfObserved()` gives for each cell and each gene
#'   the likelihood of the observed zero/one data
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `calculateLikelihoodOfObserved()` returns a `data.frame` with the
#'   likelihood of the observed zero/one
#'
#' @export
#'
#' @examples
#' lh <- calculateLikelihoodOfObserved(objCOTAN)
#'
#' @rdname CalculatingCOEX
#'
calculateLikelihoodOfObserved <- function(objCOTAN) {
  zeroOne <- getZeroOneProj(objCOTAN)

  probZero <- getProbabilityOfZero(objCOTAN)

  # estimate the likelihood of observed result
  return((1.0 - zeroOne) * probZero + zeroOne * (1.0 - probZero))
}


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
#' @param optimizeForSpeed Boolean; deprecated: always TRUE
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
  probZero <- getProbabilityOfZero(objCOTAN)
  gc()

  assert_that(!anyNA(probZero),
              msg = paste0("Error: some NA in matrix of probability",
                           " of zero UMI counts"))

  if (isTRUE(actOnCells)) {
    # dimension m x m (m number of cells)
    expectedNN <- Crossprod(probZero, probZero)
    expectedN  <- colsums(probZero, parallel = TRUE)
    rownames(expectedNN) <- colnames(expectedNN) <- getCells(objCOTAN)
  } else {
    # dimension n x n (n number of genes)
    expectedNN <- Tcrossprod(probZero, probZero)
    expectedN  <- rowsums(probZero, parallel = TRUE)
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
#' @param optimizeForSpeed Boolean; deprecated: always TRUE
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
    probZero <- getProbabilityOfZero(objCOTAN)
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
    expectedNN <- Crossprod(probZero,
                            probZero[, columnsSubset, drop = FALSE])
    expectedN  <- colsums(probZero, parallel = TRUE)
    rownames(expectedNN) <- getCells(objCOTAN)
    colnames(expectedNN) <- getCells(objCOTAN)[columnsSubset]
  } else {
    # dimension n x n (n number of genes)
    expectedNN <- Tcrossprod(probZero,
                             probZero[columnsSubset, , drop = FALSE])
    expectedN  <- rowsums(probZero, parallel = TRUE)
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
#' @param optimizeForSpeed Boolean; deprecated: always TRUE
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
                                asDspMatrices = FALSE,
                                optimizeForSpeed = TRUE)

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
#' @param optimizeForSpeed Boolean; deprecated: always TRUE
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
                                       optimizeForSpeed = TRUE)

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




# --------------- COEX and related matrices -----------

# legacy code for cases when the torch library is not available
calculateCoex_Legacy <- function(objCOTAN, actOnCells, returnPPFract) {
  startTime <- Sys.time()

  if (isTRUE(actOnCells)) {
    kind <- "cells'"
  } else {
    kind <- "genes'"
  }
  logThis(paste("Calculate", kind, "coex (legacy): START"), logLevel = 1L)

  logThis(paste("Retrieving expected", kind, "contingency table"),
          logLevel = 3L)

  # four estimators:
  c(expectedNN, expectedNY, expectedYN, expectedYY) %<-%
    expectedContingencyTables(objCOTAN,
                              actOnCells = actOnCells,
                              asDspMatrices = TRUE,
                              optimizeForSpeed = TRUE)

  expectedTime <- Sys.time()
  logThis(paste("Expected", kind, "contingency table elapsed time:",
                difftime(expectedTime, startTime, units = "secs")),
          logLevel = 3L)

  logThis(paste("Calculating", kind, "coex normalization factor"),
          logLevel = 3L)

  if (isTRUE(actOnCells)) {
    allNames <- getCells(objCOTAN)
    normFact <- 1.0 / sqrt(getNumGenes(objCOTAN)) # divided by sqrt(n)
  } else {
    allNames <- getGenes(objCOTAN)
    normFact <- 1.0 / sqrt(getNumCells(objCOTAN)) # divided by sqrt(m)
  }

  thresholdForPP <- 0.5
  problematicPairsFraction <- NA
  if (isTRUE(returnPPFract)) {
    problematicPairsFraction <-
      sum(pmin(expectedYY@x, expectedYN@x, expectedNY@x, expectedNN@x) <
            thresholdForPP) / length(expectedYY@x)

    logThis(paste("Fraction of", kind, "with very low",
                  "expected contingency tables:",
                  problematicPairsFraction), logLevel = 3L)
  }

  coex <- normFact * sqrt(1.0 / pmax(1.0, expectedYY@x) +
                          1.0 / pmax(1.0, expectedYN@x) +
                          1.0 / pmax(1.0, expectedNY@x) +
                          1.0 / pmax(1.0, expectedNN@x))

  rm(expectedYN, expectedNY, expectedNN)

  sqrtTime <- Sys.time()
  logThis(paste("Calculate", kind, "normalization factor elapsed time:",
                difftime(sqrtTime, expectedTime, units = "secs")),
          logLevel = 3L)

  logThis(paste("Retrieving observed", kind, "yes/yes contingency table"),
          logLevel = 3L)

  c(observedYY, .) %<-%
    observedContingencyTablesYY(objCOTAN,
                                actOnCells = actOnCells,
                                asDspMatrices = TRUE)

  observedTime <- Sys.time()
  logThis(paste("Observed", kind, "contingency table elapsed time:",
                difftime(observedTime, sqrtTime, units = "secs")),
          logLevel = 3L)

  # coex estimation
  logThis(paste("Estimating", kind, "coex"), logLevel = 3L)

  coex <- coex * (observedYY@x - expectedYY@x)

  assert_that(2L * length(coex) == length(allNames) * (length(allNames) + 1L),
              msg = "Output coex@x has the wrong size")

  coex <- new("dspMatrix", Dim = dim(expectedYY), x = coex)
  rownames(coex) <- colnames(coex) <- allNames

  gc()

  endTime <- Sys.time()
  logThis(paste("Calculate", kind, "coex elapsed time:",
                difftime(endTime, observedTime, units = "secs")),
          logLevel = 3L)

  logThis(paste("Total calculations elapsed time:",
                difftime(endTime, startTime, units = "secs")),
          logLevel = 2L)

  logThis(paste("Calculate", kind, "coex (legacy): DONE"), logLevel = 1L)

  return(list("coex" = coex, "ppf" = problematicPairsFraction))
}


# torch based code
calculateCoex_Torch <- function(objCOTAN, returnPPFract, deviceStr) {
  logThis(paste0("Calculate genes coex (torch) on device ", deviceStr,
                 ": START"), logLevel = 1L)

  startTime <- Sys.time()

  gc()
  torch::cuda_empty_cache()

  # TODO: 16 bits are OK here?
  if (deviceStr == "cpu") {
    dtypeForCalc <- torch::torch_float32()
    halfDtypeForCalc <- torch::torch_float32()
  } else {
    dtypeForCalc <- torch::torch_float32()
    halfDtypeForCalc <- torch::torch_float16()
  }

  device <- torch::torch_device(deviceStr)

  probOne <- function(nu, lambda, a) {
    zero     <- torch::torch_tensor( 0.0, device = device,
                                    dtype = torch::torch_float64())
    minusOne <- torch::torch_tensor(-1.0, device = device,
                                    dtype = torch::torch_float64())

    # Calculate terms based on conditions
    term1 <- (a <= zero) *
      torch::torch_exp(torch::torch_ger(nu, lambda) *
                         (torch::torch_minimum(a, zero) + minusOne))

    term2 <- (a >  zero) *
      torch::torch_pow(one + torch::torch_ger(nu, lambda)
                       * torch::torch_maximum(a, zero), minusOne / a)

    invisible(term1$add_(term2)$add_(minusOne)$neg_())

    return(term1)
  }

  logThis("Retrieving expected genes contingency table", logLevel = 3L)

  one <- torch::torch_tensor(1.0, device = device,
                             dtype = dtypeForCalc)
  m   <- torch::torch_tensor(getNumCells(objCOTAN),
                             device = device, dtype = dtypeForCalc)

  expectedYY <- torch::torch_tensor(probOne(
    torch::torch_tensor(getNu(objCOTAN),
                        dtype = torch::torch_float64(), device = device),
    torch::torch_tensor(getLambda(objCOTAN),
                        dtype = torch::torch_float64(), device = device),
    torch::torch_tensor(getDispersion(objCOTAN),
                        dtype = torch::torch_float64(), device = device)),
    device = device, dtype = dtypeForCalc)
  torch::cuda_empty_cache()

  expectedY <- torch::torch_sum(expectedYY, 1L, dtype = dtypeForCalc)

  expectedYY <- torch::torch_mm(torch::torch_t(expectedYY), expectedYY)

  expectedTime <- Sys.time()
  logThis(paste("Expected genes contingency table elapsed time:",
                difftime(expectedTime, startTime, units = "secs")),
          logLevel = 3L)

  logThis("Calculating genes coex normalization factor", logLevel = 3L)

  coex <- torch::torch_maximum(expectedYY, one)$reciprocal_()

  thresholdForPP <- 0.5
  problematicPairs <- NULL
  if (isTRUE(returnPPFract)) {
    problematicPairs <-
      torch::torch_tensor(expectedYY < thresholdForPP,
                          dtype = torch::torch_bool(), device = device)
  }

  # expectedYN
  #TODO: use in-place torch_maximum_() as soon as it becomes available
  tmp <- expectedY$view(c(-1L, 1L)) - expectedYY
  if (isTRUE(returnPPFract)) {
    invisible(problematicPairs$bitwise_or_(tmp < thresholdForPP))
  }
  invisible(coex$add_(tmp$sub_(one)$relu_()$add_(one)$reciprocal_()))
  rm(tmp)
  torch::cuda_empty_cache()

  # expectedNY
  tmp <- expectedY$view(c(1L, -1L)) - expectedYY
  if (isTRUE(returnPPFract)) {
    invisible(problematicPairs$bitwise_or_(tmp < thresholdForPP))
  }
  invisible(coex$add_(tmp$sub_(one)$relu_()$add_(one)$reciprocal_()))
  rm(tmp)
  torch::cuda_empty_cache()

  # expectedNN
  tmp <- expectedYY - expectedY$view(c(-1L, 1L))
  invisible(expectedY$sub_(m))
  invisible(tmp$sub_(expectedY$view(c(1L, -1L))))
  if (isTRUE(returnPPFract)) {
    invisible(problematicPairs$bitwise_or_(tmp < thresholdForPP))
  }
  invisible(coex$add_(tmp$sub_(one)$relu_()$add_(one)$reciprocal_()))
  rm(tmp, expectedY)
  gc()
  torch::cuda_empty_cache()

  invisible(coex$divide_(m)$sqrt_())

  # count problematic pairs
  problematicPairsFraction <- NA
  if (isTRUE(returnPPFract)) {
    numDiagPP <- torch::torch_diag(problematicPairs)$sum()$item()
    problematicPairsFraction <-
      ((problematicPairs$sum()$item() + numDiagPP) /
         (problematicPairs$size(1L) * (problematicPairs$size(1L) + 1L)))

    logThis(paste("Fraction of genes with very low",
                  "expected contingency tables:",
                  problematicPairsFraction), logLevel = 3L)

    rm(problematicPairs)
    torch::cuda_empty_cache()
  }

  sqrtTime <- Sys.time()
  logThis(paste("Calculating genes coex normalization factor elapsed time:",
                difftime(sqrtTime, expectedTime, units = "secs")),
          logLevel = 3L)

  logThis("Retrieving observed genes yes/yes contingency table", logLevel = 3L)

  # observedYY
  observedYY <- torch::torch_tensor(
    as.matrix(getRawData(objCOTAN)),
    device = device, dtype = torch::torch_int16())
  observedYY <- torch::torch_tensor(observedYY != 0L, dtype = halfDtypeForCalc)

  observedYY <- torch::torch_mm(observedYY, torch::torch_t(observedYY))

  observedTime <- Sys.time()
  logThis(paste("Observed genes contingency table elapsed time:",
                difftime(observedTime, sqrtTime, units = "secs")),
          logLevel = 3L)

  invisible(expectedYY$subtract_(observedYY))
  invisible(coex$multiply_(expectedYY$neg()))

  rm(expectedYY, observedYY)
  torch::cuda_empty_cache()

  ret <- pack(forceSymmetric(as(as.matrix(coex$cpu()), "denseMatrix")))
  rownames(ret) <- colnames(ret) <- getGenes(objCOTAN)

  rm(coex, device)
  gc()
  torch::cuda_empty_cache()

  endTime <- Sys.time()
  logThis(paste("Calculate genes coex elapsed time:",
                difftime(endTime, observedTime, units = "secs")),
          logLevel = 3L)
  logThis(paste("Total calculations elapsed time:",
                difftime(endTime, startTime, units = "secs")),
          logLevel = 2L)

  logThis(paste0("Calculate genes coex (torch) on device ", deviceStr,
                 ": DONE"), logLevel = 1L)

  return(list("coex" = ret, "ppf" = problematicPairsFraction))
}


## ---- Calculate COEX ----

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
#' @param returnPPFract Boolean; when `TRUE` the function returns the fraction
#'   of genes/cells pairs for which the *expected contingency table* is smaller
#'   than \eqn{0.5}. Default is FALSE
#' @param optimizeForSpeed Boolean; when `TRUE` `COTAN` tries to use the `torch`
#'   library to run the matrix calculations. Otherwise, or when the library is
#'   not available will run the slower legacy code
#' @param deviceStr On the `torch` library enforces which device to use to run
#'   the calculations. Possible values are `"cpu"` to us the system *CPU*,
#'   `"cuda"` to use the system *GPUs* or something like `"cuda:0"` to restrict
#'   to a specific device
#'
#' @returns `calculateCoex()` returns the updated `COTAN` object
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @importFrom methods new
#'
#' @export
#'
#' @rdname CalculatingCOEX
#'
setMethod(
  "calculateCoex",
  "COTAN",
  function(objCOTAN, actOnCells = FALSE, returnPPFract = FALSE,
           optimizeForSpeed = TRUE, deviceStr = "cuda") {
    coex <- NULL
    problematicPairsFraction <- NA

    if (isTRUE(actOnCells)) {
      if (isTRUE(optimizeForSpeed)) {
        warning("The 'torch' package is not supported yet for cells' COEX",
                " Falling back to legacy code.")
      }
      c(coex, problematicPairsFraction) %<-%
        calculateCoex_Legacy(objCOTAN, actOnCells = TRUE,
                             returnPPFract = returnPPFract)

      objCOTAN@cellsCoex <- coex
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["csync"]], TRUE)
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["cbad"]],
                                             problematicPairsFraction)
    } else {
      c(useTorch, deviceStr) %<-% canUseTorch(optimizeForSpeed, deviceStr)

      if (useTorch) {
        c(coex, problematicPairsFraction) %<-%
          calculateCoex_Torch(objCOTAN, deviceStr = deviceStr,
                              returnPPFract = returnPPFract)
      } else {
        c(coex, problematicPairsFraction) %<-%
          calculateCoex_Legacy(objCOTAN, actOnCells = FALSE,
                               returnPPFract = returnPPFract)
      }

      objCOTAN@genesCoex <- coex
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["gsync"]], TRUE)
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["gbad"]],
                                             problematicPairsFraction)
    }

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
#' @param optimizeForSpeed Boolean; deprecated: always TRUE
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
                                       optimizeForSpeed = TRUE)

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
