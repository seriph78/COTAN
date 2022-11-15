# COTAN parameters' estimates methods

#' estimateLambdaLinear
#'
# linear estimator of lambda (genes' averages)
#' @param objCOTAN a COTAN object
#' @return the updated COTAN object
#'
#' @importFrom Matrix mean
#' @importFrom Matrix rowMeans
#' @export
#'
#' @rdname estimateLambdaLinear
setMethod(
  "estimateLambdaLinear",
  "COTAN",
  function(objCOTAN) {
    objCOTAN@lambda <- rowMeans(getRawData(objCOTAN), dims = 1, na.rm = TRUE)

    return(objCOTAN)
  }
)


#' estimateNuLinear
#'
#' linear estimator of nu (normalised cells' averages)
#' @param objCOTAN a COTAN object
#' @return the updated COTAN object
#'
#' @importFrom Matrix colMeans
#' @export
#'
#' @rdname estimateNuLinear
setMethod(
  "estimateNuLinear",
  "COTAN",
  function(objCOTAN) {
    # raw column averages divided by global_mean
    objCOTAN@nu <- colMeans(getRawData(objCOTAN), dims = 1, na.rm = TRUE)
    objCOTAN@nu <- objCOTAN@nu / mean(objCOTAN@nu)

    return(objCOTAN)
  }
)


#' estimateMu
#'
#' estimate vector mu
#' @param objCOTAN a COTAN object
#' @return the estimated Mu matrix
#'
#' @importFrom Matrix t
#' @export
#'
#' @rdname estimateMu
setMethod(
  "estimateMu",
  "COTAN",
  function(objCOTAN) {
    muEstimator <- getLambda(objCOTAN) %*% t(getNu(objCOTAN))

    rownames(muEstimator) <- getGenes(objCOTAN)

    return(muEstimator)
  }
)

#' expectedContingencyTables
#'
#' method for estimating the expected values of contingency tables
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells,
#' otherwise for the genes
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
  function(objCOTAN, cells) {
    # zeroOne matrix : formed by row data matrix changed to 0-1 matrix
    zeroOne <- getZeroOneProj(objCOTAN)

    mu <- estimateMu(objCOTAN)[!getGenes(objCOTAN) %in% getHousekeepingGenes(objCOTAN), ]

    # estimate Probabilities of 0 with internal function funProbZero
    probZero <- funProbZero(getDispersion(objCOTAN), mu)

    errMgs <- "Error: some NA in matrix of probability of zero UMI counts. "
    stopifnot(errMgs = !anyNA(probZero))

    if(cells) {
      # dimension m x m (m number of cells)
      expectedNN <- t(probZero) %*% probZero
      expectedNY <- t(1 - probZero) %*% probZero
      expectedYN <- t(expectedNY)
      expectedYY <- t(1 - probZero) %*% (1 - probZero)
    } else {
      # dimension n x n (n number of genes)
      expectedNN <- probZero %*% t(probZero)
      expectedNY <- probZero %*% t(1 - probZero)
      expectedYN <- t(expectedNY)
      expectedYY <- (1 - probZero) %*% t(1 - probZero)
    }

    out <- list( "expectedNN" = as.matrix(expectedNN),
                 "expectedNY" = as.matrix(expectedNY),
                 "expectedYN" = as.matrix(expectedYN),
                 "expectedYY" = as.matrix(expectedYY) )

    return(out)
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
#' @importClassesFrom Matrix symmetricMatrix
#' @export
#'
#' @rdname observedContingencyYY
setMethod(
  "observedContingencyYY",
  "COTAN",
  function(objCOTAN, cells) {
    zeroOne <- getZeroOneProj(objCOTAN)

    if(isTRUE(cells)){
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


#' estimateNormalisedData
#'
#' This function estimates and initializes the normalized count table.
#'
#' @param objCOTAN A COTAN object
#' @return the updated COTAN object
#'
#' @importFrom rlang is_empty
#' @importFrom Matrix t
#' @export
#'
#' @rdname estimateNormalisedData
setMethod(
  "estimateNormalisedData",
  "COTAN",
  function(objCOTAN) {
    if (is_empty(getRawData(objCOTAN))) {
      stop("empty raw")
    }

    if (is_empty(getNu(objCOTAN))) {
      stop("nu must not be empty, estimate it")
    }

    objCOTAN@rawNorm <- t(t(getRawData(objCOTAN)) * (1/(getNu(objCOTAN))))

    return(objCOTAN)
  }
)


#' runEstimatesLinear
#'
#' Internal function to estimate the cell efficiency
#' @param objCOTAN COTAN object
#' @return a list of object (dist_cells, pca_cells, objCOTAN)
#'
#' @importFrom Matrix t
#' @importFrom Matrix mean
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colMeans
#'
#' @importFrom irlba prcomp_irlba
#'
#' @importFrom stats dist
#'
#' @rdname runEstimatesLinear
setMethod(
  "runEstimatesLinear",
  "COTAN",
  function(objCOTAN) {

    print("Start estimation mu with linear method")
    print(paste0("Working on [", getNumGenes(objCOTAN), "] genes and [", getNumCells(objCOTAN), "] cells"))

    objCOTAN <- estimateLambdaLinear(objCOTAN)
    objCOTAN <- estimateNuLinear(objCOTAN)
    objCOTAN <- estimateNormalisedData(objCOTAN)

    # genesMeans <- getNu(objCOTAN)
    # genesRng <- round(getNumGenes(objCOTAN)     / 2, digits = 0)
    #           : round(getNumGenes(objCOTAN) * 3 / 4, digits = 0)
    # genesMax <- names(sort(genesMeans, decreasing = TRUE)[genesRng])

    gc()

    print("Start PCA")

    pca_cells <- irlba::prcomp_irlba(t(getNormalizedData(objCOTAN)), n = 5)[["x"]]
    rownames(pca_cells) <- getCells(objCOTAN)

    gc()

    ppp <- pca_cells
    ppp <- scale(ppp)
    dist_cells <- stats::dist(ppp, method = "euclidean") # mhalanobis

    rm(ppp)
    gc()

    pca_cells <- as.data.frame(pca_cells)

    print("PCA DONE")

    return( list(dist_cells, pca_cells, objCOTAN) )
  }
)
