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

    print("Linear estimations: START")
    print(paste0("Working on [", getNumGenes(objCOTAN), "] genes and [", getNumCells(objCOTAN), "] cells"))

    objCOTAN <- estimateLambdaLinear(objCOTAN)
    objCOTAN <- estimateNuLinear(objCOTAN)
    objCOTAN <- estimateNormalisedData(objCOTAN)

    # genesMeans <- getNu(objCOTAN)
    # genesRng <- round(getNumGenes(objCOTAN)     / 2, digits = 0)
    #           : round(getNumGenes(objCOTAN) * 3 / 4, digits = 0)
    # genesMax <- names(sort(genesMeans, decreasing = TRUE)[genesRng])

    print("Linear estimations: DONE")

    gc()

    print("PCA: START")

    pcaCells <- irlba::prcomp_irlba(t(getNormalizedData(objCOTAN)), n = 5)[["x"]]
    rownames(pcaCells) <- getCells(objCOTAN)

    gc()

    distCells <- stats::dist(scale(pcaCells), method = "euclidean") # mhalanobis

    gc()

    pcaCells <- as.data.frame(pcaCells)

    print("PCA: DONE")

    return( list(distCells, pcaCells, objCOTAN) )
  }
)


#' estimateDispersion
#'
#' This is the main function that estimates the dispersion vector
#' to store all the negative binomial dispersion factors.
#' It needs to be run after \code{\link{clean}}
#'
#' @param objCOTAN A COTAN object
#' @param cores number of cores to use. Default is 11.
#' @return the updated COTAN object
#'
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' ERCC.cotan <- clean(ERCC.cotan)
#' ERCC.cotan <- estimateDispersion(ERCC.cotan)
#'
#' @rdname estimateDispersion
#'
setMethod(
  "estimateDispersion",
  "COTAN",
  function(objCOTAN, cores = 1) {
    print("Estimate dispersion: START")

    if (Sys.info()['sysname'] == "Windows" && cores != 1) {
      warning(paste0("On windows the numebr of cores used will be 1!",
                     " Multicore is not supported."))
      cores <- 1
    }

    # exclude the effectively ubiquitous genes and saved in a separate file
    objCOTAN <- findHousekeepingGenes(objCOTAN)

    genesPos <- which(!getGenes(objCOTAN) %in% getHousekeepingGenes(objCOTAN))

    genes <- getGenes(objCOTAN)[genesPos]
    zeroOneMatrix <- getZeroOneProj(objCOTAN)[genesPos, ]
    muEst <- estimateMu(objCOTAN)[genesPos, ]

    # muEst <- as.matrix(muEst)
    # print("start a minimization")

    dispList <- list()

    numGenes <- length(genes)
    pBegin <- 1
    while (pBegin <= numGenes) {
      pEnd <- pBegin + 200
      if (pEnd >= numGenes) {
        print("Final genes' trance!")
        pEnd = numGenes
      }

      res <- parallel::mclapply(
        genes[pBegin:pEnd],
        dispersionBisection,
        zeroOneMatrix = zeroOneMatrix,
        muEstimator = muEst,
        mc.cores = cores)

      dispList <- append(dispList, res)
      rm(res)

      pBegin <- pEnd + 1
      if ((pBegin %% 10) == 0) {
        print(paste("Next gene:", genes[pBegin], "number", pBegin))
      }
    }

    gc()
    print("Estimate dispersion: DONE")

    dispDf <- dispList[[1]]
    for (i in 2:length(dispList)) {
      dispDf<- rbind(dispDf, dispList[[i]])
    }

    rm(dispList)
    gc()

    objCOTAN@dispersion <- dispDf[["dispersion"]]
    names(objCOTAN@dispersion) <- genes

    rm(dispDf)
    gc()

    print(paste("dispersion",
                "| min:", min(getDispersion(objCOTAN)),
                "| max:", max(getDispersion(objCOTAN)),
                "| % negative:", (sum(getDispersion(objCOTAN) <0) * 100 /
                                  length(getDispersion(objCOTAN))) ))

    return(objCOTAN)
  }
)

