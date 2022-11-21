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


#' estimateNormalisedData
#'
#' This function estimates and initializes the normalized count table.
#'
#' @param objCOTAN A COTAN object
#' @return the updated COTAN object
#'
#' @importFrom rlang is_empty
#' @importFrom Matrix t
#'
#' @export
#'
#' @rdname estimateNormalisedData
#'
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
#' @param objCOTAN a COTAN object
#' @return the updated COTAN object
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
#'
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

    return(objCOTAN)
  }
)


#' estimateDispersion
#'
#' This is the main function that estimates the dispersion vector
#' to store all the negative binomial dispersion factors.
#' It needs to be run after \code{\link{clean}}
#'
#' @param objCOTAN A COTAN object
#' @param step number of genes to solve in batch in a single core. Default is 256.
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param cores number of cores to use. Default is 1.
#'
#' @return the updated COTAN object
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
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
  function(objCOTAN, step = 256, threshold = 0.001,
           maxIterations = 1000, cores = 1) {
    print("Estimate dispersion: START")

    if (Sys.info()['sysname'] == "Windows" && cores != 1) {
      warning(paste0("On windows the numebr of cores used will be 1!",
                     " Multicore is not supported."))
      cores <- 1
    }

    genes <- getGenes(objCOTAN)
    zeroOneMatrix <- getZeroOneProj(objCOTAN)
    muEstimator <- estimateMu(objCOTAN)

    # exclude the effectively ubiquitous genes and saved in a separate file
    objCOTAN <- findHousekeepingGenes(objCOTAN)

    housekeepingGenes <- getHousekeepingGenes(objCOTAN)

    dispList <- list()

    spIdx <- parallel::splitIndices(length(genes), ceiling(length(genes) / step))

    spGenes = lapply(spIdx, function(x) genes[x])

    numSplits <- length(spGenes)
    splitStep <- max(16, cores * 2)

    pBegin <- 1
    while (pBegin <= numSplits) {
      pEnd <- min(pBegin + splitStep - 1, numSplits)

      print(paste0("Executing ", (pEnd - pBegin + 1), " genes batches from",
                   " [", spIdx[pBegin], "] to [", spIdx[pEnd], "]"))

      res <- parallel::mclapply(
              spGenes[pBegin:pEnd],
              parallelDispersionBisection,
              zeroOneMatrix = zeroOneMatrix,
              muEstimator = muEstimator,
              housekeepingGenes = housekeepingGenes,
              threshold = threshold,
              maxIterations = maxIterations,
              mc.cores = cores)

      dispList <- append(dispList, res)
      rm(res)

      pBegin <- pEnd + 1
    }

    gc()
    print("Estimate dispersion: DONE")

    objCOTAN@dispersion <- unlist(dispList, recursive = TRUE, use.names = FALSE)
    names(objCOTAN@dispersion) <- genes

    print(paste("dispersion",
                "| min:", min(getDispersion(objCOTAN), na.rm = TRUE),
                "| max:", max(getDispersion(objCOTAN), na.rm = TRUE),
                "| % negative:", (sum(getDispersion(objCOTAN) < 0, na.rm = TRUE) * 100 /
                                  length(getDispersion(objCOTAN))) ))

    # drop housekeeping values [as they are all NA]
    # TODO: keep all values and provide full and filterd accessors
    objCOTAN@dispersion <- objCOTAN@dispersion[!genes %in% housekeepingGenes]

    return(objCOTAN)
  }
)


#'estimateNuBisection
#'
#' Estimates the 'nu' field of a COTAN object by bisection.
#' Assumes dispersion being already calculated
#'
#' @param objCOTAN a COTAN object
#' @param step number of genes to solve in batch in a single core. Default is 256.
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param cores number of cores to use. Default is 1.
#'
#' @return the updated COTAN object
#'
#' @importFrom rlang is_empty
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @export
#'
#' @rdname estimateNuBisection
#'
setMethod(
  "estimateNuBisection",
  "COTAN",
  function(objCOTAN, step = 256, threshold = 0.001,
           maxIterations = 1000, cores = 1) {
    print("Estimate nu: START")

    if (Sys.info()['sysname'] == "Windows" && cores != 1) {
      warning(paste0("On windows the numebr of cores used will be 1!",
                     " Multicore is not supported."))
      cores <- 1
    }


    # parameters estimation
    if (is_empty(getLambda(objCOTAN))) {
      objCOTAN <- estimateLambdaLinear(objCOTAN)
    }

    if (is_empty(getHousekeepingGenes(objCOTAN))) {
      objCOTAN <- findHousekeepingGenes(objCOTAN)
    }

    if (is_empty(getNu(objCOTAN))) {
      # nu vector is empty: estimate it linearly as initial guesses
      objCOTAN <- estimateNuLinear(objCOTAN)
    }

    if (is_empty(getDispersion(objCOTAN))) {
      objCOTAN <- estimateDispersion(objCOTAN, step = step,
                                     threshold = threshold,
                                     maxIterations = maxIterations,
                                     cores = cores)
    }

    # only genes not in housekeeping are used (to align to dispersion vector)
    cells <- getCells(objCOTAN)
    zeroOneMatrix <- getZeroOneProj(objCOTAN)
    lambda <- getLambda(objCOTAN)[flagNotHousekeepingGenes(objCOTAN)]
    dispersion <- getDispersion(objCOTAN)
    nu <- getNu(objCOTAN)

    nuList <- list()

    spIdx <- parallel::splitIndices(length(cells), ceiling(length(cells) / step))

    spCells = lapply(spIdx, function(x) cells[x])

    numSplits <- length(spCells)
    splitStep <- max(16, cores * 2)

    pBegin <- 1
    while (pBegin <= numSplits) {
      pEnd <- min(pBegin + splitStep - 1, numSplits)

      print(paste0("Executing ", (pEnd - pBegin + 1), " cells batches from",
                   " [", spIdx[pBegin], "] to [", spIdx[pEnd], "]"))

      res <- parallel::mclapply(
        spCells[pBegin:pEnd],
        parallelNuBisection,
        zeroOneMatrix = zeroOneMatrix,
        lambda = lambda,
        dispersion = dispersion,
        initialGuess = nu,
        threshold = threshold,
        maxIterations = maxIterations,
        mc.cores = cores)

      nuList <- append(nuList, res)
      rm(res)

      pBegin <- pEnd + 1
    }

    gc()
    print("Estimate nu: DONE")

    objCOTAN@nu <- unlist(nuList, recursive = TRUE, use.names = FALSE)
    names(objCOTAN@nu) <- cells

    # TODO: remove once code has been tested
    nuChange <- abs(nu - getNu(objCOTAN))
    print(paste("nu change",
                "| max:",     max(nuChange),
                "| median: ", median(nuChange),
                "| mean: ",   mean(nuChange) ))

    return(objCOTAN)
  }
)
