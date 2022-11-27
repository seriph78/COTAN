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
    lambda <- rowMeans(getRawData(objCOTAN), dims = 1)

    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes, lambda,
                                        "lambda", getGenes(objCOTAN))

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
    nu <- colMeans(getRawData(objCOTAN), dims = 1)
    nu <- nu / mean(nu)

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, nu,
                                        "nu", getCells(objCOTAN))

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
  function(objCOTAN, step = 512, threshold = 0.001,
           maxIterations = 1000, cores = 1) {
    print("Estimate dispersion: START")

    if (Sys.info()['sysname'] == "Windows" && cores != 1) {
      warning(paste0("On windows the numebr of cores used will be 1!",
                     " Multicore is not supported."))
      cores <- 1
    }

    genes <- getGenes(objCOTAN)
    zeroOneMatrix <- getZeroOneProj(objCOTAN)
    muEstimator <- calculateMu(objCOTAN)

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

      if (cores != 1) {
        res  <- parallel::mclapply(
                  spGenes[pBegin:pEnd],
                  parallelDispersionBisection,
                  zeroOneMatrix = zeroOneMatrix,
                  muEstimator = muEstimator,
                  threshold = threshold,
                  maxIterations = maxIterations,
                  mc.cores = cores)
      }
      else {
        res  <- lapply(
                  spGenes[pBegin:pEnd],
                  parallelDispersionBisection,
                  zeroOneMatrix = zeroOneMatrix,
                  muEstimator = muEstimator,
                  threshold = threshold,
                  maxIterations = maxIterations)
      }

      dispList <- append(dispList, res)
      rm(res)

      pBegin <- pEnd + 1
    }
    print("Estimate dispersion: DONE")

    gc()

    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes,
                                        unlist(dispList, recursive = TRUE, use.names = FALSE),
                                        "dispersion", genes)

    goodPos <- is.finite(getDispersion(objCOTAN))
    print(paste("dispersion",
                "| min:", min(getDispersion(objCOTAN)[goodPos]),
                "| max:", max(getDispersion(objCOTAN)[goodPos]),
                "| % negative:", (sum(getDispersion(objCOTAN) < 0) * 100 /
                                  getNumGenes(objCOTAN)) ))

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
  function(objCOTAN, step = 512, threshold = 0.001,
           maxIterations = 1000, cores = 1) {
    print("Estimate nu: START")

    if (Sys.info()['sysname'] == "Windows" && cores != 1) {
      warning(paste0("On windows the numebr of cores used will be 1!",
                     " Multicore is not supported."))
      cores <- 1
    }

    # parameters estimation
    if (is_empty(getNu(objCOTAN))) {
      stop("nu must not be empty, estimate it")
    }

    if (is_empty(getDispersion(objCOTAN))) {
      stop("dispersion must not be empty, estimate it")
    }

    # only genes not in housekeeping are used (to align to dispersion vector)
    cells <- getCells(objCOTAN)
    zeroOneMatrix <- getZeroOneProj(objCOTAN)
    lambda <- getLambda(objCOTAN)
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

      if (cores != 1) {
        res  <- parallel::mclapply(
                  spCells[pBegin:pEnd],
                  parallelNuBisection,
                  zeroOneMatrix = zeroOneMatrix,
                  lambda = lambda,
                  dispersion = dispersion,
                  initialGuess = nu,
                  threshold = threshold,
                  maxIterations = maxIterations,
                  mc.cores = cores)
      }
      else {
        res  <- lapply(
                  spCells[pBegin:pEnd],
                  parallelNuBisection,
                  zeroOneMatrix = zeroOneMatrix,
                  lambda = lambda,
                  dispersion = dispersion,
                  initialGuess = nu,
                  threshold = threshold,
                  maxIterations = maxIterations)
      }

      nuList <- append(nuList, res)
      rm(res)

      pBegin <- pEnd + 1
    }

    gc()
    print("Estimate nu: DONE")

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells,
                                        unlist(nuList, recursive = TRUE, use.names = FALSE),
                                        "nu", cells)

    # TODO: remove once code has been tested
    nuChange <- abs(nu - getNu(objCOTAN))
    print(paste("nu change",
                "| max:",     max(nuChange),
                "| median: ", median(nuChange),
                "| mean: ",   mean(nuChange) ))

    return(objCOTAN)
  }
)


#'estimateDispersionNuBisection
#'
#' Estimates the 'dispersion' and 'nu' field of a COTAN object by running
#' sequentially a bisection for each parameter. This allows to enforce
#' \code{mean('nu') == 1} assumption
#' Assumes \code{\link{estimateNuLinear}} being run already
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
#' @export
#'
#' @rdname estimateDispersionNuBisection
#'
setMethod(
  "estimateDispersionNuBisection",
  "COTAN",
  function(objCOTAN, step = 512, threshold = 0.001,
           maxIterations = 1000, cores = 1) {
    print("Estimate dispersion/nu: START")

    if (is_empty(getNu(objCOTAN))) {
      stop("Nu vector needs to be initialised as initial guess")
    }

    iter <- 1
    repeat {
      objCOTAN <- estimateDispersion(objCOTAN, step = step,
                                     threshold = threshold,
                                     maxIterations = maxIterations,
                                     cores = cores)

      objCOTAN <- estimateNuBisection(objCOTAN, step = step,
                                      threshold = threshold,
                                      maxIterations = maxIterations,
                                      cores = cores)

      meanNu <- mean(getNu(objCOTAN))

      print(paste("Nu mean:", meanNu))

      objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells,
                                          getNu(objCOTAN) / meanNu,
                                          "nu")

      if (abs(meanNu-1) < threshold / 100) {
        break
      }

      if (iter >= maxIterations / 100) {
        stop("Max number of outer iterations reached while finding the solution")
      }
    }

    print("Estimate dispersion/nu: DONE")

    return(objCOTAN)
  }
)
