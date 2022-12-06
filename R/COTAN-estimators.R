# COTAN parameters' estimates methods

#' estimateLambdaLinear
#'
#' @description Linear estimator of lambda (genes' counts averages)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns The updated `COTAN` object
#'
#' @importFrom Matrix mean
#' @importFrom Matrix rowMeans
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- estimateLambda(objCOTAN)
#' lambda <- getLambda(objCOTAN)
#'
#' @rdname estimateLambdaLinear
#'
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
#' @description linear estimator of nu (normalised cells' counts averages)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns The updated `COTAN` object
#'
#' @importFrom Matrix colMeans
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- estimateNuLinear(objCOTAN)
#' nu <- getNu(objCOTAN)
#'
#' @rdname estimateNuLinear
#'
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


#' estimateDispersionBisection
#'
#' @description Estimates the negative binomial dispersion factor for each gene.
#'
#' @details Determines the `dispersion` such that, for each gene, the
#'   probability of zero count matches the number of observed zeros. It assumes
#'   [es()] being already run.
#'
#' It needs to
#'   be run after [clean()]
#'
#' @param objCOTAN A `COTAN` object
#' @param step number of genes to solve in batch in a single core. Default is
#'   256.
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param cores number of cores to use. Default is 1.
#'
#' @returns The updated `COTAN` object
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#' dispersion <- getDispersion(objCOTAN)
#'
#' @rdname estimateDispersionBisection
#'
setMethod(
  "estimateDispersionBisection",
  "COTAN",
  function(objCOTAN, step = 512, threshold = 0.001,
           maxIterations = 1000, cores = 1) {
    logThis("Estimate dispersion: START", logLevel = 2)

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

      logThis(paste0("Executing ", (pEnd - pBegin + 1), " genes batches from",
                     " [", spIdx[pBegin], "] to [", spIdx[pEnd], "]"),
              logLevel = 3)

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
    logThis("Estimate dispersion: DONE", logLevel = 2)

    gc()

    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes,
                                        unlist(dispList, recursive = TRUE, use.names = FALSE),
                                        "dispersion", genes)

    goodPos <- is.finite(getDispersion(objCOTAN))
    logThis(paste("dispersion",
                  "| min:", min(getDispersion(objCOTAN)[goodPos]),
                  "| max:", max(getDispersion(objCOTAN)[goodPos]),
                  "| % negative:", (sum(getDispersion(objCOTAN) < 0) * 100 /
                                    getNumGenes(objCOTAN)) ),
            logLevel = 1)

    return(objCOTAN)
  }
)


#' estimateNuBisection
#'
#' @description Estimates the `nu` vector of a `COTAN` object by bisection.
#'
#' @details Determines the `nu` parameters such that, for each cell, the
#'   probability of zero count matches the number of observed zeros. It assumes
#'   [estimateDispersionBisection()] being already run.
#'
#' @param objCOTAN a `COTAN` object
#' @param step number of genes to solve in batch in a single core. Default is
#'   256.
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param cores number of cores to use. Default is 1.
#'
#' @return the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @export
#'
#' @examples
#' objCOTAN <- COTAN(raw = data("ERCC.cotan"))
#' objCOTAN <- clean(objCOTAN, calcExtraData = FALSE)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#' objCOTAN <- estimateNuBisection(objCOTAN, cores = 12)
#' nu <- getNu(objCOTAN)
#'
#' @rdname estimateNuBisection
#'
setMethod(
  "estimateNuBisection",
  "COTAN",
  function(objCOTAN, step = 512, threshold = 0.001,
           maxIterations = 1000, cores = 1) {
    logThis("Estimate nu: START", logLevel = 2)

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

      logThis(paste0("Executing ", (pEnd - pBegin + 1), " cells batches from",
                     " [", spIdx[pBegin], "] to [", spIdx[pEnd], "]"),
              logLevel = 3)

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
    logThis("Estimate nu: DONE", logLevel = 2)

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells,
                                        unlist(nuList, recursive = TRUE, use.names = FALSE),
                                        "nu", cells)

    {
      nuChange <- abs(nu - getNu(objCOTAN))
      logThis(paste("nu change",
                    "| max:",     max(nuChange),
                    "| median: ", median(nuChange),
                    "| mean: ",   mean(nuChange) ),
              logLevel = 1)
    }

    return(objCOTAN)
  }
)


#'estimateDispersionNuBisection
#'
#' Estimates the 'dispersion' and 'nu' field of a `COTAN` object by running
#' sequentially a bisection for each parameter. This allows to enforce
#' `mean('nu') == 1` assumption
#' Assumes [estimateNuLinear()] being run already
#'
#' @param objCOTAN a `COTAN` object
#' @param step number of genes to solve in batch in a single core. Default is 256.
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param cores number of cores to use. Default is 1.
#'
#' @return the updated `COTAN` object
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
    logThis("Estimate 'dispersion'/'nu': START", logLevel = 2)

    if (is_empty(getNu(objCOTAN))) {
      stop("'nu' vector needs to be initialised as initial guess")
    }

    iter <- 1
    repeat {
      objCOTAN <- estimateDispersionBisection(objCOTAN, step = step,
                                              threshold = threshold,
                                              maxIterations = maxIterations,
                                              cores = cores)

      objCOTAN <- estimateNuBisection(objCOTAN, step = step,
                                      threshold = threshold,
                                      maxIterations = maxIterations,
                                      cores = cores)

      meanNu <- mean(getNu(objCOTAN))

      if (!is.finite(meanNu)) {
        stopifnot("Cannot have infinite mean 'nu' only after the first loop" <- (iter == 1))
        warning(paste0("Infinite 'nu' found: one of the cells expressed all genes\n",
                       " Returning values after a single loop"))
        return(objCOTAN)
      }

      logThis(paste("Nu mean:", meanNu), logLevel = 3)

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

    logThis("Estimate dispersion/nu: DONE", logLevel = 2)

    return(objCOTAN)
  }
)
