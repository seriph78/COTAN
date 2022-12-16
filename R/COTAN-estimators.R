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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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

    {
      oldLambda <- getMetadataGenes(objCOTAN)[["lambda"]]
      if (!identical(lambda, oldLambda)) {
        # flag the coex slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[5], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[6], FALSE)
      }
    }

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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
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

    {
      oldNu <-  getMetadataCells(objCOTAN)[["nu"]]
      if (!identical(nu, oldNu)) {
        # flag the coex slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[5], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[6], FALSE)
      }
    }

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
#'   [estimateNuLinear()] being already run.
#'
#'   It needs to be run after [clean()]
#'
#' @param objCOTAN A `COTAN` object
#' @param threshold minimal solution precision
#' @param cores number of cores to use. Default is 1.
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param chunkSize number of genes to solve in batch in a single core. Default
#'   is 1024.
#'
#' @returns The updated `COTAN` object
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @export
#'
#' @examples
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#' dispersion <- getDispersion(objCOTAN)
#'
#' @rdname estimateDispersionBisection
#'
setMethod(
  "estimateDispersionBisection",
  "COTAN",
  function(objCOTAN, threshold = 0.001, cores = 1,
           maxIterations = 100, chunkSize = 1024) {
    logThis("Estimate dispersion: START", logLevel = 2)

    if (Sys.info()['sysname'] == "Windows" && cores != 1) {
      warning(paste0("On windows the numebr of cores used will be 1!",
                     " Multicore is not supported."))
      cores <- 1
    }

    genes <- getGenes(objCOTAN)
    sumZeros <- set_names(getNumCells(objCOTAN) -
                            rowSums(getZeroOneProj(objCOTAN)[genes, , drop = FALSE]),
                          genes)
    lambda <- getLambda(objCOTAN)
    nu <- getNu(objCOTAN)

    dispList <- list()

    spIdx <- parallel::splitIndices(length(genes), ceiling(length(genes) / chunkSize))

    spGenes = lapply(spIdx, function(x) genes[x])

    numSplits <- length(spGenes)
    splitStep <- max(16, cores * 2)

    gc()

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
                  sumZeros = sumZeros,
                  lambda = lambda,
                  nu = nu,
                  threshold = threshold,
                  maxIterations = maxIterations,
                  mc.cores = cores)
      }
      else {
        res  <- lapply(
                  spGenes[pBegin:pEnd],
                  parallelDispersionBisection,
                  sumZeros = sumZeros,
                  lambda = lambda,
                  nu = nu,
                  threshold = threshold,
                  maxIterations = maxIterations)
      }

      dispList <- append(dispList, res)
      rm(res)

      pBegin <- pEnd + 1
    }
    logThis("Estimate dispersion: DONE", logLevel = 2)

    gc()

    dispersion <- unlist(dispList, recursive = TRUE, use.names = FALSE)
    {
      oldDispersion <-  getMetadataGenes(objCOTAN)[["dispersion"]]
      if (!identical(dispersion, oldDispersion)) {
        # flag the coex slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[5], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[6], FALSE)
      }
    }

    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes, dispersion,
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
#' @param threshold minimal solution precision
#' @param cores number of cores to use. Default is 1.
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param chunkSize number of genes to solve in batch in a single core. Default
#'   is 1024.
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
#' data("raw.dataset")
#' objCOTAN <- COTAN(raw = raw.dataset)
#' objCOTAN <- clean(objCOTAN)
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 12)
#' objCOTAN <- estimateNuBisection(objCOTAN, cores = 12)
#' nu <- getNu(objCOTAN)
#'
#' @rdname estimateNuBisection
#'
setMethod(
  "estimateNuBisection",
  "COTAN",
  function(objCOTAN, threshold = 0.001, cores = 1,
           maxIterations = 100, chunkSize = 1024) {
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
    sumZeros <- set_names(getNumGenes(objCOTAN) -
                            colSums(zeroOneMatrix[, cells, drop = FALSE]),
                          cells)
    lambda <- getLambda(objCOTAN)
    dispersion <- getDispersion(objCOTAN)
    initialGuess <- getNu(objCOTAN)

    nuList <- list()

    spIdx <- parallel::splitIndices(length(cells), ceiling(length(cells) / chunkSize))

    spCells = lapply(spIdx, function(x) cells[x])

    numSplits <- length(spCells)
    splitStep <- max(16, cores * 2)

    gc()

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
                  sumZeros = sumZeros,
                  lambda = lambda,
                  dispersion = dispersion,
                  initialGuess = initialGuess,
                  threshold = threshold,
                  maxIterations = maxIterations,
                  mc.cores = cores)
      }
      else {
        res  <- lapply(
                  spCells[pBegin:pEnd],
                  parallelNuBisection,
                  sumZeros = sumZeros,
                  lambda = lambda,
                  dispersion = dispersion,
                  initialGuess = initialGuess,
                  threshold = threshold,
                  maxIterations = maxIterations)
      }

      nuList <- append(nuList, res)
      rm(res)

      pBegin <- pEnd + 1
    }

    gc()
    logThis("Estimate nu: DONE", logLevel = 2)

    nu <- unlist(nuList, recursive = TRUE, use.names = FALSE)
    {
      if (!identical(nu, initialGuess)) {
        # flag the coex slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[5], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[6], FALSE)
      }
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, nu,
                                        "nu", cells)

    {
      nuChange <- abs(nu - initialGuess)
      logThis(paste("nu change (abs)",
                    "| max:",     max(nuChange),
                    "| median: ", median(nuChange),
                    "| mean: ",   mean(nuChange) ),
              logLevel = 1)
    }

    return(objCOTAN)
  }
)


#' stimateDispersionNuBisection
#'
#' Estimates the 'dispersion' and 'nu' field of a `COTAN` object by running
#' sequentially a bisection for each parameter. This allows to enforce
#' `mean('nu') == 1` assumption Assumes [estimateNuLinear()] being run already
#'
#' @param objCOTAN a `COTAN` object
#' @param threshold minimal solution precision
#' @param cores number of cores to use. Default is 1.
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param chunkSize number of genes to solve in batch in a single core. Default
#'   is 1024.
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
  function(objCOTAN, threshold = 0.001, cores = 1,
           maxIterations = 100, chunkSize = 1024, enforceNuAverageToOne = FALSE) {
    logThis("Estimate 'dispersion'/'nu': START", logLevel = 2)

    # getNu() would show a warning when no 'nu' present
    if (is_empty(getMetadataCells(objCOTAN)[["nu"]])) {
      objCOTAN <- estimateNuLinear(objCOTAN)
    }

    sumZeros <- getNumCells(objCOTAN) - rowSums(getZeroOneProj(objCOTAN))

    iter <- 1
    repeat {
      # a smalle threshold is used in order to ensure the global convergence
      objCOTAN <- estimateDispersionBisection(objCOTAN,
                                              threshold = threshold / 10,
                                              cores = cores,
                                              maxIterations = maxIterations,
                                              chunkSize = chunkSize)

      gc()

      objCOTAN <- estimateNuBisection(objCOTAN,
                                      threshold = threshold / 10,
                                      cores = cores,
                                      maxIterations = maxIterations,
                                      chunkSize = chunkSize)

      gc()

      meanNu <- mean(getNu(objCOTAN))
      logThis(paste("Nu mean:", meanNu),
              logLevel = (if(enforceNuAverageToOne) {1} else {3}))

      if (isTRUE(enforceNuAverageToOne)) {
        if (!is.finite(meanNu)) {
          stopifnot("Cannot have infinite mean 'nu' only after the first loop" <- (iter == 1))
          warning(paste0("Infinite 'nu' found: one of the cells expressed all genes\n",
                         " Setting 'enforceNuAverageToOne <- FALSE'"))
          enforceNuAverageToOne <- FALSE
        } else {
          objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells,
                                              getNu(objCOTAN) / meanNu,
                                              "nu")
        }
      }

      if (iter >= maxIterations) {
        warning(paste("Max number of outer iterations", maxIterations,
                      "reached while finding the solution"))
        break
      }

      genesMarginals <- rowSums(funProbZero(getDispersion(objCOTAN),
                                            calculateMu(objCOTAN)))

      marginalsErrors <- abs(genesMarginals - sumZeros)

      logThis(paste("Marginal errors | max:", max(marginalsErrors),
                    "| median", median(marginalsErrors),
                    "| mean:", mean(marginalsErrors)), logLevel = 2)

      gc()

      if (max(marginalsErrors) < threshold &&
          (isFALSE(enforceNuAverageToOne) || abs(meanNu - 1) < threshold)) {
        break
      }

      iter <- iter + 1
    }

    logThis("Estimate dispersion/nu: DONE", logLevel = 2)

    return(objCOTAN)
  }
)
