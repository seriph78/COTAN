# ----- `COTAN` parameters' estimates methods -----

#'
#' @title Estimation of the `COTAN` model's parameters
#'
#' @description These functions are used to estimate the `COTAN` model's
#'   parameters. That is the average count for each gene (`lambda`) the average
#'   count for each cell (`nu`) and the `dispersion` parameter for each gene to
#'   match the probability of zero.
#'
#'   The estimator methods are named `Linear` if they can be calculated as a
#'   linear statistic of the raw data or `Bisection` if they are found via a
#'   parallel bisection solver.
#'
#' @name ParametersEstimations
NULL


### ------ estimateLambdaLinear -----

#'
#' @aliases estimateLambdaLinear
#'
#' @details `estimateLambdaLinear()` does a linear estimation of `lambda`
#'   (genes' counts averages)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `estimateLambdaLinear()` returns the updated `COTAN` object
#'
#' @importFrom Matrix mean
#' @importFrom Matrix rowMeans
#'
#' @export
#'
#' @examples
#' data("test.dataset")
#' objCOTAN <- COTAN(raw = test.dataset)
#'
#' objCOTAN <- estimateLambdaLinear(objCOTAN)
#' lambda <- getLambda(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "estimateLambdaLinear",
  "COTAN",
  function(objCOTAN) {
    lambda <- rowMeans(getRawData(objCOTAN), dims = 1L)

    if (TRUE) {
      oldLambda <- getMetadataGenes(objCOTAN)[["lambda"]]
      if (!identical(lambda, oldLambda)) {
        # flag the COEX slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["gsync"]], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["csync"]], FALSE)
      }
    }

    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes, lambda,
                                        "lambda", getGenes(objCOTAN))

    return(objCOTAN)
  }
)


### ------ estimateNuLinear -----

#' @aliases estimateNuLinear
#'
#' @details `estimateNuLinear()` does a linear estimation of `nu` (normalized
#'   cells' counts averages)
#'
#' @param objCOTAN a `COTAN` object
#'
#' @returns `estimateNuLinear()` returns the updated `COTAN` object
#'
#' @importFrom Matrix colSums
#'
#' @export
#'
#' @examples
#' objCOTAN <- estimateNuLinear(objCOTAN)
#' nu <- getNu(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "estimateNuLinear",
  "COTAN",
  function(objCOTAN) {
    # raw column sums divided by global average
    nu <- colSums(getRawData(objCOTAN), dims = 1L)
    nu <- nu / mean(nu)

    if (TRUE) {
      oldNu <- getMetadataCells(objCOTAN)[["nu"]]
      if (!identical(nu, oldNu)) {
        # flag the COEX slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["gsync"]], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["csync"]], FALSE)
      }
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, nu,
                                        "nu", getCells(objCOTAN))

    return(objCOTAN)
  }
)


# ------- `COTAN` clusterization data accessors ------

#' @title Handling cells' *clusterization* and related functions
#'
#' @description These functions manage the *clusterizations* and their
#'   associated *cluster* `COEX` `data.frame`s.
#'
#'   A *clusterization* is any partition of the cells where to each cell it is
#'   assigned a **label**; a group of cells with the same label is called
#'   *cluster*.
#'
#'   For each *cluster* is also possible to define a `COEX` value for each gene,
#'   indicating its increased or decreased expression in the *cluster* compared
#'   to the whole background. A `data.frame` with these values listed in a
#'   column for each *cluster* is stored separately for each *clusterization* in
#'   the `clustersCoex` member.
#'
#'   The formulae for this *In/Out* `COEX` are similar to those used in the
#'   [calculateCoex()] method, with the **role** of the second gene taken by the
#'   *In/Out* status of the cells with respect to each *cluster*.
#'
#' @name HandlingClusterizations
NULL


### ------ estimateNuLinearByCluster -----

#' @aliases estimateNuLinearByCluster
#'
#' @details `estimateNuLinearByCluster()` does a linear estimation of `nu`:
#'   cells' counts averages normalized *cluster* by *cluster*
#'
#' @param objCOTAN a `COTAN` object
#' @param clName The name of the *clusterization*. If not given the last
#'   available *clusterization* will be used, as it is probably the most
#'   significant!
#' @param clusters A *clusterization* to use. If given it will take precedence
#'   on the one indicated by `clName`
#'
#' @returns `estimateNuLinearByCluster()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#' @importFrom rlang set_names
#'
#' @importFrom Matrix colSums
#'
#' @export
#'
#' @rdname HandlingClusterizations
#'
setMethod(
  "estimateNuLinearByCluster",
  "COTAN",
  function(objCOTAN, clName = "", clusters = NULL) {
    c(clName, clusters) %<-%
      normalizeNameAndLabels(objCOTAN, name = clName,
                             labels = clusters, isCond = FALSE)

    # raw column sums divided by cluster average
    nu <- colSums(getRawData(objCOTAN), dims = 1L)

    for (cl in levels(clusters)) {
      c <- clusters == cl
      nu[c] <- nu[c] / mean(nu[c])
    }

    if (TRUE) {
      oldNu <- getMetadataCells(objCOTAN)[["nu"]]
      if (!identical(nu, oldNu)) {
        # flag the COEX slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["gsync"]], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["csync"]], FALSE)
      }
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, nu,
                                        "nu", getCells(objCOTAN))

    return(objCOTAN)
  }
)


# local utility wrapper for parallel estimation of dispersion
runDispSolver <- function(genesBatches, sumZeros, lambda, nu,
                          threshold, maxIterations, cores) {

  ## forward caller that allows for logging progress
  worker <- function(genesBatch, sumZeros, lambda, nu,
                     threshold, maxIterations) {
    tryCatch({
      return(parallelDispersionBisection(genes         = genesBatch,
                                         sumZeros      = sumZeros,
                                         lambda        = lambda,
                                         nu            = nu,
                                         threshold     = threshold,
                                         maxIterations = maxIterations))
    }, error = function(e) {
      structure(list(e), class = "try-error")
    })
  }

  if (cores != 1L) {
    ## create a tiny private env `mini` and a worker within it:
    ## this reduces spawning memory foot-print
    mini <- new.env(parent = baseenv())
    mini$parallelDispersionBisection <- parallelDispersionBisection
    mini$rowsums <- Rfast::rowsums
    mini$funProbZero <- funProbZero
    mini$logThis <- logThis
    mini$worker <- worker

    ## set on the copy in `mini` to only use `mini`
    environment(mini$parallelDispersionBisection) <- mini
    environment(mini$worker) <- mini

    res <- parallel::mclapply(
      genesBatches,
      worker,
      sumZeros = sumZeros,
      lambda = lambda,
      nu = nu,
      threshold = threshold,
      maxIterations = maxIterations,
      mc.cores = cores,
      mc.preschedule = FALSE)

    # spawned errors are stored as try-error classes
    resError <- unlist(lapply(res, inherits, "try-error"))
    if (any(resError)) {
      stop(res[[which(resError)[[1L]]]], call. = FALSE)
    }

    return(res)
  } else {
    res <- lapply(genesBatches,
                  parallelDispersionBisection,
                  sumZeros = sumZeros,
                  lambda = lambda,
                  nu = nu,
                  threshold = threshold,
                  maxIterations = maxIterations)

    return(res)
  }
}


### ------ estimateDispersionBisection -----

#' @aliases estimateDispersionBisection
#'
#' @details `estimateDispersionBisection()` estimates the negative binomial
#'   dispersion factor for each gene (`dispersion`). Determines the value such
#'   that, for each gene, the probability of zero count matches the number of
#'   observed zeros. It assumes [estimateNuLinear()] being already run.
#'
#' @param objCOTAN a `COTAN` object
#' @param threshold minimal solution precision
#' @param cores number of cores to use. Default is 1.
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param chunkSize number of elements to solve in batch in a single core.
#'   Default is 1024.
#'
#' @returns `estimateDispersionBisection()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_null
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @importFrom parallelly supportsMulticore
#' @importFrom parallelly availableCores
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @export
#'
#' @examples
#' objCOTAN <- estimateDispersionBisection(objCOTAN, cores = 6L)
#' dispersion <- getDispersion(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "estimateDispersionBisection",
  "COTAN",
  function(objCOTAN, threshold = 0.001, cores = 1L,
           maxIterations = 100L, chunkSize = 1024L) {
    logThis("Estimate `dispersion`: START", logLevel = 2L)

    cores <- handleMultiCore(cores)

    genes <- getGenes(objCOTAN)
    sumZeros <- getNumCells(objCOTAN) - getNumOfExpressingCells(objCOTAN)[genes]

    lambda <- suppressWarnings(getLambda(objCOTAN))
    assert_that(!is_empty(lambda),
                msg = "`lambda` must not be empty, estimate it")

    nu <- suppressWarnings(getNu(objCOTAN))
    assert_that(!is_empty(nu),
                msg = "`nu` must not be empty, estimate it")

    dispList <- list()

    spIdx <- parallel::splitIndices(length(genes),
                                    ceiling(length(genes) / chunkSize))

    spGenes <- lapply(spIdx, function(x) genes[x])

    cores <- min(cores, length(spGenes))

    logThis(paste0("Executing ", length(spGenes), " genes batches"),
            logLevel = 3L)

    dispList <- runDispSolver(spGenes,
                              sumZeros      = sumZeros,
                              lambda        = lambda,
                              nu            = nu,
                              threshold     = threshold,
                              maxIterations = maxIterations,
                              cores         = cores)

    gc()

    dispersion <- unlist(dispList, recursive = TRUE, use.names = FALSE)
    if (TRUE) {
      oldDispersion <-  getMetadataGenes(objCOTAN)[["dispersion"]]
      if (!identical(dispersion, oldDispersion)) {
        # flag the COEX slots are out of sync (if any)!
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["gsync"]], FALSE)
        objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                               datasetTags()[["csync"]], FALSE)
      }
    }

    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes, dispersion,
                                        "dispersion", genes)

    goodPos <- is.finite(getDispersion(objCOTAN))
    logThis(paste("`dispersion`",
                  "| min:", min(getDispersion(objCOTAN)[goodPos]),
                  "| max:", max(getDispersion(objCOTAN)[goodPos]),
                  "| % negative:", (sum(getDispersion(objCOTAN) < 0.0) * 100.0 /
                                    getNumGenes(objCOTAN))),
            logLevel = 1L)

    return(objCOTAN)
  }
)



# local utility wrapper for parallel estimation of nu
runNuSolver <- function(cellsBatches, sumZeros, lambda, dispersion,
                        initialGuess, threshold, maxIterations, cores) {

  ## forward caller that allows for logging progress
  worker <- function(batch, sumZeros, lambda, dispersion,
                     initialGuess, threshold, maxIterations) {
    tryCatch({
      return(parallelNuBisection(batch,
                                 sumZeros      = sumZeros,
                                 lambda        = lambda,
                                 dispersion    = dispersion,
                                 initialGuess  = initialGuess,
                                 threshold     = threshold,
                                 maxIterations = maxIterations))
    }, error = function(e) {
      structure(list(e), class = "try-error")
    })
  }

  if (cores != 1L) {
    # create a tiny private env and a worker within:
    # it reduces spawning memory foot-print
    mini <- new.env(parent = baseenv())
    mini$parallelNuBisection <- parallelNuBisection
    mini$assert_that <- assertthat::assert_that
    mini$colsums <- Rfast::colsums
    mini$funProbZero <- funProbZero
    mini$logThis <- logThis
    mini$worker <- worker

    ## set on the copy in `mini` to only use `mini`
    environment(mini$parallelNuBisection) <- mini
    environment(mini$worker) <- mini

    res <- parallel::mclapply(
      cellsBatches,
      parallelNuBisection,
      sumZeros = sumZeros,
      lambda = lambda,
      dispersion = dispersion,
      initialGuess = initialGuess,
      threshold = threshold,
      maxIterations = maxIterations,
      mc.cores = cores,
      mc.preschedule = FALSE)

    # spawned errors are stored as try-error classes
    resError <- unlist(lapply(res, inherits, "try-error"))
    if (any(resError)) {
      stop(res[[which(resError)[[1L]]]], call. = FALSE)
    }

    return(res)
  } else {
    res <- lapply(cellsBatches,
                  parallelNuBisection,
                  sumZeros = sumZeros,
                  lambda = lambda,
                  dispersion = dispersion,
                  initialGuess = initialGuess,
                  threshold = threshold,
                  maxIterations = maxIterations)

    return(res)
  }
}


### ------ estimateNuBisection -----

#' @aliases estimateNuBisection
#'
#' @details `estimateNuBisection()` estimates the `nu` vector of a `COTAN`
#'   object by bisection. It determines the `nu` parameters such that, for each
#'   cell, the probability of zero counts matches the number of observed zeros.
#'   It assumes [estimateDispersionBisection()] being already run. Since this
#'   breaks the assumption that the average `nu` is one, it is recommended not
#'   to run this in isolation but use [estimateDispersionNuBisection()] instead.
#'
#' @param objCOTAN a `COTAN` object
#' @param threshold minimal solution precision
#' @param cores number of cores to use. Default is 1.
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param chunkSize number of elements to solve in batch in a single core.
#'   Default is 1024.
#'
#' @returns `estimateNuBisection()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_null
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom parallel mclapply
#' @importFrom parallel splitIndices
#'
#' @importFrom parallelly supportsMulticore
#' @importFrom parallelly availableCores
#'
#' @importFrom stats median
#'
#' @importFrom zeallot %<-%
#' @importFrom zeallot %->%
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "estimateNuBisection",
  "COTAN",
  function(objCOTAN, threshold = 0.001, cores = 1L,
           maxIterations = 100L, chunkSize = 1024L) {
    logThis("Estimate `nu`: START", logLevel = 2L)

    cores <- handleMultiCore(cores)

    # parameters estimation


    cells <- getCells(objCOTAN)
    sumZeros <- getNumGenes(objCOTAN) - getNumExpressedGenes(objCOTAN)

    lambda <- suppressWarnings(getLambda(objCOTAN))
    assert_that(!is_empty(lambda),
                msg = "`lambda` must not be empty, estimate it")

    nu <- suppressWarnings(getNu(objCOTAN))
    assert_that(!is_empty(nu),
                msg = "`nu` must not be empty, estimate it")

    dispersion <- suppressWarnings(getDispersion(objCOTAN))
    assert_that(!is_empty(dispersion),
                msg = "`dispersion` must not be empty, estimate it")

    initialGuess <- nu

    nuList <- list()

    spIdx <- parallel::splitIndices(length(cells),
                                    ceiling(length(cells) / chunkSize))

    spCells <- lapply(spIdx, function(x) cells[x])

    cores <- min(cores, length(spCells))

    logThis(paste0("Executing ", length(spCells), " cells batches"),
            logLevel = 3L)

    nuList <- runNuSolver(spCells,
                          sumZeros = sumZeros,
                          lambda = lambda,
                          dispersion = dispersion,
                          initialGuess = initialGuess,
                          threshold = threshold,
                          maxIterations = maxIterations,
                          cores = cores)

    gc()
    logThis("Estimate `nu`: DONE", logLevel = 2L)

    nu <- unlist(nuList, recursive = TRUE, use.names = FALSE)
    if (!identical(nu, initialGuess)) {
      # flag the COEX slots are out of sync (if any)!
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["gsync"]], FALSE)
      objCOTAN@metaDataset <- updateMetaInfo(objCOTAN@metaDataset,
                                             datasetTags()[["csync"]], FALSE)
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, nu,
                                        "nu", cells)

    if (TRUE) {
      nuChange <- abs(nu - initialGuess)
      logThis(paste("`nu` change (abs)",
                    "| max:",     max(nuChange),
                    "| median: ", median(nuChange),
                    "| mean: ",   mean(nuChange)),
              logLevel = 1L)
    }

    return(objCOTAN)
  }
)


### ------ estimateDispersionNuBisection -----

#' @aliases estimateDispersionNuBisection
#'
#' @details `estimateDispersionNuBisection()` estimates the `dispersion` and
#'   `nu` field of a `COTAN` object by running sequentially a bisection for each
#'   parameter.
#'
#' @param objCOTAN a `COTAN` object
#' @param threshold minimal solution precision
#' @param cores number of cores to use. Default is 1.
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param chunkSize number of elements to solve in batch in a single core.
#'   Default is 1024.
#' @param enforceNuAverageToOne a Boolean on whether to keep the average `nu`
#'   equal to 1
#'
#' @returns `estimateDispersionNuBisection()` returns the updated `COTAN` object
#'
#' @export
#'
#' @importFrom rlang is_empty
#'
#' @importFrom Rfast rowsums
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#' objCOTAN <- estimateDispersionNuBisection(objCOTAN, cores = 6L,
#'                                           enforceNuAverageToOne = TRUE)
#' nu <- getNu(objCOTAN)
#' dispersion <- getDispersion(objCOTAN)
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "estimateDispersionNuBisection",
  "COTAN",
  function(objCOTAN, threshold = 0.001, cores = 1L,
           maxIterations = 100L, chunkSize = 1024L,
           enforceNuAverageToOne = TRUE) {
    logThis("Estimate `dispersion`/`nu`: START", logLevel = 2L)

    # getNu() would show a warning when no `nu` present
    if (is_empty(suppressWarnings(getNu(objCOTAN)))) {
      objCOTAN <- estimateNuLinear(objCOTAN)
    }

    sumZeros <- getNumCells(objCOTAN) - getNumOfExpressingCells(objCOTAN)

    iter <- 1L
    repeat {
      # a smaller threshold is used in order to ensure the global convergence
      objCOTAN <- estimateDispersionBisection(objCOTAN,
                                              threshold = threshold / 10.0,
                                              cores = cores,
                                              maxIterations = maxIterations,
                                              chunkSize = chunkSize)

      gc()

      objCOTAN <- estimateNuBisection(objCOTAN,
                                      threshold = threshold / 10.0,
                                      cores = cores,
                                      maxIterations = maxIterations,
                                      chunkSize = chunkSize)

      gc()

      meanNu <- mean(getNu(objCOTAN))
      logThis(paste("`nu` mean:", meanNu), logLevel = 2L)

      if (isTRUE(enforceNuAverageToOne)) {
        if (is.finite(meanNu)) {
          objCOTAN@metaCells <-
            setColumnInDF(objCOTAN@metaCells, getNu(objCOTAN) / meanNu, "nu")
        } else {
          assert_that(iter == 1L,
                      msg = paste("It can happen to have infinite mean `nu`",
                                  "only after the first loop"))
          warning("Infinite `nu` found: one of the cells expressed all genes\n",
                  " Setting 'enforceNuAverageToOne <- FALSE'")
          enforceNuAverageToOne <- FALSE
        }
      }

      if (iter >= maxIterations) {
        warning("Max number of outer iterations", maxIterations,
                "reached while finding the solution")
        break
      }

      genesMarginals <- rowsums(getProbabilityOfZero(objCOTAN),
                                parallel = TRUE)

      marginalsErrors <- abs(genesMarginals - sumZeros)

      logThis(paste("Marginal errors | max:", max(marginalsErrors),
                    "| median", median(marginalsErrors),
                    "| mean:", mean(marginalsErrors)), logLevel = 2L)

      gc()

      if (max(marginalsErrors) < threshold &&
          (isFALSE(enforceNuAverageToOne) || abs(meanNu - 1.0) < threshold)) {
        break
      }

      iter <- iter + 1L
    }

    logThis("Estimate `dispersion`/`nu`: DONE", logLevel = 2L)

    return(objCOTAN)
  }
)


### ------ estimateDispersionNuNlminb -----

#' @aliases estimateDispersionNuNlminb
#'
#' @details `estimateDispersionNuNlminb()` estimates the `nu` and
#'   `dispersion` parameters to minimize the discrepancy between the observed
#'   and expected probability of zero. It uses the [stats::nlminb()] solver,
#'   but since the joint parameters have too high dimensionality, it converges
#'   too slowly to be actually useful in real cases.
#'
#' @param objCOTAN a `COTAN` object
#' @param threshold minimal solution precision
#' @param maxIterations max number of iterations (avoids infinite loops)
#' @param chunkSize number of elements to solve in batch in a single core.
#'   Default is 1024.
#' @param enforceNuAverageToOne a Boolean on whether to keep the average `nu`
#'   equal to 1
#'
#' @returns `estimateDispersionNuNlminb()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_null
#' @importFrom rlang is_empty
#'
#' @importFrom assertthat assert_that
#'
#' @importFrom Rfast rowsums
#' @importFrom Rfast colsums
#'
#' @importFrom stats nlminb
#'
#' @rdname ParametersEstimations
#'
setMethod(
  "estimateDispersionNuNlminb",
  "COTAN",
  function(objCOTAN, threshold = 0.001,
           maxIterations = 50L, chunkSize = 1024L,
           enforceNuAverageToOne = TRUE) {
    logThis("Estimate `dispersion`/`nu`: START", logLevel = 2L)

    lambda <- suppressWarnings(getLambda(objCOTAN))
    assert_that(!is_empty(lambda),
                msg = "`lambda` must not be empty, estimate it")

    nu <- suppressWarnings(getNu(objCOTAN))
    if (is_empty(nu)) {
      objCOTAN <- estimateNuLinear(objCOTAN)
      nu <- getNu(objCOTAN)
    }

    zeroGenes <- getNumCells(objCOTAN) - getNumOfExpressingCells(objCOTAN)
    zeroCells <- getNumGenes(objCOTAN) - getNumExpressedGenes(objCOTAN)

    assert_that(all(zeroGenes != 0L), all(zeroCells != 0L),
                msg = paste("Method cannot handle fully-expressed genes",
                            "or fully-expressing cells yet"))


    errorFunction <- function(params, lambda, zeroGenes,
                              zeroCells, enforceNuAverageToOne) {
      numGenes <- length(lambda)
      dispersion <- params[1L:numGenes]
      nu <- exp(params[(numGenes + 1L):length(params)])

      probZero <- funProbZero(dispersion, lambda %o% nu)

      diffGenes <- rowsums(probZero, parallel = TRUE) - zeroGenes
      diffCells <- colsums(probZero, parallel = TRUE) - zeroCells

      diffNu <- 0.0
      if (enforceNuAverageToOne) {
        diffNu <- mean(nu) - 1.0
      }

      diff <- c(diffGenes, diffCells, diffNu)
      error <- sqrt(sum(diff^2L) / length(diff))
      return(error)
    }

    solution <- nlminb(
      start = c(rep(0.0, getNumGenes(objCOTAN)), log(getNu(objCOTAN))),
      lambda = lambda, zeroGenes = zeroGenes, zeroCells = zeroCells,
      enforceNuAverageToOne = enforceNuAverageToOne,
      objective = errorFunction,
      control = list(iter.max = maxIterations, abs.tol = threshold,
                     trace = 2L, step.min = 0.001, step.max = 0.1))

    numGenes <- length(lambda)
    dispersion <- solution[["par"]][1L:numGenes]
    nu <- exp(solution[["par"]][(numGenes + 1L):length(solution[["par"]])])

    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes, dispersion,
                                        "dispersion", getGenes(objCOTAN))

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, nu,
                                        "nu", getCells(objCOTAN))

    logThis("Estimate `dispersion`/`nu`: DONE", logLevel = 2L)

    return(objCOTAN)
  }
)
