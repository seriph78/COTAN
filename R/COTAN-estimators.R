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
#' @param cores number of cores to use. Default is 1.
#' @param step number of genes to solve for at once. Default is 200.
#'
#' @return the updated COTAN object
#'
#' @importFrom parallel mclapply
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
  function(objCOTAN, cores = 1, step = 200) {
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

    numGenes <- length(genes)
    pBegin <- 1
    while (pBegin <= numGenes) {
      pEnd <- pBegin + step - 1 # pEnd is always a multiple of step
      if (pEnd >= numGenes) {
        print("Final genes' trance!")
        pEnd = numGenes
      }

      res <- parallel::mclapply(
        genes[pBegin:pEnd],
        dispersionBisection,
        zeroOneMatrix = zeroOneMatrix,
        muEstimator = muEstimator,
        housekeepingGenes = housekeepingGenes,
        mc.cores = cores)

      dispList <- append(dispList, res)
      rm(res)

      pBegin <- pEnd + 1
      if (((pBegin - 1) %% (10 * step)) == 0) { # print a message each 10*step genes
        print(paste("Next gene:", genes[pBegin], "number", pBegin))
      }
    }

    gc()
    print("Estimate dispersion: DONE")

    objCOTAN@dispersion <- unlist(dispList, recursive = FALSE, use.names = FALSE)
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

