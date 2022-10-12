#' estimates by bisection the 'dispersion' field of a COTAN object
#' @param objCOTAN a COTAN object
#' @param cores number of cores to use. Default is 1
#' @param step number of genes taken into account at each application. 
#' Default is 200
#' of 'dispersionBisection' for parallel computing
#' @return A COTAN object
#' @export
setMethod(
  "estimateDispersionBisection",
  "COTAN",
  function(objCOTAN, cores, step) {
    if (Sys.info()["sysname"] == "Windows") {
      message("On windows the numebr of cores used will be 1!
              Multicore is not supported.")
      cores <- 1
    }

    # house keeping genes
    # this genes are removed from zeroOneMatrix and muEstimator
    if (is_empty(objCOTAN@hkGenes)) {
      objCOTAN <- housekeepingGenes(objCOTAN)
    }

    # taken the information if the values are greater than 0 or not
    zeroOneMatrix <- objCOTAN@raw
    zeroOneMatrix[zeroOneMatrix > 0] <- 1
    zeroOneMatrix[zeroOneMatrix <= 0] <- 0
    zeroOneMatrix <-
      zeroOneMatrix[!rownames(zeroOneMatrix) %in% objCOTAN@hkGenes, ]

    # estimator of mu
    muEstimator <- estimateMu(objCOTAN)
    muEstimator <- muEstimator[!rownames(muEstimator) %in% objCOTAN@hkGenes, ]
    muEstimator <- as.matrix(muEstimator)

    # tot accumulates the results calculated in parallel
    tot <- list()
    
    # idex
    i <- 1
    limit <- length(rownames(muEstimator))

    while (i <= limit) {
      if ((i + step) <= limit) {
        # parallel computing
        toAppend <- parallel::mclapply(
          rownames(muEstimator)[i:(i + step)], # vector
          dispersionBisection, # function
          zeroOneMatrix = zeroOneMatrix, # parameter
          muEstimator = muEstimator, # parameter
          mc.cores = cores # cores
        )
      } else {
        # last block
        toAppend <- parallel::mclapply(
          rownames(muEstimator)[i:limit], # vector
          dispersionBisection, # function
          zeroOneMatrix = zeroOneMatrix, # parameter
          muEstimator = muEstimator, # parameter
          mc.cores = cores # cores
        )
      }
      # raccogliere i risultati
      tot <- append(tot, toAppend)
      i <- i + step + 1
    }
    
    tot2 <- tot[[1]]
    for (j in 2:length(tot)) {
      tot2 <- rbind(tot2, tot[[j]])
    }
    
    objCOTAN@dispersion <- tot2$dispersion
    names(objCOTAN@dispersion) <- rownames(tot2)
    
    return(objCOTAN)
  }
)
