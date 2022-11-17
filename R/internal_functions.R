
funProbZero <- function(disp, mu) {
  ad <- abs(disp)
  (disp <= 0) * (exp(-(1 + ad) * mu)) +
  (disp >  0) * (1 + ad * mu)^(-1 / ad)
}

#' dispersionBisection
#'
#' private function invoked by 'estimateDispersion' for the estimation
#' of 'dispersion' field of a COTAN object via a bisection solver
#'
#' the goal is to find dispersion value that produces a difference between
#' the number of estimated and counted zeros close to 0
#'
#' @param gene name of the relevant gene
#' @param zeroOneMatrix raw data matrix changed to 0-1 matrix
#' @param muEstimator a matrix estimator of vector mu
#' @param threshold minimal solution precision
#'
#' @return r, data.frame(a, u)
dispersionBisection <- function(gene,
                                zeroOneMatrix,
                                muEstimator,
                                housekeepingGenes,
                                threshold = 0.001) {
  if(gene %in% housekeepingGenes) {
    return(NA)
  }

  sumZeros <- sum(zeroOneMatrix[gene, ] == 0)
  muEstimator <- muEstimator[gene, ]

  # we look for two dispersion values where the first leads to a
  # diffZeros negative and the second positive
  disp1 <- 0
  diff1 <- sum(funProbZero(disp1, muEstimator)) - sumZeros
  if (abs(diff1) <= threshold) {
    return(disp1)
  }

  disp2 <- -1 * sign(diff1) # we assume error is an increasing function of disp
  repeat {
    diff2 <- sum(funProbZero(disp2, muEstimator)) - sumZeros

    if (diff2 * diff1 < 0) {
      break
    }

    disp1 <- disp2 # disp is closer to producing 0
    diff1 <- diff2

    disp2 <- 2 * disp2 # we double at each step
  }

  # once we have found the two bounds to the dispersion value, we use bisection
  repeat {
    disp <- (disp1 + disp2) / 2
    diff <- sum(funProbZero(disp, muEstimator)) - sumZeros

    if (abs(diff) <= threshold) {
      return(disp)
    }

    # drop same sign diff point
    if (diff * diff2 > 0) {
      disp2 <- disp
      diff2 <- diff
    } else {
      disp1 <- disp
      diff1 <- diff
    }
  }
}


setGeneric("get.S", function(object) standardGeneric("get.S"))
setMethod("get.S","scCOTAN",
  function(object) {
      print("function to generate S ")
      if (is(class(object@coex)[1], "dtCMatrix") ||
          as.vector(class(object@coex)) %in% "dtCMatrix") {
          print("COTAN object in the old format! Converting...")
          object <- calculateCoex(object)
      }
      S <- (object@coex$values)^2 * getNumCells(object)
      return( list("genes" = object@coex$genes, "values" = S) )
  }
)


setGeneric("obs_ct", function(object) standardGeneric("obs_ct"))
setMethod("obs_ct","scCOTAN",
  function(object) {
    #---------------------------------------------------
    # Cells matrix : formed by row data matrix changed to 0-1 matrix
    cells <- getZeroOneProj(object)
    print("Generating contingency tables for observed data")
    cellsRowSums <- as.matrix(rowSums(cells))
    si_any <- do.call("cbind", replicate(length(rownames(cellsRowSums)),
                                         cellsRowSums, simplify = FALSE))

    colnames(si_any) = rownames(si_any)

    si_si <- observedContingencyYY(object)
    si_no <- si_any - si_si

    si_any <- t(si_any)
    no_si <- si_any - si_si

    no_no <- length(colnames(cells)) - (si_si + no_si + si_no)
    out <- list("yes_yes"=si_si,"no_yes"=no_si,"yes_no"=si_no,"no_no"=no_no)
    return(out)
  }
)


setGeneric("get.G", function(object) standardGeneric("get.G"))
setMethod(
  "get.G", "scCOTAN",
  function(object) {
    print("function to generate G ")
    noHKFlags <- flagNotHousekeepingGenes(object)
    ll <- obs_ct(object)

    ll$no_yes  <- ll$no_yes [noHKFlags, noHKFlags]
    ll$no_no   <- ll$no_no  [noHKFlags, noHKFlags]
    ll$yes_yes <- ll@yes_yes[noHKFlags, noHKFlags]
    ll$yes_no  <- ll$yes_no [noHKFlags, noHKFlags]

    est <- expectedContingencyTables(object, FALSE)
    for (i in est) {
      stopifnot("Some expected values are 0!" = !any(i == 0))
    }

    print("G estimation")

    t1 <- as.matrix(ll$yes_yes) * log(as.matrix(ll$yes_yes) /
      as.matrix(est$expectedYY))
    t1[which(as.matrix(ll$yes_yes) == 0)] <- 0

    t2 <- as.matrix(ll$no_no) * log(as.matrix(ll$no_no) /
      as.matrix(est$expectedNN))
    t2[which(as.matrix(ll$no_no) == 0)] <- 0

    t3 <- as.matrix(ll$yes_no) * log(as.matrix(ll$yes_no) /
      as.matrix(est$expectedYN))
    t3[which(as.matrix(ll$yes_no) == 0)] <- 0

    t4 <- as.matrix(ll$no_yes) * log(as.matrix(ll$no_yes) /
      as.matrix(est$expectedNY))
    t4[which(as.matrix(ll$no_yes) == 0)] <- 0

    G <- 2 * (t1 + t2 + t3 + t4)
    return(G)
  }
)
