

setGeneric("get.S", function(object) standardGeneric("get.S"))
setMethod("get.S","scCOTAN",
  function(object) {
      print("function to generate S ")
      if(is(class(object@coex)[1], "dtCMatrix") | (as.vector(class(object@coex)) %in% "dtCMatrix")){
          print("COTAN object in the old format! Converting...")
          object <- get.coex(object)
      }
      S <- (object@coex$values)^2 * getNumCells(object)
      S <- list("genes"=object@coex$genes,"values"=S)
      return(S)
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
    hk <- object@hk
    ll <- obs_ct(object)

    ll$no_yes  <- ll$no_yes [!rownames(ll$no_yes)  %in% hk, !colnames(ll$no_yes)  %in% hk]
    ll$no_no   <- ll$no_no  [!rownames(ll$no_no)   %in% hk, !colnames(ll$no_no)   %in% hk]
    ll$yes_yes <- ll@yes_yes[!rownames(ll@yes_yes) %in% hk, !colnames(ll@yes_yes) %in% hk]
    ll$yes_no  <- ll$yes_no [!rownames(ll$yes_no)  %in% hk, !colnames(ll$yes_no)  %in% hk]

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
