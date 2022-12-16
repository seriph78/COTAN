
#' DEAOnClusters
#'
#' @description This function is used to run the differential expression
#'   analysis using the `COTAN` contingency tables on each cluster in the given
#'   clusterization
#'
#' @param objCOTAN A COTAN object
#' @param clusterization a `vector` or `factor` with the cell clusterization to
#'   be analysed. If empty, the function will use the last clusterization stored
#'   in the `objCOTAN`
#'
#' @return a list with two objects:
#'   * the correlation matrix for the genes in each cluster,
#'   * the corresponding p-values matrix.
#'
#' @export
#'
#' @examples
#'
#' @rdname DEAOnClusters
#'
DEAOnClusters <- function(objCOTAN, clusterization = NULL) {
  if (is_empty(clusterization)) {
    clusterization <- getClusterizationData(objCOTAN)[["clusters"]]
  }

  stopifnot("Passed/retrieved clusterization has the wrong size" <-
              (length(clusterization) == getNumCells(objCOTAN)))

  stopifnot("Some cells in the clusterization are not part of the 'COTAN' object" <-
              all(names(clusterization) %in% getCells(objCOTAN)))

  logThis("Differential Expression Analysis - START", logLevel = 2)

  clustersList <- toClustersList(clusterization)

  noHKFlags <- flagNotHousekeepingGenes(objCOTAN)

  zeroOne <- getZeroOneProj(objCOTAN)[noHKFlags, ]

  muEstimator <- calculateMu(objCOTAN)[noHKFlags, ]

  M <- funProbZero(getDispersion(objCOTAN)[noHKFlags], muEstimator) # matrix of 0 probabilities

  clusters_coex <- data.frame()
  clusters_pval <- data.frame()

  for (cl in names(clustersList)) {
    gc()
    logThis(paste0("Analysis of cluster: '", cl, "' - START"), logLevel = 3)

    cellsIn <- getCells(objCOTAN) %in%  clustersList[[cl]]

                  yes_in  <- rowSums(zeroOne[, cellsIn])
                  yes_out <- rowSums(zeroOne[,!cellsIn])
                  no_in   <- sum( cellsIn) - yes_in
                  no_out  <- sum(!cellsIn) - yes_out

                  estimator_no_in   <- rowSums(M[, cellsIn])
                  estimator_no_out  <- rowSums(M[,!cellsIn])

                  estimator_yes_in  <- sum( cellsIn) - estimator_no_in
                  estimator_yes_out <- sum(!cellsIn) - estimator_no_out

                  exp_yes <- rowSums(zeroOne)
                  exp_no <- ncol(zeroOne) - exp_yes

                  new_estimator_no_in <- estimator_no_in
                  new_estimator_no_in[new_estimator_no_in < 1] <- 1
                  new_estimator_no_out <- estimator_no_out
                  new_estimator_no_out[new_estimator_no_out < 1] <- 1

                  new_estimator_yes_in <- estimator_yes_in
                  new_estimator_yes_in[new_estimator_yes_in < 1] <- 1
                  new_estimator_yes_out <- estimator_yes_out
                  new_estimator_yes_out[new_estimator_yes_out < 1] <- 1

                  dif_no_in <- (as.matrix(no_in) - estimator_no_in)**2/new_estimator_no_in
                  dif_no_out <- (as.matrix(no_out) - estimator_no_out)**2/new_estimator_no_out
                  dif_yes_in <- (as.matrix(yes_in) - estimator_yes_in)**2/new_estimator_yes_in
                  dif_yes_out <- (as.matrix(yes_out) - estimator_yes_out)**2/new_estimator_yes_out

                  S <- dif_no_in + dif_no_out + dif_yes_in + dif_yes_out

                  rm(dif_yes_out,dif_no_in,dif_no_out,dif_yes_in)
                  gc()

                  if (any(is.na(S))) {
                      print(paste("Error: some NA in matrix S",
                                  which(is.na(S), arr.ind = T)))
                      break()
                  }

                  logThis(paste0("Analysis of cluster: '", cl, "' - calculating p-values"), logLevel = 3)

                  p_value <- as.data.frame(pchisq(as.matrix(S), df<-1, lower.tail<-F))

                  gc()

                  coex <- ((as.matrix(yes_in) - as.matrix(estimator_yes_in))/as.matrix(new_estimator_yes_in)) +
                      ((as.matrix(no_out) - as.matrix(estimator_no_out))/as.matrix(new_estimator_no_out)) -
                      ((as.matrix(yes_out) - as.matrix(estimator_yes_out))/as.matrix(new_estimator_yes_out)) -
                      ((as.matrix(no_in) - as.matrix(estimator_no_in))/as.matrix(new_estimator_no_in))

                  rm(yes_in, yes_out,no_out,no_in)
                  rm(estimator_yes_out,estimator_yes_in, estimator_no_out, estimator_no_in)
                  gc()

                  coex <- coex / sqrt(1/new_estimator_yes_in + 1/new_estimator_no_in + 1/new_estimator_yes_out + 1/new_estimator_no_out)

                  coex <- coex/sqrt(getNumCells(objCOTAN))
                  coex <- as.data.frame(coex)

                  clusters_coex <- setColumnInDF(clusters_coex, colToSet = coex,
                                                 colName = cl, rowNames = rownames(zeroOne))
                  clusters_pval <- setColumnInDF(clusters_pval, colToSet = p_value,
                                                 colName = cl, rowNames = rownames(zeroOne))

    logThis(paste0("Analysis of cluster: '", cl, "' - DONE"), logLevel = 3)
  }

  logThis("Differential Expression Analysis - DONE", logLevel = 2)

  return(list(clusters_coex, clusters_pval))
}

