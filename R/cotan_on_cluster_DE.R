#' DEA_on_clusters
#'
#' Next method is used to do the differential expression analysis using
#' COTAN contingency tables on the cluster checked with the previous method
#'
#' @param obj A COTAN object
#' @param cell.clusterization a vector of cell clusterization. If empty, it takes the last clusterization in obj
#'
#' @return a list with two objects: the first is a scCOTAN with the new
#' object with also the correlation matrix for the genes in each cluster (obj@cluster_data),
#' the second is the p-values matrix.
#'
#' @export
#'
#' @examples
setGeneric("DEA_on_clusters", function(obj, cell.clusterization = NULL)
  standardGeneric("DEA_on_clusters"))
#' @rdname DEA_on_clusters
setMethod("DEA_on_clusters","COTAN",
          function(obj, cell.clusterization = NULL) {

            if (is.null(cell.clusterization)) {
              #cl.name <- getClusterizations(obj)[length(getClusterizations(obj))]
              clusters.names <- unique(getClusterizationData(obj)[[1]])
              cell.clusterization <- vector("list", length(clusters.names))
              names(cell.clusterization) <- clusters.names
              for (c in clusters.names) {
                tmp <- getClusterizationData(obj)[[1]]
                tmp <- tmp[tmp %in% c]
                cell.clusterization[[as.character(c)]] <- names(tmp)
              }
            }else{
              clusters.names <- unique(cell.clusterization)
              cell.clusterization.list <- vector("list", length(clusters.names))
              names(cell.clusterization.list) <- clusters.names
              for (c in clusters.names) {
                tmp <- cell.clusterization
                tmp <- tmp[tmp %in% c]
                cell.clusterization.list[[as.character(c)]] <- names(tmp)
              }
              cell.clusterization <- cell.clusterization.list
            }

              cells <- getRawData(obj)

              mu_estimator <- calculateMu(obj)

              noHKFlags <- flagNotHousekeepingGenes(obj)

              mu_estimator <- mu_estimator[noHKFlags,]
              cells <- sign(cells[noHKFlags, ]) # now a zero-one matrix
              M <- funProbZero(getDispersion(obj), mu_estimator[,colnames(cells)]) # matrix of 0 probabilities
              N <- 1-M

              cluster_data <- data.frame()
              cluster_pval <- data.frame()
              # cells_set  > set of cell code corresponding to cluster
              for (condition in names(cell.clusterization)) {
                gc()
                  print(paste("cluster", condition))

                  cells_set <- unlist(cell.clusterization[condition])

                  cells_in <- colnames(cells) %in% cells_set

                  stopifnot("ERROR. Some cells are not present!"<-  all(cells_set %in% colnames(cells)) )

                  # Cells matrix : formed by row data matrix changed to 0-1 matrix


                  yes_in  <- rowSums(cells[, cells_in])
                  yes_out <- rowSums(cells[,!cells_in])
                  no_in   <- sum( cells_in) - yes_in
                  no_out  <- sum(!cells_in) - yes_out

                  estimator_no_in   <- rowSums(M[, cells_in])
                  estimator_no_out  <- rowSums(M[,!cells_in])

                  estimator_yes_in  <- rowSums(N[, cells_in])
                  estimator_yes_out <- rowSums(N[,!cells_in])

                  if( (sum(((yes_in+no_in)-(estimator_no_in+estimator_yes_in))**2)+
                       sum(((yes_out+no_out)-(estimator_yes_out+estimator_no_out))**2)) > 0.0001){
                      print("Problems with estimators!")
                  }


                  exp_yes <- rowSums(cells)
                  exp_no <- ncol(cells) - exp_yes
                  if( sum(yes_in+yes_out-exp_yes) != 0 |  sum(no_in+no_out-exp_no) != 0 ){
                      print("Problems with observed counts!")
                  }


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


                  print("Calculating p values")
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

                  coex <- coex/sqrt(getNumCells(obj))
                  coex <- as.data.frame(coex)

                  cluster_data <- setColumnInDF(cluster_data, colToSet = coex,
                                                colName = condition, rowNames = getGenes(obj))
                  cluster_pval <- setColumnInDF(cluster_pval, colToSet = p_value,
                                                colName = condition, rowNames = getGenes(obj))
              }

              return(list(cluster_data,cluster_pval))
          }
)
