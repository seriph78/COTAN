setMethod(
  "plotCellsHeatmap",
  "COTAN",
  function(objCOTAN, cellsName, clusters) {
    # controls
    if (is_empty(objCOTAN@cellsCoex)) {
      stop("cellsCoex filed of the COTAN object is empty")
    }

    # if(xor(missing(cellsName), missing(clusters))){
    #   stop("exactly one between the parameters 'cellsName' and 'clusters' must not be empty")
    # }

    # coex matrix preparation
    coexMat <- vec2mat_rfast(objCOTAN@cellsCoex)
    diag(coexMat) <- integer(length = length(diag(coexMat)))

    # if clustering is needed
    if (!missing(clusters)) {
      # identifier for each cluster
      clustersIdentifier <- unique(sort(clusters))

      # size of each cluster
      clustersSize <- integer(length(clustersIdentifier))
      for (i in seq(1, length(clustersIdentifier))) {
        clustersSize[i] <- sum(clusters == clustersIdentifier[i])
      }

      # cell names ordered by the identifier of the cluster to which they belong
      nameSort <- names(sort(clusters))
      coexMat <- coexMat[nameSort, nameSort]

      colnames(coexMat) <- clusters[nameSort]
      rownames(coexMat) <- clusters[nameSort]

      Heatmap(coexMat,
        border = TRUE,
        column_split = factor(rep(clustersIdentifier, clustersSize), 
                              levels = title
        ),
        cluster_rows = FALSE,
        cluster_columns = FALSE
      )
    } else {
      Heatmap(coexMat[cellsName, cellsName],
        cluster_rows = FALSE,
        cluster_columns = FALSE
      )
    }
  }
)
