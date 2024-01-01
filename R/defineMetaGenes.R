#' Meta-genes: definition and usage
#'
#' @description A *meta-gene* is a small group of genes deemed similar according
#'   to some given distance. This grouping is intended to make less noisy and
#'   more robust any procedure that would try to cluster cells based on coex
#'   information
#'
#' @details `defineMetaGenes()` groups together *similar* genes into
#'   *meta-genes* according to the given distance. Usually the distance is based
#'   on some function of the PCA of the genes
#'
#' @param points a `matrix` where each point is presented as a row
#' @param maxDist a maximum distance between two points inside a single
#'   *meta-gene*
#' @param minSize: the minimum size of a *meta-gene*
#' @param maxSize: the maximum size of a *meta-gene*
#' @param distance the distance method to use. Default is `"euclidian"`
#' @param permute: a Boolean. When `FALSE` the points are used in the given
#'   order as potential meta-genes centers, otherwise a random a permutation is
#'   performed
#'
#' @return `defineMetaGenes()`returns a `list` of *meta-genes*, each named after
#'   its center gene. A *meta-gene* is an array of gene names, thus making this
#'   similar to a [ClusterList]
#'
#' @export
#'
#' @importFrom parallelDist parDist
#'
#' @rdname MetaGenes
#'
defineMetaGenes <- function(points,
                            maxDist = 15.0,
                            minSize = 25L,
                            maxSize = 100L,
                            distance = "cosine",
                            permutation = FALSE) {
  logThis("Defining meta-genes - START", logLevel = 2L)

  metaGenes <- vector(mode = "list")

  possibleCenters <- rownames(points)
  if (isTRUE(permutation)) {
    possibleCenters <- sample(possibleCenters)
  }

  pointsDist <- as.matrix(parDist(points, method = distance,
                                  diag = TRUE, upper = TRUE))

  repeat{
    center <- possibleCenters[[1L]]

    logThis(paste0("Working on center: ", center, " - "),
            logLevel = 3L, appendLF = FALSE)

    distances <- sort(pointsDist[center, ], decreasing = FALSE)

    numNeighbors <- sum(distances <= maxDist)

    if (numNeighbors < minSize) {
      # this center is too isolated
      possibleCenters <- possibleCenters[2L:length(possibleCenters)]
      logThis("dropped", logLevel = 3L)
    } else {
      # the closest points are inserted into the set
      namesOfClosest  <- names(head(distances, n = min(maxSize, numNeighbors)))
      possibleCenters <- setdiff(possibleCenters, namesOfClosest)

      metaGenes[[center]] <- namesOfClosest
      logThis(paste0("found ", length(namesOfClosest), " near points"),
                     logLevel = 3L)
    }

    if (is_empty(possibleCenters)) {
      # no points left to process
      break
    }
  }

  logThis("Defining meta-genes - DONE", logLevel = 2L)

  return(metaGenes)
}
