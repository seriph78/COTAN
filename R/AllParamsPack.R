# ----------------- execution options --------------------

#' @title Execution options for heavy COTAN computations
#'
#' @description A small parameter object bundling execution-related controls
#' such as multi-core usage, `torch` device selection, and solver batching
#' options. Some fields are only used by specific functions.
#'
#' @slot cores Integer scalar. Requested CPU core count; resolved by
#'   [handleMultiCore()].
#' @slot optimizeForSpeed Logical scalar. Whether to try \pkg{torch}
#'   acceleration; resolved by [canUseTorch()].
#' @slot deviceStr Character scalar. Requested \pkg{torch} device string;
#'   resolved by [canUseTorch()].
#' @slot chunkSize Integer scalar. Solver batch size for solver paths that
#'   support batching.
#'
#' @name ExecutionOptions-class
#'
#' @exportClass ExecutionOptions
#'
#' @rdname ExecutionOptions
#'
setClass(
  "ExecutionOptions",
  slots = c(
    cores = "integer",
    optimizeForSpeed = "logical",
    deviceStr = "character",
    chunkSize = "integer"
  ),
  prototype = list(
    cores = 1L,
    optimizeForSpeed = TRUE,
    deviceStr = "cuda",
    chunkSize = 1024L
  ),
  validity = function(object) {
    if (length(object@cores) != 1L || is.na(object@cores) ||
        object@cores < 0L) {
      return("`cores` must be a non-negative integer scalar")
    }

    if (length(object@optimizeForSpeed) != 1L ||
        is.na(object@optimizeForSpeed)) {
      return("`optimizeForSpeed` must be a non-missing logical scalar")
    }

    if (length(object@deviceStr) != 1L || is.na(object@deviceStr)) {
      return("`deviceStr` must be a character scalar")
    }

    if (length(object@chunkSize) != 1L || is.na(object@chunkSize) ||
        object@chunkSize < 1L) {
      return("`chunkSize` must be a positive integer scalar")
    }

    return(TRUE)
  }
)

#' @title Build execution options
#'
#' @param cores Requested number of CPU cores. The effective value is bounded by
#'   [handleMultiCore()], using the current session capabilities reported by
#'   \pkg{parallelly}.
#' @param optimizeForSpeed Whether to try accelerated computation through
#'   \pkg{torch}. See [canUseTorch()] for the runtime checks and fallback rules.
#' @param deviceStr Requested \pkg{torch} device string, for example `"cpu"`,
#'   `"cuda"`, or `"cuda:0"`. See [canUseTorch()] for fallback behavior when
#'   the requested device is unavailable.
#' @param chunkSize Integer scalar controlling solver batching where supported,
#'   in particular dispersion / p-value solver paths such as
#'   [estimateDispersionViaSolver()].
#'
#' @returns An object of class `ExecutionOptions`
#'
#' @seealso [handleMultiCore()], [canUseTorch()], [calculateCoex()],
#'   [estimateDispersionViaSolver()]
#'
#' @export
#'
#' @examples
#'   exec <- ExecutionOptions(
#'     cores = 4L,
#'     optimizeForSpeed = TRUE,
#'     deviceStr = "cuda",
#'     chunkSize = 1024L
#'   )
#'
#' @rdname ExecutionOptions
ExecutionOptions <- function(cores = 1L,
                             optimizeForSpeed = TRUE,
                             deviceStr = "cuda",
                             chunkSize = 1024L) {
  methods::new(
    "ExecutionOptions",
    cores = as.integer(cores),
    optimizeForSpeed = as.logical(optimizeForSpeed),
    deviceStr = as.character(deviceStr),
    chunkSize = as.integer(chunkSize)
  )
}


# internal helper: build the object from legacy loose parameters
legacyExecutionOptions <- function(cores = 1L,
                                   optimizeForSpeed = TRUE,
                                   deviceStr = "cuda",
                                   chunkSize = 1024L) {
  ExecutionOptions(
    cores = cores,
    optimizeForSpeed = optimizeForSpeed,
    deviceStr = deviceStr,
    chunkSize = chunkSize
  )
}


# internal helper: apply runtime feasibility checks
resolveExecutionOptions <- function(executionOptions) {
  assertthat::assert_that(
    methods::is(executionOptions, "ExecutionOptions"),
    msg = "`executionOptions` must be an `ExecutionOptions` object"
  )

  effCores <- handleMultiCore(executionOptions@cores)

  torchInfo <- canUseTorch(
    optimizeForSpeed = executionOptions@optimizeForSpeed,
    deviceStr = executionOptions@deviceStr
  )

  return(list(
    "cores" = effCores,
    "useTorch" = torchInfo[["useTorch"]],
    "deviceStr" = torchInfo[["deviceStr"]],
    "chunkSize" = executionOptions@chunkSize
  ))
}


# ----------------- reduction options --------------------

#' @title Reduction options for COTAN dimensionality reduction
#'
#' @description A small parameter object bundling controls used to build the
#'   reduced data matrix used by clusterizations and UMAP plots.
#'
#'   This object intentionally stores only dimensionality-reduction policy
#'   parameters. UMAP layout options and clustering-specific options should stay
#'   outside this class.
#'
#' @slot useCoexEigen Logical scalar. Whether to use the first COEX
#'   eigen-vectors instead of PCA on a selected gene matrix.
#' @slot dataMethod Character scalar. Data matrix method passed to
#'   [getDataMatrix()].
#' @slot numComp Integer scalar. Number of reduced components to calculate.
#' @slot genesSel Character vector. Selector name or explicit gene list passed
#'   to [getSelectedGenes()]. Empty string is usually allowed as default.
#' @slot numGenes Integer scalar. Number of genes to select when `genesSel`
#'   names a selection method.
#'
#' @name ReductionOptions-class
#'
#' @exportClass ReductionOptions
#'
#' @rdname ReductionOptions
#'
setClass(
  "ReductionOptions",
  slots = c(
    useCoexEigen = "logical",
    dataMethod = "character",
    numComp = "integer",
    genesSel = "character",
    numGenes = "integer"
  ),
  prototype = list(
    useCoexEigen = TRUE,
    dataMethod = "LogLikelihood",
    numComp = 25L,
    genesSel = "HGDI",
    numGenes = 2000L
  ),
  validity = function(object) {
    if (length(object@useCoexEigen) != 1L ||
        is.na(object@useCoexEigen)) {
      return("`useCoexEigen` must be a non-missing logical scalar")
    }

    if (length(object@dataMethod) != 1L || is.na(object@dataMethod)) {
      return("`dataMethod` must be a non-missing character scalar")
    }

    if (length(object@numComp) != 1L || is.na(object@numComp) ||
        object@numComp < 1L) {
      return("`numComp` must be a positive integer scalar")
    }

    if (length(object@genesSel) < 1L || anyNA(object@genesSel)) {
      return("`genesSel` must be a non-missing character vector")
    }

    if (length(object@genesSel) > 1L &&
        any(sapply(object@genesSel, isEmptyName))) {
      return(paste(
        "`genesSel` can be an explicit vector of gene names,",
        "but in that case all names must be non-empty strings"
      ))
    }

    usesGeneSelector <-
      length(object@genesSel) == 1L && !isEmptyName(object@genesSel)

    if (length(object@numGenes) != 1L || is.na(object@numGenes) ||
        (isFALSE(object@useCoexEigen) && usesGeneSelector &&
         object@numGenes < 1L)) {
      return(paste("`numGenes` must be a positive integer scalar",
                   "when a gene selector is used"))
    }

    return(TRUE)
  }
)

#' @title Build reduction options
#'
#' @param useCoexEigen Whether to use the first COEX eigenvectors instead of PCA
#'   on a selected gene matrix. See [calculateReducedDataMatrix()] for the exact
#'   reduction path.
#' @param dataMethod Data matrix method passed to [getDataMatrix()]. See
#'   [getDataMatrix()] for accepted aliases such as `"LogNormalized"` and
#'   `"LogLikelihood"`.
#' @param numComp Number of reduced components to calculate. See
#'   [calculateReducedDataMatrix()] for how this is interpreted by the
#'   COEX-eigen and PCA branches.
#' @param genesSel Gene-selection method or explicit gene vector passed to
#'   [getSelectedGenes()]. See [getSelectedGenes()] for accepted selector names
#'   such as `"HGDI"`, `"HVG_Seurat"`, and `"HVG_Scanpy"`.
#' @param numGenes Number of genes to select when `genesSel` names a selector
#'   method. Ignored when `genesSel` is an explicit vector of gene names.
#'
#' @returns An object of class `ReductionOptions`
#'
#' @seealso [calculateReducedDataMatrix()], [getDataMatrix()],
#'   [getSelectedGenes()], [cellsUMAPPlot()], [cellsUniformClustering()]
#'
#' @export
#'
#' @examples
#'   redOpt <- ReductionOptions()
#'
#'   pcaRedOpt <- ReductionOptions(
#'     useCoexEigen = FALSE,
#'     dataMethod = "LogNormalized",
#'     numComp = 25L,
#'     genesSel = "HGDI",
#'     numGenes = 2000L
#'   )
#'
#' @rdname ReductionOptions
#'
ReductionOptions <- function(useCoexEigen = TRUE,
                             dataMethod = "LogLikelihood",
                             numComp = 25L,
                             genesSel = "HGDI",
                             numGenes = 2000L) {
  methods::new(
    "ReductionOptions",
    useCoexEigen = as.logical(useCoexEigen),
    dataMethod = as.character(dataMethod),
    numComp = as.integer(numComp),
    genesSel = as.character(genesSel),
    numGenes = as.integer(numGenes)
  )
}


# internal helper: build the object from legacy loose parameters
legacyReductionOptions <- function(useCoexEigen,
                                   dataMethod,
                                   numComp,
                                   genesSel,
                                   numGenes) {
  ReductionOptions(
    useCoexEigen = useCoexEigen,
    dataMethod = dataMethod,
    numComp = numComp,
    genesSel = genesSel,
    numGenes = numGenes
  )
}




# ----------------- cleaning options --------------------

#' @title Cleaning options
#'
#' @description Parameter object bundling the thresholds used by [clean()]:
#'   low-expression cutoffs used to drop genes/cells and high-expression
#'   thresholds used to mark fully-expressed genes or fully-expressing cells.
#'
#' @slot cellsCutoff Numeric scalar. Genes expressed in at most this fraction
#'   of cells are dropped.
#' @slot genesCutoff Numeric scalar. Cells expressing at most this fraction of
#'   genes are dropped.
#' @slot cellsThreshold Numeric scalar. Genes expressed in more than this
#'   fraction of cells are marked as fully-expressed.
#' @slot genesThreshold Numeric scalar. Cells expressing more than this fraction
#'   of genes are marked as fully-expressing.
#'
#' @name CleaningOptions-class
#'
#' @exportClass CleaningOptions
#'
#' @rdname CleaningOptions
#'
setClass(
  "CleaningOptions",
  slots = c(
    cellsCutoff = "numeric",
    genesCutoff = "numeric",
    cellsThreshold = "numeric",
    genesThreshold = "numeric"
  ),
  prototype = list(
    cellsCutoff = 0.003,
    genesCutoff = 0.002,
    cellsThreshold = 0.99,
    genesThreshold = 0.99
  ),
  validity = function(object) {
    numericSlots <- c("cellsCutoff", "genesCutoff",
                      "cellsThreshold", "genesThreshold")

    for (slotName in numericSlots) {
      value <- methods::slot(object, slotName)
      if (length(value) != 1L || is.na(value) || value < 0.0) {
        return(paste0("`", slotName,
                      "` must be a non-missing non-negative numeric scalar"))
      }
    }

    return(TRUE)
  }
)

#' @title Build cleaning options
#'
#' @param cellsCutoff Fraction of cells used as the low-expression cutoff for
#'   genes. Consumed by [clean()].
#' @param genesCutoff Fraction of genes used as the low-expression cutoff for
#'   cells. Consumed by [clean()].
#' @param cellsThreshold Fraction of cells used to mark fully-expressed genes.
#'   Consumed by [clean()] and related fully-expressed-gene utilities.
#' @param genesThreshold Fraction of genes used to mark fully-expressing cells.
#'   Consumed by [clean()] and related fully-expressing-cell utilities.
#'
#' @returns An object of class `CleaningOptions`
#'
#' @seealso [clean()], [findFullyExpressedGenes()], [findFullyExpressingCells()]
#'
#' @export
#'
#' @examples
#'   cleanOpt <- CleaningOptions()
#'
#'   stricterCleanOpt <- CleaningOptions(
#'     cellsCutoff = 0.005,
#'     genesCutoff = 0.003,
#'     cellsThreshold = 0.98,
#'     genesThreshold = 0.98
#'   )
#'
#' @rdname CleaningOptions
#'
CleaningOptions <- function(cellsCutoff = 0.003,
                            genesCutoff = 0.002,
                            cellsThreshold = 0.99,
                            genesThreshold = 0.99) {
  methods::new(
    "CleaningOptions",
    cellsCutoff = as.numeric(cellsCutoff),
    genesCutoff = as.numeric(genesCutoff),
    cellsThreshold = as.numeric(cellsThreshold),
    genesThreshold = as.numeric(genesThreshold)
  )
}

legacyCleaningOptions <- function(cellsCutoff = 0.003,
                                  genesCutoff = 0.002,
                                  cellsThreshold = 0.99,
                                  genesThreshold = 0.99) {
  CleaningOptions(
    cellsCutoff = cellsCutoff,
    genesCutoff = genesCutoff,
    cellsThreshold = cellsThreshold,
    genesThreshold = genesThreshold
  )
}

resolveCleaningOptions <- function(cellsCutoff = 0.003,
                                   genesCutoff = 0.002,
                                   cellsThreshold = 0.99,
                                   genesThreshold = 0.99,
                                   cleaningOptions = NULL) {
  if (is.null(cleaningOptions)) {
    return(legacyCleaningOptions(
      cellsCutoff = cellsCutoff,
      genesCutoff = genesCutoff,
      cellsThreshold = cellsThreshold,
      genesThreshold = genesThreshold
    ))
  }

  assertthat::assert_that(
    methods::is(cleaningOptions, "CleaningOptions"),
    msg = "`cleaningOptions` must be a `CleaningOptions` object"
  )

  assertthat::assert_that(
    identical(cellsCutoff, 0.003),
    identical(genesCutoff, 0.002),
    identical(cellsThreshold, 0.99),
    identical(genesThreshold, 0.99),
    msg = paste(
      "Do not mix `cleaningOptions` with the legacy cleaning arguments",
      "`cellsCutoff`, `genesCutoff`, `cellsThreshold`, and",
      "`genesThreshold`."
    )
  )

  return(cleaningOptions)
}


# ----------------- cluster distance options --------------------

#' @title Cluster distance options
#'
#' @description Parameter object bundling the policy used to calculate
#'   distances between cell clusters.
#'
#' @slot useDEA Logical scalar. Whether to use DEA profiles to calculate
#'   distances between clusters. When `FALSE`, distances are calculated from
#'   average Zero-One counts.
#' @slot distance Character scalar. Distance method passed to
#'   [parallelDist::parDist()], or `""` for the COTAN default selected by
#'   [distancesBetweenClusters()].
#'
#' @name ClusterDistanceOptions-class
#'
#' @exportClass ClusterDistanceOptions
#'
#' @rdname ClusterDistanceOptions
#'
setClass(
  "ClusterDistanceOptions",
  slots = c(
    useDEA = "logical",
    distance = "character"
  ),
  prototype = list(
    useDEA = TRUE,
    distance = ""
  ),
  validity = function(object) {
    if (length(object@useDEA) != 1L || is.na(object@useDEA)) {
      return("`useDEA` must be a non-missing logical scalar")
    }

    if (length(object@distance) != 1L || is.na(object@distance)) {
      return("`distance` must be a non-missing character scalar")
    }

    return(TRUE)
  }
)

#' @title Build cluster distance options
#'
#' @param useDEA Whether to calculate distances from DEA profiles. When `FALSE`,
#'   distances are calculated from average Zero-One counts. See
#'   [distancesBetweenClusters()].
#' @param distance Distance method passed to [parallelDist::parDist()]. Use the
#'   empty string `""` to keep the **COTAN** function-level default: `"cosine"`
#'   when `useDEA = TRUE`, `"euclidean"` when `useDEA = FALSE`.
#'
#' @returns An object of class `ClusterDistanceOptions`
#'
#' @seealso [distancesBetweenClusters()], [parallelDist::parDist()]
#'
#' @export
#'
#' @examples
#'   clDistOpt <- ClusterDistanceOptions()
#'
#'   zeroOneDistOpt <- ClusterDistanceOptions(
#'     useDEA = FALSE,
#'     distance = "euclidean"
#'   )
#'
#' @rdname ClusterDistanceOptions
#'
ClusterDistanceOptions <- function(useDEA = TRUE,
                                   distance = "") {
  if (is.null(distance)) {
    distance <- ""
  }

  methods::new(
    "ClusterDistanceOptions",
    useDEA = as.logical(useDEA),
    distance = as.character(distance)
  )
}

legacyClusterDistanceOptions <- function(useDEA = TRUE,
                                         distance = NULL) {
  if (is.null(distance)) {
    distance <- ""
  }

  ClusterDistanceOptions(
    useDEA = useDEA,
    distance = distance
  )
}

resolveClusterDistanceOptions <- function(useDEA = TRUE,
                                          distance = NULL,
                                          clusterDistanceOptions = NULL) {
  if (is.null(clusterDistanceOptions)) {
    return(legacyClusterDistanceOptions(
      useDEA = useDEA,
      distance = distance
    ))
  }

  assertthat::assert_that(
    methods::is(clusterDistanceOptions, "ClusterDistanceOptions"),
    msg = "`clusterDistanceOptions` must be a `ClusterDistanceOptions` object"
  )

  assertthat::assert_that(
    identical(useDEA, TRUE),
    is.null(distance),
    msg = paste(
      "Do not mix `clusterDistanceOptions` with the legacy distance arguments",
      "`useDEA` and `distance`."
    )
  )

  return(clusterDistanceOptions)
}


# ----------------- cluster tree options --------------------

#' @title Cluster tree options
#'
#' @description Parameter object bundling cell-cluster distance policy together
#'   with hierarchical-tree construction policy. It extends
#'   `ClusterDistanceOptions`, so it also carries `useDEA` and `distance`.
#'
#' @details Inherited slots `useDEA` and `distance` have the same meaning as in
#'   [ClusterDistanceOptions()].
#'
#' @slot hclustMethod Character scalar. Clustering method passed to
#'   [stats::hclust()].
#'
#' @name ClusterTreeOptions-class
#'
#' @exportClass ClusterTreeOptions
#'
#' @rdname ClusterTreeOptions
#'
setClass(
  "ClusterTreeOptions",
  contains = "ClusterDistanceOptions",
  slots = c(
    hclustMethod = "character"
  ),
  prototype = list(
    hclustMethod = "ward.D2"
  ),
  validity = function(object) {
    if (length(object@hclustMethod) != 1L ||
        is.na(object@hclustMethod) ||
        isEmptyName(object@hclustMethod)) {
      return("`hclustMethod` must be a non-empty character scalar")
    }

    return(TRUE)
  }
)

#' @title Build cluster tree options
#'
#' @param useDEA Whether to calculate distances from DEA profiles. Inherited
#'   from `ClusterDistanceOptions`; see [distancesBetweenClusters()].
#' @param distance Distance method passed to [parallelDist::parDist()]. Use `""`
#'   to keep the COTAN function-level default. Inherited from
#'   `ClusterDistanceOptions`; see [distancesBetweenClusters()].
#' @param hclustMethod Clustering method passed to [stats::hclust()]. See
#'   [stats::hclust()] for accepted values such as `"ward.D2"`, `"complete"`,
#'   `"average"`, and `"single"`.
#'
#' @returns An object of class `ClusterTreeOptions`
#'
#' @seealso [ClusterDistanceOptions()], [distancesBetweenClusters()],
#'   [reorderClusterization()], [clustersTreePlot()], [stats::hclust()]
#'
#' @export
#'
#' @examples
#'   clTreeOpt <- ClusterTreeOptions()
#'
#'   zeroOneTreeOpt <- ClusterTreeOptions(
#'     useDEA = FALSE,
#'     distance = "euclidean",
#'     hclustMethod = "ward.D2"
#'   )
#'
#' @rdname ClusterTreeOptions
#'
ClusterTreeOptions <- function(useDEA = TRUE,
                               distance = "",
                               hclustMethod = "ward.D2") {
  if (is.null(distance)) {
    distance <- ""
  }

  methods::new(
    "ClusterTreeOptions",
    useDEA = as.logical(useDEA),
    distance = as.character(distance),
    hclustMethod = as.character(hclustMethod)
  )
}

legacyClusterTreeOptions <- function(useDEA = TRUE,
                                     distance = NULL,
                                     hclustMethod = "ward.D2") {
  ClusterTreeOptions(
    useDEA = useDEA,
    distance = distance,
    hclustMethod = hclustMethod
  )
}

resolveClusterTreeOptions <- function(useDEA = TRUE,
                                      distance = NULL,
                                      hclustMethod = "ward.D2",
                                      clusterTreeOptions = NULL) {
  if (is.null(clusterTreeOptions)) {
    return(legacyClusterTreeOptions(
      useDEA = useDEA,
      distance = distance,
      hclustMethod = hclustMethod
    ))
  }

  assertthat::assert_that(
    methods::is(clusterTreeOptions, "ClusterTreeOptions"),
    msg = "`clusterTreeOptions` must be a `ClusterTreeOptions` object"
  )

  assertthat::assert_that(
    identical(useDEA, TRUE),
    is.null(distance),
    identical(hclustMethod, "ward.D2"),
    msg = paste(
      "Do not mix `clusterTreeOptions` with the legacy tree arguments",
      "`useDEA`, `distance`, and `hclustMethod`."
    )
  )

  return(clusterTreeOptions)
}
