
# internal implementation

.cleanImpl <- function(objCOTAN,
                       cleaningOptions) {
  cellsCutoff <- cleaningOptions@cellsCutoff
  genesCutoff <- cleaningOptions@genesCutoff
  cellsThreshold <- cleaningOptions@cellsThreshold
  genesThreshold <- cleaningOptions@genesThreshold

  # We want to discard genes having less than given cutoff
  # default: less than 3 non-zero counts per 1000 cells
  cutoff <- round(getNumCells(objCOTAN) * cellsCutoff, digits = 0L)
  genesToDrop <-
    getGenes(objCOTAN)[getNumOfExpressingCells(objCOTAN) <= cutoff]

  if (!is_empty(genesToDrop)) {
    objCOTAN <- dropGenesCells(objCOTAN, genes = genesToDrop)
  }

  # We want to discard cells having less than given cutoff
  # default: less than 2 non-zero counts per 1000 genes
  cutoff <- round(getNumGenes(objCOTAN) * genesCutoff, digits = 0L)
  cellsToDrop <-
    getCells(objCOTAN)[getNumExpressedGenes(objCOTAN) <= cutoff]

  if (!is_empty(cellsToDrop)) {
    objCOTAN <- dropGenesCells(objCOTAN, cells = cellsToDrop)
  }

  logThis(paste0("Genes/cells selection done:",
                 " dropped [", length(genesToDrop), "] genes",
                 " and [", length(cellsToDrop), "] cells"),
          logLevel = 2L)

  logThis(paste0("Working on [", getNumGenes(objCOTAN), "]",
                 " genes and [", getNumCells(objCOTAN), "] cells"),
          logLevel = 1L)

  objCOTAN <- findFullyExpressedGenes(objCOTAN,
                                      cellsThreshold = cellsThreshold)
  objCOTAN <- findFullyExpressingCells(objCOTAN,
                                       genesThreshold = genesThreshold)

  objCOTAN <- estimateNuLinear(objCOTAN)

  return(objCOTAN)
}


#' @aliases clean
#'
#' @details `clean()` is the main method that can be used to check and clean the
#'   dataset. It will discard any genes that has less than 3 non-zero counts per
#'   thousand cells and all cells expressing less than 2 per thousand genes.
#'   also produces and stores the estimators for `nu`
#'
#' @param objCOTAN a `COTAN` object
#' @param cellsCutoff `clean()` will delete from the `raw` data any gene that is
#'   expressed in less cells than threshold times the total number of cells.
#'   Default cutoff is \eqn{0.003 \; (0.3\%)}
#' @param genesCutoff `clean()` will delete from the `raw` data any cell that is
#'   expressing less genes than threshold times the total number of genes.
#'   Default cutoff is \eqn{0.002 \; (0.2\%)}
#' @param cellsThreshold any gene that is expressed in more cells than threshold
#'   times the total number of cells will be marked as **fully-expressed**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
#' @param genesThreshold any cell that is expressing more genes than threshold
#'   times the total number of genes will be marked as **fully-expressing**.
#'   Default threshold is \eqn{0.99 \; (99.0\%)}
#' @param cleaningOptions A `CleaningOptions` object bundling cleaning cutoffs
#'   and fully-expressed / fully-expressing thresholds. This is the preferred
#'   interface for new code. It must not be mixed with the legacy cleaning
#'   arguments `cellsCutoff`, `genesCutoff`, `cellsThreshold`, and
#'   `genesThreshold`.
#'
#' @returns `clean()` returns the updated `COTAN` object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname RawDataCleaning
#'
setMethod(
  "clean",
  signature(
    objCOTAN = "COTAN",
    cleaningOptions = "missing"),
  function(objCOTAN, cellsCutoff = 0.003, genesCutoff = 0.002,
           cellsThreshold = 0.99, genesThreshold = 0.99,
           cleaningOptions = NULL) {

    cleaningOptions <- resolveCleaningOptions(
      cellsCutoff = cellsCutoff,
      genesCutoff = genesCutoff,
      cellsThreshold = cellsThreshold,
      genesThreshold = genesThreshold,
      cleaningOptions = cleaningOptions
    )

    return(.cleanImpl(objCOTAN = objCOTAN, cleaningOptions = cleaningOptions))
  }
)


#' @details Alternative interface using a `CleaningOptions` object.
#'
#' @rdname RawDataCleaning
#' @aliases clean,COTAN,CleaningOptions-method
setMethod(
  "clean",
  signature(
    objCOTAN = "COTAN",
    cleaningOptions = "CleaningOptions"),
  function(objCOTAN, cellsCutoff = 0.003, genesCutoff = 0.002,
           cellsThreshold = 0.99, genesThreshold = 0.99,
           cleaningOptions) {
    assert_that(
      identical(cellsCutoff, 0.003),
      identical(genesCutoff, 0.002),
      identical(cellsThreshold, 0.99),
      identical(genesThreshold, 0.99),
      msg = paste("Do not mix `cleaningOptions` with",
                  "legacy cleaning arguments",
                  "(`cellsCutoff`, `genesCutoff`,",
                  "`cellsThreshold`, `genesThreshold`)."))

    return(.cleanImpl(objCOTAN = objCOTAN, cleaningOptions = cleaningOptions))
  }
)


#' @details Invalid non-missing `cleaningOptions` values are rejected
#'   explicitly.
#'
#' @rdname RawDataCleaning
#' @aliases clean,COTAN,ANY-method
setMethod(
  "clean",
  signature(
    objCOTAN = "COTAN",
    cleaningOptions = "ANY"),
  function(objCOTAN, cellsCutoff = 0.003, genesCutoff = 0.002,
           cellsThreshold = 0.99, genesThreshold = 0.99,
           cleaningOptions) {
    assert_that(
      methods::is(cleaningOptions, "CleaningOptions"),
      msg = "`cleaningOptions` must be a `CleaningOptions` object"
    )
  }
)
