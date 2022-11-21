#' initializeMetaDataset
#'
#' initialize meta-data data-set
#'
#' @param objCOTAN the COTAN object
#' @param GEO a code reporting the GEO identification
#'            or other specific dataset code
#' @param sequencingMethod a string reporting the method used for the sequencing
#' @param sampleCondition a string reporting the specific sample condition or time point
#'
#' @return the given COTAN object with updated metaDataset
#' @export
#' @examples
#'
#' data("raw.dataset")
#' obj <- COTAN(raw = raw.dataset)
#' obj <- initRaw(obj, GEO = "code", sequencingMethod = "10X",
#'                     sampleCondition = "mouse dataset")
#'
#' @rdname initializeMetaDataset
setMethod(
  "initializeMetaDataset",
  "COTAN",
  function(objCOTAN, GEO, sequencingMethod, sampleCondition) {
    print("Initializing COTAN meta-data")

    numCells = getNumCells(objCOTAN)
    objCOTAN@metaDataset[1,seq_len(2)] = c("GEO:", GEO)
    objCOTAN@metaDataset[2,seq_len(2)] = c("scRNAseq method:", sequencingMethod)
    objCOTAN@metaDataset[3,seq_len(2)] = c("starting n. of cells:", numCells)
    objCOTAN@metaDataset[4,seq_len(2)] = c("Condition sample:", sampleCondition)

    #TODO: remove this!
    objCOTAN@metaCells <- data.frame(clusters = rep(NA, numCells),
                                     row.names = getCells(objCOTAN))

    return(objCOTAN)
  }
)

#' findHousekeepingGenes
#'
#' determines the housekeeping genes vector of a COTAN object
#'
#' @param objCOTAN the COTAN object
#' @return the given COTAN object with updated housekeepingGenes
#'
#' @importFrom Matrix rowSums
#' @export
#'
#' @rdname findHousekeepingGenes
setMethod(
  "findHousekeepingGenes",
  "COTAN",
  function(objCOTAN) {
    nCells <- getNumCells(objCOTAN)

    if (nCells == 0) {
      stop("zero cells given")
    }

    # determine positive UMI
    cells <- getZeroOneProj(objCOTAN)

    # name of the genes with positive UMI count in every single cell
    objCOTAN@hkGenes <- names(which(rowSums(cells) == nCells))

    return(objCOTAN)
  }
)


#' dropGenesCells
#'
#' This function remove an array of genes and/or cells from the current COTAN
#' object.
#'
#' @param object a COTAN object
#' @param genes an array of gene names
#' @param cells an array of cell names
#'
#' @return the original object but with genes/cells expunged as indicated.
#' @importFrom rlang is_empty
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' genes.to.rem <- getGenes(ERCC.cotan)[grep('^MT', getGenes(ERCC.cotan))]
#' cells.to.rem <- getCells(ERCC.cotan)[which(getCellsSize(ERCC.cotan) == 0)]
#' ERCC.cotan <- dropGenesCells(ERCC.cotan, genes.to.rem, cells.to.rem)
#' @rdname dropGenesCells
setMethod(
  "dropGenesCells",
  "COTAN",
  function(objCOTAN, genes, cells) {
    stopifnot(validObject(objCOTAN))

    if (is_empty(objCOTAN@raw)) {
      stop("No raw count data was given")
    }

    genesPosToKeep <- which(!(getGenes(objCOTAN) %in% genes))
    cellsPosToKeep <- which(!(getCells(objCOTAN) %in% cells))

    anyGenesDropped <- length(genesPosToKeep) != getNumGenes((objCOTAN))
    anyCellsDropped <- length(cellsPosToKeep) != getNumCells((objCOTAN))

    if (!anyGenesDropped && !anyCellsDropped) {
      # nothing to do
      return(objCOTAN)
    }

    objCOTAN@raw <- objCOTAN@raw[genesPosToKeep, cellsPosToKeep]

    if (!is_empty(objCOTAN@coex) || !is_empty(objCOTAN@cellsCoex)) {
      stop("Cannot drop genes/cells once 'coex' matrices have been initialised")
    }

    if (!is_empty(objCOTAN@nu) && anyCellsDropped) {
      objCOTAN@nu <- objCOTAN@nu[cellsPosToKeep]
    }

    if (!is_empty(objCOTAN@lambda) && anyGenesDropped) {
      objCOTAN@lambda <- objCOTAN@lambda[genesPosToKeep]
    }

    if (!is_empty(objCOTAN@dispersion) || !is_empty(objCOTAN@hkGenes)) {
      stop("Cannot drop genes/cells once 'dispersion' or 'hkGenes' have been initialised")
    }

    if (!is_empty(objCOTAN@metaCells) && anyCellsDropped) {
      colNames <- colnames(objCOTAN@metaCells)
      objCOTAN@metaCells <- as.data.frame(objCOTAN@metaCells[cellsPosToKeep,],
                                          row.names = getCells(objCOTAN))
      colnames(objCOTAN@metaCells) <- colNames
    }

    if (!is_empty(objCOTAN@clustersCoex)) {
      stop("Cannot drop genes/cells once 'clustersCoex' has been initialised")
    }

    validObject(objCOTAN)

    return(objCOTAN)
  }
)


#' standardizeCoex
#'
#' @param object a COTAN object
#'
#' @return the original object but with coex and cellsCoex slots
#' in standard format.
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @rdname standardizeCoex
#'
setMethod(
  "standardizeCoex",
  "COTAN",
  function(objCOTAN) {

    if (!is_empty(objCOTAN@coex) &&
        (is(class(objCOTAN@coex)[1], "dtCMatrix") ||
         as.vector(class(objCOTAN@coex)) %in% "dtCMatrix") ) {
      print("Coex slot in the old format! Converting...")
      objCOTAN@coex <- mat2vec_rfast(objCOTAN@coex)
    }

    if (!is_empty(objCOTAN@cellsCoex) &&
        (is(class(objCOTAN@cellsCoex)[1], "dtCMatrix") ||
         as.vector(class(objCOTAN@cellsCoex)) %in% "dtCMatrix") ) {
      print("CellsCoex in the old format! Converting...")
      objCOTAN@cellsCoex <- mat2vec_rfast(objCOTAN@cellsCoex)
    }

    return(objCOTAN)
  }
)
