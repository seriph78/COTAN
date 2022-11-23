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

    objCOTAN@metaDataset[1,seq_len(2)] = c("GEO:", GEO)
    objCOTAN@metaDataset[2,seq_len(2)] = c("scRNAseq method:", sequencingMethod)
    objCOTAN@metaDataset[3,seq_len(2)] = c("starting n. of cells:", getNumCells(objCOTAN))
    objCOTAN@metaDataset[4,seq_len(2)] = c("Condition sample:", sampleCondition)

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
#' @return a new object with the raw data where the with genes/cells were
#' expunged as indicated. The meta-data for the dataset are kept, while the
#' rest is dropped as no more relevant with the restricted matrix
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
      # nothing dropped
      return(objCOTAN)
    }
    else {
      # as all estimates would be wrong, a completely new object is created
      # with the same meta data for the data-set as the original
      output <- COTAN(objCOTAN@raw[genesPosToKeep, cellsPosToKeep])
      output@metaDataset <- getMetadataDataset(objCOTAN)

      return(output)
    }
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

    if (!is_empty(objCOTAN@coex) && isa(objCOTAN@coex, "dtCMatrix")) {
      print("Coex slot in the old format! Converting...")
      objCOTAN@coex <- mat2vec_rfast(objCOTAN@coex)
    }

    if (!is_empty(objCOTAN@cellsCoex) && isa(objCOTAN@cellsCoex, "dtCMatrix")) {
      print("CellsCoex in the old format! Converting...")
      objCOTAN@cellsCoex <- mat2vec_rfast(objCOTAN@cellsCoex)
    }

    return(objCOTAN)
  }
)
