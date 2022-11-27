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


#' addElementToMetaDataset
#'
#' This function is used to add a line of information to the information data frame (metadata).
#'
#' @param objCOTAN a COTAN object
#' @param tag the new information tag
#' @param value a value (or an array) containing the information
#'
#' @return the updated COTAN object
#'
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' ERCC.cotan <- addElementToMetaDataset(ERCC.cotan, "Test", c("These are ", "some values"))
#' getMetadataDataset(ERCC.cotan)
#' @rdname addElementToMetaDataset
setMethod(
  "addElementToMetaDataset",
  "COTAN",
  function(objCOTAN, tag, value) {
    newLine <- c(tag, value)

    objCOTAN@metaDataset[(nrow(objCOTAN@metaDataset) + 1), seq_along(newLine)] <- newLine

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
    # determine positive UMI
    cells <- getZeroOneProj(objCOTAN)

    # flag the genes with positive UMI count in every single cell
    objCOTAN@metaGenes <- setColumnInDF(objCOTAN@metaGenes,
                                        rowSums(cells) == getNumCells(objCOTAN),
                                        "hkGenes", getGenes(objCOTAN))

    return(objCOTAN)
  }
)


#' dropGenesCells
#'
#' This function remove an array of genes and/or cells from the current COTAN
#' object.
#'
#' @param objCOTAN a COTAN object
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
#' @param objCOTAN a COTAN object
#'
#' @return the original object but with 'genesCoex' and 'cellsCoex' slots
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

    if (!is_empty(objCOTAN@genesCoex) && isa(objCOTAN@genesCoex, "dtCMatrix")) {
      print("'genesCoex' slot in the old format! Converting...")
      objCOTAN@genesCoex <- mat2vec_rfast(objCOTAN@genesCoex)
    }

    if (!is_empty(objCOTAN@cellsCoex) && isa(objCOTAN@cellsCoex, "dtCMatrix")) {
      print("'cellsCoex' in the old format! Converting...")
      objCOTAN@cellsCoex <- mat2vec_rfast(objCOTAN@cellsCoex)
    }

    return(objCOTAN)
  }
)
