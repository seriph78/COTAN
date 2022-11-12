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

    anyGenesDropped <- length(genesPosToKeep) == getNumGenes((objCOTAN))
    anyCellsDropped <- length(cellsPosToKeep) == getNumCells((objCOTAN))

    if (!anyGenesDropped && !anyCellsDropped) {
      # nothing to do
      return(objCOTAN)
    }

    objCOTAN@raw <- objCOTAN@raw[genesPosToKeep, cellsPosToKeep]

    if(!is_empty(objCOTAN@rawNorm)) {
      objCOTAN@rawNorm <- objCOTAN@rawNorm[genesPosToKeep, cellsPosToKeep]
    }

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

    return(objCOTAN)
  }
)
