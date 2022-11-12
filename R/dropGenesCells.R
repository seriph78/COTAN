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
#' @importFrom rlang is_missing
#' @importFrom rlang missing_arg
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

    genes.pos.to.keep <- which(!(getGenes(objCOTAN) %in% genes))
    if (length(genes.pos.to.keep) == getNumGenes((objCOTAN))) {
      # no genes to drop
      # genes.pos.to.keep <- missing_arg()
    }

    cells.pos.to.keep <- which(!(getCells(objCOTAN) %in% cells))
    if (length(cells.pos.to.keep) == getNumCells((objCOTAN))) {
      # no cells to drop
      # cells.pos.to.keep <- missing_arg()
    }

    if (is_missing(genes.pos.to.keep) && is_missing(cells.pos.to.keep)) {
      # nothing to do
      return(objCOTAN)
    }

    objCOTAN@raw <- objCOTAN@raw[genes.pos.to.keep, cells.pos.to.keep]

    if(!is_empty(objCOTAN@rawNorm)) {
      objCOTAN@rawNorm <- objCOTAN@rawNorm[genes.pos.to.keep, cells.pos.to.keep]
    }

    if (!is_empty(objCOTAN@coex) || !is_empty(objCOTAN@cellsCoex)) {
      stop("Cannot drop genes/cells once 'coex' matrices have been initialised")  
    }
    
    if (!is_empty(objCOTAN@nu) && !is_missing(cells.pos.to.keep)) {
      objCOTAN@nu <- objCOTAN@nu[cells.pos.to.keep]
    }

    if (!is_empty(objCOTAN@lambda) && !is_missing(genes.pos.to.keep)) {
      objCOTAN@lambda <- objCOTAN@lambda[genes.pos.to.keep]
    }

    if (!is_empty(objCOTAN@dispersion) || !is_empty(objCOTAN@hkGenes)) {
      stop("Cannot drop genes/cells once 'dispersion' or 'hkGenes' have been initialised")  
    }

    if (!is_empty(objCOTAN@metaCells) && !is_missing(cells.pos.to.keep)) {
      colNames <- colnames(objCOTAN@metaCells)
      objCOTAN@metaCells <- as.data.frame(objCOTAN@metaCells[cells.pos.to.keep,],
                                          row.names = getCells(objCOTAN))
      colnames(objCOTAN@metaCells) <- colNames
    }

    if (!is_empty(objCOTAN@clustersCoex)) {
      stop("Cannot drop genes/cells once 'clustersCoex' has been initialised")  
    }

    return(objCOTAN)
  }
)
