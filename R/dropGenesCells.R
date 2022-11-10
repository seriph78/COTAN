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
    if (!is_empty(objCOTAN@rawNorm)) {
      stop("dropGenesCells is not supported after initialisation, yet")
    }
  
    genes.pos.to.keep <- which(!(getGenes(objCOTAN) %in% genes))
    cells.pos.to.keep <- which(!(getCells(objCOTAN) %in% cells))
    
    objCOTAN@raw <- objCOTAN@raw[genes.pos.to.keep, cells.pos.to.keep]

    return(objCOTAN)
  }
)
