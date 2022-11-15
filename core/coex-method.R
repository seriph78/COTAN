#' This function estimates and stores the coex matrix in the coex field if the
#' parameter 'cells' is FALSE, otherwise it estimates and stores the cellsCoex
#' matrix in the cellsCoex field
#' @param objCOTAN A COTAN object
#' @param cells Boolean, if true, the function works for the cells,
#' otherwise for the genes
#'
#' @return It returns a COTAN object
#' @export
setMethod(
  "coex",
  "COTAN",
  function(objCOTAN, cells) {
    objCOTAN <- findHousekeepingGenes(objCOTAN)

    # observed yes/yes contingency table
    observedYY <- observedContingencyYY(objCOTAN, cells)

    if(!cells){
      # remove the genes in housekeepingGenes
      observedYY <- observedYY[
        !rownames(observedYY) %in% objCOTAN@hkGenes,
        !colnames(observedYY) %in% objCOTAN@hkGenes
      ]
    }

    observedYY <- mat2vec_rfast(as.matrix(observedYY))

    # four estimators: expectedYY, expectedYN, expectedNY, expectedNN
    expectedTable <- expectedContingencyTables(objCOTAN, cells)

    expectedYY <- mat2vec_rfast(expectedTable$expectedYY)
    expectedYY$values[expectedYY$values < 1] <- 1

    expectedYN <- mat2vec_rfast(expectedTable$expectedYN)
    expectedYN$values[expectedYN$values < 1] <- 1

    expectedNY <- mat2vec_rfast(expectedTable$expectedNY)
    expectedNY$values[expectedNY$values < 1] <- 1

    expectedNN <- mat2vec_rfast(expectedTable$expectedNN)
    expectedNN$values[expectedNN$values < 1] <- 1

    # sum for division
    sumForDiv <-
      1 / expectedYY$values +
      1 / expectedYN$values +
      1 / expectedNY$values +
      1 / expectedNN$values

    # coex estimation
    coex <- (observedYY$values -
      mat2vec_rfast(as.matrix(expectedTable$expectedYY))$values)
    coex <- coex * sqrt(sumForDiv)

    if (cells){
      coex <- coex / sqrt(nrow(objCOTAN@raw)) # divided by radq(n)
    } else {
      coex <- coex / sqrt(getNumCells(objCOTAN)) # divided by radq(m)
    }

    coex <- list("genes" = observedYY$genes, "values" = coex)

    if(cells){
      objCOTAN@cellsCoex <- coex
    }else{
      objCOTAN@coex <- coex
    }

    return(objCOTAN)
  }
)
