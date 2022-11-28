#' initializeMetaDataset
#'
#' initialize meta-data data-set
#'
#' @param objCOTAN the COTAN object
#' @param GEO a code reporting the GEO identification or other specific dataset code
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
    validObject(objCOTAN)

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


#' addClusterization
#'
#' This function adds a clusterization to the current COTAN object, by adding a
#' new column in the \code{metaCells} data.frame and adding a new element
#' in the \code{clustersCoex} list using the passed in coex data.frame.
#'
#' @param objCOTAN a COTAN object
#' @param clusterizationName the name of a new clusterization.
#' It cannot match an existing name; in case use \link{\code{dropClusterization}} first.
#' @param clusters a factors array of cluster names
#' @param coexDF a data.frame where each column indicates the coex
#' for each (or some) of the clusters of the clusterization
#'
#' @return the updated COTAN object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' ERCC.cotan <- addClusterization(ERCC.cotan, "merged", clusters)
#'
#' @rdname addClusterization
setMethod(
  "addClusterization",
  "COTAN",
  function(objCOTAN, clusterizationName, clusters, coexDF = data.frame()) {
    clName <- clusterizationName
    if (!startsWith(clName, "CL_")) {
      clName <- paste0("CL_", clName)
    }

    if (clName %in% colnames(getMetadataCells(objCOTAN))) {
      stop(paste0("A clusterization with name '", clusterizationName, "' already exists."))
    }
    if (length(clusters) != getNumCells(objCOTAN)) {
      stop(paste0("The passed clusterization has the wrong number of elements [",
                  length(clusters), "] instead of the expected number of cells [",
                  getNumCells(objCOTAN), "]."))
    }
    if (!isa(coexDF, "data.frame")) {
      stop(paste0("'clusterCoex' is supposedly composed of data.frames.",
                  " A '", class(coexDF), "' was given instead for clusterization '",
                  clusterizationName, "'."))
    }

    objCOTAN@metaCells <- setColumnInDF(objCOTAN@metaCells, clusters,
                                        clName, getCells(objCOTAN))

    objCOTAN@clustersCoex <- append(objCOTAN@clustersCoex, coexDF)
    names(objCOTAN@clustersCoex)[length(objCOTAN@clustersCoex)] <- clName

    validObject(objCOTAN)

    return(objCOTAN)
  }
)

#' addClusterizationCoex
#'
#' This function adds a clusterization coex data.frame to the current COTAN object.
#' It requires the clusterization to be already present, see \link{\code{addClusterization}}
#'
#' @param objCOTAN a COTAN object
#' @param clusterizationName the name of an existing clusterization
#' @param coexDF a data.frame where each column indicates the coex
#' for each (or some) of the clusters of the clusterization
#'
#' @return the updated COTAN object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#'
#' data("ERCC.cotan")
#' ERCC.cotan <- addClusterizationCoex(ERCC.cotan, "merge", coex)
#'
#' @rdname addClusterizationCoex
setMethod(
  "addClusterizationCoex",
  "COTAN",
  function(objCOTAN, clusterizationName, coexDF) {
    if (!isa(coexDF, "data.frame")) {
      stop(paste0("'clusterCoex' is supposedly composed of data.frames.",
                  " A '", class(coexDF), "' was given instead for clusterization '",
                  clusterizationName, "'."))
    }

    clName <- clusterizationName
    if (!startsWith(clName, "CL_")) {
      clName <- paste0("CL_", clName)
    }

    if (!clName %in% names(getClustersCoex(objCOTAN))) {
      stop(paste0("A clusterization with name '", clusterizationName, "' does not exists."))
    }
    objCOTAN@clustersCoex[[clName]] <- coexDF

    return(objCOTAN)
  }
)


#' dropClusterization
#'
#' This function drops a clusterization from the current COTAN object, by removing
#' the corresponding column in the \code{metaCells} data.frame and in case dropping
#' the corresponding coex data.frame from the \code{clustersCoex} list.
#'
#' @param objCOTAN a COTAN object
#' @param clusterizationName the name of an existing clusterization.
#'
#' @return the updated COTAN object
#'
#' @importFrom rlang is_empty
#'
#' @export
#'
#' @examples
#' data("ERCC.cotan")
#' ERCC.cotan <- dropClusterization(ERCC.cotan, "merged")
#'
#' @rdname dropClusterization
setMethod(
  "dropClusterization",
  "COTAN",
  function(objCOTAN, clusterizationName) {
    clName <- clusterizationName
    if (!startsWith(clName, "CL_")) {
      clName <- paste0("CL_", clName)
    }

    if (!clName %in% names(getClustersCoex(objCOTAN))) {
      stop(paste0("A clusterization with name '", clusterizationName, "' does not exists."))
    }

    keptCols <- !colnames(objCOTAN@metaCells) %in% clName
    objCOTAN@metaCells <- objCOTAN@metaCells[, keptCols]

    keptCols <- !colnames(objCOTAN@clustersCoex) %in% clName
    objCOTAN@clustersCoex <- objCOTAN@clustersCoex[keptCols]

    return(objCOTAN)
  }
)
