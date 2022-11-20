
setGeneric(
  "clean",
  function(objCOTAN, calcExtraData = TRUE) standardGeneric("clean")
)

setGeneric(
  "cleanPlots",
  function(objCOTAN, pcaCells, D) standardGeneric("cleanPlots")
)

setGeneric(
  "plotCellsHeatmap",
  function(objCOTAN, cellsName, clusters) standardGeneric("plotCellsHeatmap")
)

setGeneric(
  "estimateNuBisection",
  function(objCOTAN, threshold = 0.001, maxIterations = 1000)
    standardGeneric("estimateNuBisection")
)

# -------------------------------------- calculate coex

setGeneric(
  "observedContingencyYY",
  function(objCOTAN, actOnCells = FALSE) standardGeneric("observedContingencyYY")
)
setGeneric(
  "observedContingency",
  function(objCOTAN, actOnCells = FALSE) standardGeneric("observedContingency")
)

setGeneric(
  "expectedContingencyTables",
  function(objCOTAN, actOnCells = FALSE) standardGeneric("expectedContingencyTables")
)

setGeneric(
  "calculateCoex",
  function(objCOTAN, actOnCells = FALSE) standardGeneric("calculateCoex")
)

setGeneric(
  "calculateS",
  function(objCOTAN) standardGeneric("calculateS")
)

setGeneric(
  "calculateG",
  function(objCOTAN) standardGeneric("calculateG")
)

setGeneric(
  "calculateGDI",
  function(objCOTAN, type = "S") standardGeneric("calculateGDI")
)

#-------------------------------------- getters

setGeneric(
  "getRawData",
  function(objCOTAN) standardGeneric("getRawData")
)

setGeneric(
  "getNumCells",
  function(objCOTAN) standardGeneric("getNumCells")
)

setGeneric(
  "getNumGenes",
  function(objCOTAN) standardGeneric("getNumGenes")
)

setGeneric(
  "getCells",
  function(objCOTAN) standardGeneric("getCells")
)

setGeneric(
  "getGenes",
  function(objCOTAN) standardGeneric("getGenes")
)

setGeneric(
  "getZeroOneProj",
  function(objCOTAN) standardGeneric("getZeroOneProj")
)

setGeneric(
  "getCellsSize",
  function(objCOTAN) standardGeneric("getCellsSize")
)

setGeneric(
  "getNormalizedData",
  function(objCOTAN) standardGeneric("getNormalizedData")
)

setGeneric(
  "getNu",
  function(objCOTAN) standardGeneric("getNu")
)

setGeneric(
  "getLambda",
  function(objCOTAN) standardGeneric("getLambda")
)

setGeneric(
  "getDispersion",
  function(objCOTAN) standardGeneric("getDispersion")
)

setGeneric(
  "getHousekeepingGenes",
  function(objCOTAN) standardGeneric("getHousekeepingGenes")
)

setGeneric(
  "flagNotHousekeepingGenes",
  function(objCOTAN) standardGeneric("flagNotHousekeepingGenes")
)

setGeneric(
  "getMetadataDataset",
  function(objCOTAN) standardGeneric("getMetadataDataset")
)

setGeneric(
  "getCoex",
  function(objCOTAN, asMatrix = TRUE, genes = "all") standardGeneric("getCoex")
)

setGeneric(
  "getCellsCoex",
  function(objCOTAN, asMatrix = TRUE, cells = "all") standardGeneric("getCellsCoex")
)

#-------------------------------------- modifiers

setGeneric(
  "initializeMetaDataset",
  function(objCOTAN, GEO, sequencingMethod = "10X", sampleCondition) {
    standardGeneric("initializeMetaDataset")
  }
)

setGeneric(
  "findHousekeepingGenes",
  function(objCOTAN) standardGeneric("findHousekeepingGenes")
)

setGeneric(
  "dropGenesCells",
  function(objCOTAN, genes = c(), cells = c()) standardGeneric("dropGenesCells")
)

setGeneric(
  "standardizeCoex",
  function(objCOTAN) standardGeneric("standardizeCoex")
)

#-------------------------------------- estimators
setGeneric(
  "estimateLambdaLinear",
  function(objCOTAN) standardGeneric("estimateLambdaLinear")
)

setGeneric(
  "estimateNuLinear",
  function(objCOTAN) standardGeneric("estimateNuLinear")
)

setGeneric(
  "estimateMu",
  function(objCOTAN) standardGeneric("estimateMu")
)

setGeneric(
  "estimateNormalisedData",
  function(objCOTAN) standardGeneric("estimateNormalisedData")
)

setGeneric(
  "runEstimatesLinear",
  function(objCOTAN) standardGeneric("runEstimatesLinear")
)

setGeneric(
  "estimateDispersion",
  function(objCOTAN, cores = 1, step = 200) standardGeneric("estimateDispersion")
)
