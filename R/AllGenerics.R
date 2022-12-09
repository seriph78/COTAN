
# -------------------------------------- calculate coex

setGeneric(
  "calculateMu",
  function(objCOTAN) standardGeneric("calculateMu")
)

setGeneric(
  "calculateCoex",
  function(objCOTAN, actOnCells = FALSE, optimizeForSpeed = TRUE)
    standardGeneric("calculateCoex")
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
  "getGenesSize",
  function(objCOTAN) standardGeneric("getGenesSize")
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
  "getMetadataDataset",
  function(objCOTAN) standardGeneric("getMetadataDataset")
)

setGeneric(
  "getMetadataElement",
  function(objCOTAN, tag) standardGeneric("getMetadataElement")
)

setGeneric(
  "getMetadataGenes",
  function(objCOTAN) standardGeneric("getMetadataGenes")
)

setGeneric(
  "getMetadataCells",
  function(objCOTAN) standardGeneric("getMetadataCells")
)

setGeneric(
  "getClustersCoex",
  function(objCOTAN) standardGeneric("getClustersCoex")
)

setGeneric(
  "flagNotHousekeepingGenes",
  function(objCOTAN) standardGeneric("flagNotHousekeepingGenes")
)

setGeneric(
  "flagNotFullyExpressedCells",
  function(objCOTAN) standardGeneric("flagNotFullyExpressedCells")
)

setGeneric(
  "getHousekeepingGenes",
  function(objCOTAN) standardGeneric("getHousekeepingGenes")
)

setGeneric(
  "getFullyExpressedCells",
  function(objCOTAN) standardGeneric("getFullyExpressedCells")
)

setGeneric(
  "getGenesCoex",
  function(objCOTAN, genes = c()) standardGeneric("getGenesCoex")
)

setGeneric(
  "getCellsCoex",
  function(objCOTAN, cells = c()) standardGeneric("getCellsCoex")
)

setGeneric(
  "getClusterizations",
  function(objCOTAN, dropNoCoex = FALSE, keepPrefix = FALSE)
    standardGeneric("getClusterizations")
)

setGeneric(
  "getClusterizationData",
  function(objCOTAN, clName = NULL)
    standardGeneric("getClusterizationData")
)

setGeneric(
  "getDims",
  function(objCOTAN) standardGeneric("getDims")
)

#-------------------------------------- modifiers

setGeneric(
  "initializeMetaDataset",
  function(objCOTAN, GEO, sequencingMethod = "10X", sampleCondition) {
    standardGeneric("initializeMetaDataset")
  }
)

setGeneric(
  "addElementToMetaDataset",
  function(objCOTAN, tag, value) standardGeneric("addElementToMetaDataset")
)

setGeneric(
  "findHousekeepingGenes",
  function(objCOTAN) standardGeneric("findHousekeepingGenes")
)

setGeneric(
  "findFullyExpressedCells",
  function(objCOTAN) standardGeneric("findFullyExpressedCells")
)

setGeneric(
  "dropGenesCells",
  function(objCOTAN, genes = c(), cells = c()) standardGeneric("dropGenesCells")
)

setGeneric(
  "clean",
  function(objCOTAN, calcExtraData = TRUE) standardGeneric("clean")
)

setGeneric(
  "addClusterization",
  function(objCOTAN, clName, clusters, coexDF = NULL) standardGeneric("addClusterization")
)

setGeneric(
  "addClusterizationCoex",
  function(objCOTAN, clName, coexDF) standardGeneric("addClusterizationCoex")
)

setGeneric(
  "dropGenesCoex",
  function(objCOTAN, clName) standardGeneric("dropGenesCoex")
)

setGeneric(
  "dropCellsCoex",
  function(objCOTAN, clName) standardGeneric("dropCellsCoex")
)

setGeneric(
  "dropClusterization",
  function(objCOTAN, clName) standardGeneric("dropClusterization")
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
  "estimateDispersionBisection",
  function(objCOTAN, step = 256, threshold = 0.001,
           maxIterations = 1000, cores = 1)
    standardGeneric("estimateDispersionBisection")
)

setGeneric(
  "estimateNuBisection",
  function(objCOTAN, step = 256, threshold = 0.001,
           maxIterations = 1000, cores = 1)
    standardGeneric("estimateNuBisection")
)

setGeneric(
  "estimateDispersionNuBisection",
  function(objCOTAN, step = 256, threshold = 0.001,
           maxIterations = 1000, cores = 1, enforceNuAverageToOne = FALSE)
    standardGeneric("estimateDispersionNuBisection")
)
