
# -------------------------------------- calculate coex

setGeneric(
  "calculateMu",
  function(objCOTAN) standardGeneric("calculateMu")
)

setGeneric(
  "observedContingencyTablesYY",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE)
    standardGeneric("observedContingencyTablesYY")
)

setGeneric(
  "observedContingencyTables",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE)
    standardGeneric("observedContingencyTables")
)

setGeneric(
  "expectedContingencyTablesNN",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE, optimizeForSpeed = TRUE)
    standardGeneric("expectedContingencyTablesNN")
)

setGeneric(
  "expectedContingencyTables",
  function(objCOTAN, actOnCells = FALSE, asDspMatrices = FALSE, optimizeForSpeed = TRUE)
    standardGeneric("expectedContingencyTables")
)

setGeneric(
  "calculateCoex",
  function(objCOTAN, actOnCells = FALSE, optimizeForSpeed = TRUE)
    standardGeneric("calculateCoex")
)

setGeneric(
  "calculateS",
  function(objCOTAN, geneSubsetCol = c(), geneSubsetRow = c())
    standardGeneric("calculateS")
)

setGeneric(
  "calculateG",
  function(objCOTAN, geneSubsetCol = c(), geneSubsetRow = c())
    standardGeneric("calculateG")
)

setGeneric(
  "calculateGDI",
  function(objCOTAN, type = "S") standardGeneric("calculateGDI")
)

setGeneric(
  "calculatePValue",
  function(objCOTAN, statType = "S", geneSubsetCol = c(), geneSubsetRow = c() )
    standardGeneric("calculatePValue")
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
  "getMetadataDataset",
  function(objCOTAN) standardGeneric("getMetadataDataset")
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
  "getHousekeepingGenes",
  function(objCOTAN) standardGeneric("getHousekeepingGenes")
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
           maxIterations = 1000, cores = 1)
    standardGeneric("estimateDispersionNuBisection")
)
