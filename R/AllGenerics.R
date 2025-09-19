
setGeneric(
  "proceedToCoex",
  function(objCOTAN, calcCoex = TRUE,
           optimizeForSpeed = TRUE, deviceStr = "cuda", cores = 1L,
           cellsCutoff = 0.003, genesCutoff = 0.002,
           cellsThreshold = 0.99, genesThreshold = 0.99,
           saveObj = FALSE, outDir = ".") {
    standardGeneric("proceedToCoex")
  }
)

# -------- calculate COEX --------

setGeneric(
  "calculateCoex",
  function(objCOTAN, actOnCells = FALSE, returnPPFract = FALSE,
           optimizeForSpeed = TRUE, deviceStr = "cuda") {
    standardGeneric("calculateCoex")
  }
)

# -------- getters --------

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
  "getNumExpressedGenes",
  function(objCOTAN) standardGeneric("getNumExpressedGenes")
)

setGeneric(
  "getGenesSize",
  function(objCOTAN) standardGeneric("getGenesSize")
)

setGeneric(
  "getNumOfExpressingCells",
  function(objCOTAN) standardGeneric("getNumOfExpressingCells")
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
  "getGDI",
  function(objCOTAN) standardGeneric("getGDI")
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
  "flagNotFullyExpressedGenes",
  function(objCOTAN) standardGeneric("flagNotFullyExpressedGenes")
)

setGeneric(
  "flagNotFullyExpressingCells",
  function(objCOTAN) standardGeneric("flagNotFullyExpressingCells")
)

setGeneric(
  "getFullyExpressedGenes",
  function(objCOTAN) standardGeneric("getFullyExpressedGenes")
)

setGeneric(
  "getFullyExpressingCells",
  function(objCOTAN) standardGeneric("getFullyExpressingCells")
)

setGeneric(
  "getGenesCoex",
  function(objCOTAN, genes = vector(mode = "character"),
           zeroDiagonal = TRUE, ignoreSync = FALSE) {
    standardGeneric("getGenesCoex")
  }
)

setGeneric(
  "getCellsCoex",
  function(objCOTAN, cells = vector(mode = "character"),
           zeroDiagonal = TRUE, ignoreSync = FALSE) {
    standardGeneric("getCellsCoex")
  }
)

setGeneric(
  "isCoexAvailable",
  function(objCOTAN, actOnCells = FALSE, ignoreSync = FALSE) {
    standardGeneric("isCoexAvailable")
  }
)

## ------- `COTAN` clusterization and conditions ------

setGeneric(
  "getClusterizations",
  function(objCOTAN, dropNoCoex = FALSE, keepPrefix = FALSE) {
    standardGeneric("getClusterizations")
  }
)

setGeneric(
  "getClusterizationName",
  function(objCOTAN, clName = "", keepPrefix = FALSE) {
    standardGeneric("getClusterizationName")
  }
)

setGeneric(
  "getClusterizationData",
  function(objCOTAN, clName = "") standardGeneric("getClusterizationData")
)

setGeneric(
  "getAllConditions",
  function(objCOTAN, keepPrefix = FALSE) standardGeneric("getAllConditions")
)

setGeneric(
  "getConditionName",
  function(objCOTAN, condName = "", keepPrefix = FALSE) {
    standardGeneric("getConditionName")
  }
)

setGeneric(
  "getCondition",
  function(objCOTAN, condName = "") standardGeneric("getCondition")
)

setGeneric(
  "getDims",
  function(objCOTAN) standardGeneric("getDims")
)

# -------- modifiers --------

setGeneric(
  "initializeMetaDataset",
  function(objCOTAN, GEO, sequencingMethod, sampleCondition) {
    standardGeneric("initializeMetaDataset")
  }
)

setGeneric(
  "addElementToMetaDataset",
  function(objCOTAN, tag, value) standardGeneric("addElementToMetaDataset")
)

setGeneric(
  "findFullyExpressedGenes",
  function(objCOTAN, cellsThreshold = 0.99) {
    standardGeneric("findFullyExpressedGenes")
  }
)

setGeneric(
  "findFullyExpressingCells",
  function(objCOTAN, genesThreshold = 0.99) {
    standardGeneric("findFullyExpressingCells")
  }
)

setGeneric(
  "setLambda",
  function(objCOTAN, lambda) standardGeneric("setLambda")
)

setGeneric(
  "setDispersion",
  function(objCOTAN, dispersion) standardGeneric("setDispersion")
)

setGeneric(
  "setNu",
  function(objCOTAN, nu) standardGeneric("setNu")
)

setGeneric(
  "dropGenesCells",
  function(objCOTAN, genes = vector(mode = "character"),
           cells = vector(mode = "character")) {
    standardGeneric("dropGenesCells")
  }
)

setGeneric(
  "clean",
  function(objCOTAN, cellsCutoff = 0.003, genesCutoff = 0.002,
           cellsThreshold = 0.99, genesThreshold = 0.99) {
    standardGeneric("clean")
  }
)

setGeneric(
  "addClusterization",
  function(objCOTAN, clName, clusters,
           coexDF = data.frame(), override = FALSE) {
    standardGeneric("addClusterization")
  }
)

setGeneric(
  "addClusterizationCoex",
  function(objCOTAN, clName, coexDF) standardGeneric("addClusterizationCoex")
)

setGeneric(
  "dropClusterization",
  function(objCOTAN, clName) standardGeneric("dropClusterization")
)

setGeneric(
  "addCondition",
  function(objCOTAN, condName, conditions, override = FALSE) {
    standardGeneric("addCondition")
  }
)

setGeneric(
  "dropCondition",
  function(objCOTAN, condName) standardGeneric("dropCondition")
)

setGeneric(
  "dropGenesCoex",
  function(objCOTAN) standardGeneric("dropGenesCoex")
)

setGeneric(
  "dropCellsCoex",
  function(objCOTAN) standardGeneric("dropCellsCoex")
)

# -------- estimators --------

setGeneric(
  "estimateLambdaLinear",
  function(objCOTAN) standardGeneric("estimateLambdaLinear")
)

setGeneric(
  "estimateNuLinear",
  function(objCOTAN) standardGeneric("estimateNuLinear")
)

setGeneric(
  "estimateNuLinearByCluster",
  function(objCOTAN, clName = "", clusters = NULL) {
    standardGeneric("estimateNuLinearByCluster")
  }
)

setGeneric(
  "estimateDispersionViaSolver",
  function(objCOTAN, threshold = 0.001, cores = 1L,
           maxIterations = 50L, chunkSize = 1024L) {
    standardGeneric("estimateDispersionViaSolver")
  }
)

setGeneric(
  "estimateNuBisection",
  function(objCOTAN, threshold = 0.001, cores = 1L,
           maxIterations = 50L, chunkSize = 1024L) {
    standardGeneric("estimateNuBisection")
  }
)

setGeneric(
  "estimateDispersionNuBisection",
  function(objCOTAN, threshold = 0.001, cores = 1L,
           maxIterations = 50L, chunkSize = 1024L,
           enforceNuAverageToOne = TRUE) {
    standardGeneric("estimateDispersionNuBisection")
  }
)

setGeneric(
  "estimateDispersionNuNlminb",
  function(objCOTAN, threshold = 0.001,
           maxIterations = 50L, chunkSize = 1024L,
           enforceNuAverageToOne = TRUE) {
    standardGeneric("estimateDispersionNuNlminb")
  }
)

setGeneric(
  "storeGDI",
  function(objCOTAN, genesGDI) standardGeneric("storeGDI")
)

# -------- Uniform Transcript Checks --------

setGeneric(
  "checkObjIsUniform",
  function(currentC, previousC = NULL, objCOTAN = NULL) {
    standardGeneric("checkObjIsUniform")
  }
)

setGeneric(
  "getCheckerThreshold",
  function(checker) standardGeneric("getCheckerThreshold")
)

setGeneric(
  "calculateThresholdShiftToUniformity",
  function(checker) standardGeneric("calculateThresholdShiftToUniformity")
)

setGeneric(
  "shiftCheckerThresholds",
  function(checker, shift) standardGeneric("shiftCheckerThresholds")
)
