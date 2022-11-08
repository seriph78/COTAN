setGeneric(
  "plotCellsHeatmap",
  function(objCOTAN, cellsName, clusters) standardGeneric("plotCellsHeatmap")
)

setGeneric(
  "estimateNuBisection",
  function(objCOTAN, threshold = 0.001) standardGeneric("estimateNuBisection")
)

setGeneric(
  "coex",
  function(objCOTAN, cells = FALSE) standardGeneric("coex")
)

setGeneric(
  "estimateDispersionBisection",
  function(objCOTAN, cores = 1, step = 200) 
    standardGeneric("estimateDispersionBisection")
)

setGeneric(
  "rawNorm",
  function(objCOTAN) standardGeneric("rawNorm")
)

setGeneric(
  "housekeepingGenes",
  function(objCOTAN) standardGeneric("housekeepingGenes")
)

setGeneric(
  "nCells",
  function(objCOTAN) standardGeneric("nCells")
)

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

#-------------------------------------- get and set
setGeneric(
  "getRawData",
  function(objCOTAN) standardGeneric("getRawData")
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
  "getGenes",
  function(objCOTAN) standardGeneric("getGenes")
)

setGeneric(
  "getMetadataDataset",
  function(objCOTAN) standardGeneric("getMetadataDataset")
)


#-------------------------------------- private
setGeneric(
  "expectedContingencyTables",
  function(objCOTAN, cells = FALSE) standardGeneric("expectedContingencyTables")
)

setGeneric(
  "observedContingencyYY", 
  function(objCOTAN, cells = FALSE) standardGeneric("observedContingencyYY")
)