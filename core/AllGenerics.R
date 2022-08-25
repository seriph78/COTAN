setGeneric(
  "coex",
  function(objCOTAN) standardGeneric("coex")
)

setGeneric(
  "estimateDispertionBisection",
  function(objCOTAN, cores = 1, step = 200) 
    standardGeneric("estimateDispertionBisection")
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
  "getLambda",
  function(objCOTAN) standardGeneric("getLambda")
)

setGeneric(
  "getNu",
  function(objCOTAN) standardGeneric("getNu")
)

#-------------------------------------- private
setGeneric(
  "expectedContingencyTables",
  function(objCOTAN) standardGeneric("expectedContingencyTables")
)

setGeneric(
  "observedContingencyYY", 
  function(objCOTAN) standardGeneric("observedContingencyYY")
)