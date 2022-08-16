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
  "houseKeepingGenes",
  function(objCOTAN) standardGeneric("houseKeepingGenes")
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


## get and set
setGeneric(
  "getLambda",
  function(objCOTAN) standardGeneric("getLambda")
)

setGeneric(
  "getNu",
  function(objCOTAN) standardGeneric("getNu")
)

# private
setGeneric(
  "dispertionBisection",
  function(genes, zeroOneMatrix, muEstimator, threshold = 0.001) 
    standardGeneric("dispertionBisection")
)

setGeneric(
  "diffZeros",
  function(disp, sumZeros, muEstimatorCell) 
    standardGeneric("diffZeros")
)