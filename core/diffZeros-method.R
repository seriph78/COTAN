setMethod(
  "diffZeros",
  "numeric",
  function(disp, sumZeros, muEstimatorCell) {
    if (disp > 0) {
      sumZerosEstimator <- sum((1 + disp * muEstimatorCell)^(-1 / disp))
    } else {
      sumZerosEstimator <- sum((exp(-(1 - disp) * muEstimatorCell)))
    }
    
    sumZerosEstimator - sumZeros
  }
)
