setMethod(
  "estimateNuBisection",
  "COTAN",
  function(objCOTAN, threshold) {
    zeroOne <- sign(objCOTAN@raw)

    # parameters estimation
    if (is_empty(objCOTAN@lambda)) {
      objCOTAN <- estimateLambdaLinear(objCOTAN)
    }

    if (is_empty(objCOTAN@nu)) {
      warning("nu vector is empty, estimated linearly initially")
      objCOTAN <- estimateNuLinear(objCOTAN)
    }

    if (is_empty(objCOTAN@nCells)) {
      objCOTAN <- nCells(objCOTAN)
    }

    # vector of zeros
    diffZero <- numeric(objCOTAN@nCells)

    # looking for the positive and negative value for bisection for each nu
    # element
    for (i in seq(1, objCOTAN@nCells)) {
      numZeroEst <- sum(funProbZero(
        objCOTAN@dispersion,
        objCOTAN@nu[i] * objCOTAN@lambda
      ))
      numZeroObs <- sum(zeroOne[, i] == 0)

      leftNu <- objCOTAN@nu[i]
      leftDiffZero <- numZeroEst - numZeroObs
      
      # if nu[i] produces a difference of zeros equals to 0, skip
      if (abs(leftDiffZero) < threshold) {
        next
      }

      rightNu <- leftNu
      rightDiffZero <- leftDiffZero

      if (numZeroEst - numZeroObs < 0) {
        # we move the two values until leftNu produces a positive leftDiffZero
        leftNu <- leftNu / 2
        leftDiffZero <- sum(funProbZero(
          objCOTAN@dispersion,
          leftNu * objCOTAN@lambda
        ))
        leftDiffZero <- leftDiffZero - numZeroObs

        while (leftDiffZero < 0) {
          # leftNu still produces a negative difference in the number of zeros,
          # but it is more accurate than rightNu
          rightNu <- leftNu
          rightDiffZero <- leftDiffZero

          leftNu <- leftNu / 2
          leftDiffZero <- sum(funProbZero(
            objCOTAN@dispersion,
            leftNu * objCOTAN@lambda
          ))
          leftDiffZero <- leftDiffZero - numZeroObs
        }
      } else {
        # we move the two values until rightNu produces a positive rightDiffZero
        leftDiffZero <- numZeroEst - numZeroObs

        rightNu <- rightNu * 2
        rightDiffZero <- sum(funProbZero(
          objCOTAN@dispersion,
          rightNu * objCOTAN@lambda
        ))
        rightDiffZero <- leftDiffZero - numZeroObs

        while (rightDiffZero > 0) {
          # rightNU still produces a positive difference in the number of zeros,
          # but it is more accurate than leftNu
          leftNu <- rightNu
          leftDiffZero <- rightDiffZero

          rightNu <- rightNu * 2
          rightDiffZero <- sum(funProbZero(
            objCOTAN@dispersion,
            rightNu * objCOTAN@lambda
          ))
          rightDiffZero <- rightDiffZero - numZeroObs
        }
      }

      # bisection
      middleNu <- (leftNu + rightNu) / 2
      middleDiffZero <- sum(funProbZero(
        objCOTAN@dispersion,
        middleNu * objCOTAN@lambda
      ))
      middleDiffZero <- middleDiffZero - numZeroObs

      while (abs(middleDiffZero) > threshold) {
        if (middleDiffZero > 0) {
          leftNu <- middleNu
          leftDiffZero <- middleDiffZero
        } else {
          rightNu <- middleNu
          rightDiffZero <- middleDiffZero
        }

        middleNu <- (leftNu + rightNu) / 2
        middleDiffZero <- sum(funProbZero(
          objCOTAN@dispersion,
          middleNu * objCOTAN@lambda
        ))
        middleDiffZero <- middleDiffZero - numZeroObs
      }

      objCOTAN@nu[i] <- middleNu
    }

    return(objCOTAN)
  }
)
