#' estimates by bisection the 'nu' field of a COTAN object
#' @param objCOTAN a COTAN object
#' @param threshold real value close to zero
#' @return A COTAN object
#' @export
setMethod(
  "estimateNuBisection",
  "COTAN",
  function(objCOTAN, threshold, maxIterations) {
    zeroOne <- sign(objCOTAN@raw)

    # parameters estimation
    if (is_empty(objCOTAN@lambda)) {
      objCOTAN <- estimateLambdaLinear(objCOTAN)
    }

    if(is_empty(objCOTAN@hkGenes)){
      objCOTAN <- housekeepingGenes(objCOTAN)
    }

    # only genes not in housekeeping are used
    lambda <- objCOTAN@lambda[flagNotHousekeepingGenes(objCOTAN)]

    if (is_empty(objCOTAN@nu)) {
      warning("nu vector is empty, estimated linearly initially")
      objCOTAN <- estimateNuLinear(objCOTAN)
    }

    # vector of zeros
    diffZero <- numeric(getNumCells(objCOTAN))

    for (i in seq(1, getNumCells(objCOTAN))) {
      # estimation
      numZeroEst <- sum(funProbZero(
        objCOTAN@dispersion,
        objCOTAN@nu[i] * lambda
      ))
      numZeroObs <- sum(zeroOne[, i] == 0)

      leftNu <- objCOTAN@nu[i]

      #difference of zeros
      leftDiffZero <- numZeroEst - numZeroObs

      # if nu[i] produces a difference of zeros equals to 0, skip
      if (abs(leftDiffZero) < threshold) {
        next
      }

      # Two values of nu are required to perform bisection, one that generates a
      # positive value and the other positive in the difference of zeros
      rightNu <- leftNu
      rightDiffZero <- leftDiffZero

      if (numZeroEst - numZeroObs < 0) {
        # we move the two values until leftNu produces a positive leftDiffZero
        leftNu <- leftNu / 2 # move leftNu

        # estimation
        leftDiffZero <- sum(funProbZero(
          objCOTAN@dispersion,
          leftNu * lambda
        ))
        leftDiffZero <- leftDiffZero - numZeroObs

        while (leftDiffZero < 0) {
          # leftNu still produces a negative difference in the number of zeros,
          # but it is more accurate than rightNu
          rightNu <- leftNu # move rightNu
          rightDiffZero <- leftDiffZero

          leftNu <- leftNu / 2 # move leftNu

          # estimation
          leftDiffZero <- sum(funProbZero(
            objCOTAN@dispersion,
            leftNu * lambda
          ))
          leftDiffZero <- leftDiffZero - numZeroObs
        }
      } else {
        # we move the two values until rightNu produces a positive rightDiffZero
        leftDiffZero <- numZeroEst - numZeroObs

        rightNu <- rightNu * 2 # move rightNu

        # estimation
        rightDiffZero <- sum(funProbZero(
          objCOTAN@dispersion,
          rightNu * lambda
        ))
        rightDiffZero <- leftDiffZero - numZeroObs

        while (rightDiffZero > 0) {
          # rightNU still produces a positive difference in the number of zeros,
          # but it is more accurate than leftNu
          leftNu <- rightNu # move leftNu
          leftDiffZero <- rightDiffZero

          rightNu <- rightNu * 2 # move rigthNu

          # estimation
          rightDiffZero <- sum(funProbZero(
            objCOTAN@dispersion,
            rightNu * lambda
          ))
          rightDiffZero <- rightDiffZero - numZeroObs
        }
      }

      # bisection starts here
      middleNu <- (leftNu + rightNu) / 2

      # estimation
      middleDiffZero <- sum(funProbZero(
        objCOTAN@dispersion,
        middleNu * lambda
      ))
      middleDiffZero <- middleDiffZero - numZeroObs

      iterCount <- 0
      while (abs(middleDiffZero) > threshold & iterCount < maxIterations) {
        if (middleDiffZero > 0) {
          leftNu <- middleNu
          leftDiffZero <- middleDiffZero
        } else {
          rightNu <- middleNu
          rightDiffZero <- middleDiffZero
        }

        middleNu <- (leftNu + rightNu) / 2

        # estimation
        middleDiffZero <- sum(funProbZero(
          objCOTAN@dispersion,
          middleNu * lambda
        ))
        middleDiffZero <- middleDiffZero - numZeroObs

        iterCount <- iterCount + 1
      }

      objCOTAN@nu[i] <- middleNu
    }

    return(objCOTAN)
  }
)
