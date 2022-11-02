#' get.GDI
#'
#' This function produce a dataframe with the GDI for each genes.
#'
#' @param objCOTAN A COTAN object
#' @param type Type of statistic to be used. Default is "S":
#' Pearson's chi-squared test statistics. "G" is G-test statistics
#'
#' @return A dataframe
#' @export
#' @importFrom stats pchisq
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#' @rdname GDI

setMethod(
  "GDI", 
  "COTAN",
  function(objCOTAN, test = "chiSquare") {
    if (test == "chiSquare") {
      test <- chiSquaredTest(objCOTAN)
    } else if (test == "Gtest") {
      test <- Gtest(objCOTAN)
    } else {
      stop("Type of test unknown. Using 'chiSquare' or 'Gtest'")
    }
    
    # prepare the test result
    test <- vec2mat_rfast(test)
    diag(test) <- 0
    
    testSorted <- apply(test, 2, sort, decreasing = TRUE)
    numberRowUsed <- round(nrow(testSorted) / 20, digits = 0) # 5% of the total
    testSorted <- testSorted[seq_len(numberRowUsed), ]
    testSorted <- pchisq(as.matrix(testSorted), df = 1, lower.tail = FALSE)
    
    GDI <- log(-log(colMeans(testSorted)))
    GDI <- as.data.frame(GDI)
    colnames(GDI) <- "meanPval"
    
    # sum.raw.norm <- log(rowSums(as.matrix(object@raw.norm)))
    # 
    # cells <- as.matrix(object@raw)
    # cells[cells > 0] <- 1
    # cells[cells <= 0] <- 0
    # 
    # exp.cells <- (rowSums(cells) / object@n_cells) * 100
    # 
    # GDI <- merge(GDI, as.data.frame(sum.raw.norm), by = "row.names", all.x = TRUE)
    # rownames(GDI) <- GDI$Row.names
    # GDI <- GDI[, 2:ncol(GDI)]
    # GDI <- merge(GDI, as.data.frame(exp.cells), by = "row.names", all.x = TRUE)
    # rownames(GDI) <- GDI$Row.names
    # GDI$log.mean.p <- -log(GDI$mean.pval)
    # GDI$GDI <- log(GDI$log.mean.p)
    # GDI <- GDI[, c("sum.raw.norm", "GDI", "exp.cells")]
    
    return(GDI)
  }
)
