#' get.GDI
#'
#' This function produce a dataframe with the GDI for each genes.
#'
#' @param object A COTAN object
#' @param type Type of statistic to be used. Default is "S":
#' Pearson's chi-squared test statistics. "G" is G-test statistics
#'
#' @return A dataframe
#' @export
#' @importFrom stats pchisq
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#' @rdname get.GDI
#' @examples
#' data("ERCC.cotan")
#' quant.p <- get.GDI(ERCC.cotan)
setGeneric("get.GDI", function(object, type = "S") standardGeneric("get.GDI"))
#' @rdname get.GDI
setMethod(
  "get.GDI", "scCOTAN",
  function(object, type = "S") {
    print("function to generate GDI dataframe")
    if (type == "S") {
      print("Using S")
      S <- get.S(object)
    } else if (type == "G") {
      print("Using G")
      S <- get.G(object)
    }

    S <- vec2mat_rfast(S)
    diag(S) <- 0
    CD.sorted <- apply(S, 2, sort, decreasing = TRUE)
    rg <- round(nrow(CD.sorted) / 20, digits = 0)
    CD.sorted <- CD.sorted[seq_len(rg), ]
    CD.sorted <- pchisq(as.matrix(CD.sorted), df = 1, lower.tail = FALSE)

    GDI <- colMeans(CD.sorted)
    GDI <- as.data.frame(GDI)
    colnames(GDI) <- "mean.pval"

    sum.raw.norm <- log(rowSums(getNormalizedData(object)))

    cells <- getRawData(object)
    cells[cells > 0] <- 1
    cells[cells <= 0] <- 0

    exp.cells <- (rowSums(cells) / getNumCells(object)) * 100

    GDI <- merge(GDI, as.data.frame(sum.raw.norm), by = "row.names", all.x = TRUE)
    rownames(GDI) <- GDI$Row.names
    GDI <- GDI[, 2:ncol(GDI)]
    GDI <- merge(GDI, as.data.frame(exp.cells), by = "row.names", all.x = TRUE)
    rownames(GDI) <- GDI$Row.names
    GDI$log.mean.p <- -log(GDI$mean.pval)
    GDI$GDI <- log(GDI$log.mean.p)
    GDI <- GDI[, c("sum.raw.norm", "GDI", "exp.cells")]

    return(GDI)
  }
)
