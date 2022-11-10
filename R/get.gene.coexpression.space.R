#' get.gene.coexpression.space
#'
#' To make the GDI more specific, it may be desirable to restrict the set of genes against which
#' GDI is computed to a selected subset V with the recommendation to include a
#' consistent fraction of cell-identity genes, and possibly focusing on markers specific
#' for the biological question of interest (for instance neural cortex layering markers).
#' In this case we denote it as local differentiation index (LDI) relative to V.
#'
#' @param object The COTAN object.
#' @param n.genes.for.marker The number of genes correlated with the primary markers that
#' we want to consider.
#' By default this is 25.
#' @param primary.markers A vector of primary marker names.
#'
#' @return A dataframe
#' @export
#' @importFrom stats quantile
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#' @rdname get.gene.coexpression.space
#' @examples
#' data("ERCC.cotan")
#' df <- get.gene.coexpression.space(ERCC.cotan,
#'   n.genes.for.marker = 10,
#'   primary.markers = getGenes(ERCC.cotan)[sample(getNumGenes(ERCC.cotan), 5)]
#' )
setGeneric("get.gene.coexpression.space", function(object, n.genes.for.marker = 25, primary.markers) {
  standardGeneric("get.gene.coexpression.space")
})
#' @rdname get.gene.coexpression.space
setMethod(
  "get.gene.coexpression.space", "scCOTAN",
  # da sistemare
  function(object, n.genes.for.marker = 25, primary.markers) {
    print("calculating gene coexpression space: output tanh of reduced coex matrix")

    p.val.matrix <- get.pval(object, gene.set.col = primary.markers)
    if (!length(primary.markers) == ncol(p.val.matrix)) {
      print(paste("Gene", primary.markers[!primary.markers %in% colnames(p.val.matrix)],
        "not present!",
        sep = " "
      ))
      primary.markers <- primary.markers[primary.markers %in% colnames(p.val.matrix)]
    }

    all.genes.to.an <- vector()
    for (m in primary.markers) {
      tm <- rownames(p.val.matrix[order(p.val.matrix[, m]), ])[seq_len(n.genes.for.marker)]
      all.genes.to.an <- c(all.genes.to.an, tm)
      all.genes.to.an <- unique(all.genes.to.an)
    }
    all.genes.to.an <- unique(c(all.genes.to.an, primary.markers))
    print(paste0("Secondary markers:", length(all.genes.to.an)))
    tmp <- p.val.matrix[all.genes.to.an, ]
    for (m in primary.markers) {
      tmp <- as.data.frame(tmp[order(tmp[, m]), ])
      tmp$rank <- c(seq_len(nrow(tmp)))
      colnames(tmp)[ncol(tmp)] <- paste("rank", m, sep = ".")
    }
    rank.genes <- tmp[, (length(primary.markers) + 1):ncol(tmp)]
    for (c in seq_along(colnames(rank.genes))) {
      colnames(rank.genes)[c] <- strsplit(colnames(rank.genes)[c],
        split = ".",
        fixed = TRUE
      )[[1]][2]
    }

    S <- get.S(object)
    S <- vec2mat_rfast(S)
    S <- S[, colnames(S) %in% all.genes.to.an]
    S <- as.matrix(S)

    CD.sorted <- t(apply(t(S), 2, sort, decreasing = TRUE))
    rg <- round(ncol(CD.sorted) / 10, digits = 0)
    CD.sorted <- CD.sorted[, seq_len(rg)] # 20
    CD.sorted <- pchisq(as.matrix(CD.sorted), df = 1, lower.tail = FALSE)

    quant.p.val2 <- rowMeans(CD.sorted)
    quant.p.val2 <- as.data.frame(quant.p.val2)
    colnames(quant.p.val2) <- "loc.GDI"

    quant.p.val2$names <- rownames(quant.p.val2)

    quant.p.val2 <- quant.p.val2[quant.p.val2$loc.GDI <=
      stats::quantile(quant.p.val2$loc.GDI, probs = 0.1), ]
    genes.raw <- quant.p.val2$names

    to.plot_cl.genes <- vec2mat_rfast(object@coex)
    to.plot_cl.genes <- to.plot_cl.genes[
      rownames(to.plot_cl.genes) %in% genes.raw,
      colnames(to.plot_cl.genes) %in% all.genes.to.an
    ]

    to.plot_cl.genes <- to.plot_cl.genes * sqrt(getNumCells(object))

    to.plot_cl.genes <- tanh(to.plot_cl.genes)
    print(paste0(
      "Columns (V set) number: ", dim(to.plot_cl.genes)[2],
      " Rows (U set) number: ", dim(to.plot_cl.genes)[1]
    ))

    return(to.plot_cl.genes)
  }
)
