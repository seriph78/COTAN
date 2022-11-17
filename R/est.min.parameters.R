#' est.min.parameters
#'
#' This function estimates the nu and a parameters to minimize the...
#'
#' @param object A COTAN object
#' @param par A vector of inizial a and log(nu)
#' @importFrom stats nlminb
#' @return vector
#' @export
#'
#' @examples
#' data("ERCC.cotan")
setGeneric("est.min.parameters", function(object, par) standardGeneric("est.min.parameters"))
#' @rdname est.min.parameters
setMethod(
  "est.min.parameters", "scCOTAN",
  function(object, par) {
    function.to.min <- function(par, object) {
      a <- par[seq_along(object@a)]
      nu <- exp(par[(length(object@a) + 1):(length(object@a) + length(object@nu))])

      noHKFlags <- flagNotHousekeepingGenes(object)
      zero.cells <- colSums(as.data.frame(object@raw == 0))
      zero.genes <- rowSums(as.data.frame(object@raw[noHKFlags, ] == 0))

      a_1 <- 1 / a
      head(a_1)

      l.nu <- outer(object@lambda[noHKFlags], nu)

      a_1.matrix <- matrix(a_1, nrow = length(a_1), ncol = ncol(l.nu))

      mat.for.sum <- (((a_1.matrix / (l.nu + a_1.matrix)))^a_1.matrix)

      to.min <- mean((zero.cells - colSums(mat.for.sum))^2) + mean((zero.genes - rowSums(mat.for.sum))^2)
    }

    test.nlminb <- nlminb(
      start = c(a = object@a, nu = log(object@nu)), object = object,
      objective = function.to.min, control = list(
        iter.max = 2,
        abs.tol = 10^(-5),
        trace = 2,
        step.min = 0.001, step.max = 0.1
      )
    )
    return(test.nlminb)
  }
)
