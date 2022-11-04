#' get.coex
#'
#' This function estimates and stores the coex matrix in the coex field.
#' It need to be run after \code{\link{cotan_analysis}}
#' @param object A COTAN object
#'
#' @return It returns a COTAN object
#' @export
#' @rdname get.coex
#' @examples
#'
#' data("ERCC.cotan")
#' obj <- get.coex(ERCC.cotan)
#'
setGeneric("get.coex", function(object) standardGeneric("get.coex"))
#' @rdname get.coex
setMethod(
  "get.coex", "scCOTAN",
  function(object) {
    print("coex dataframe creation")
    hk <- object@hk

    yes_yes <- observedContingencyYY(object)
    si_si <- yes_yes[!rownames(yes_yes) %in% hk,
                     !colnames(yes_yes) %in% hk]
    est <- expected_ct(object)
    gc()

    new_estimator_yes_yes_n <- mat2vec_rfast(est$estimator_yes_yes)
    new_estimator_yes_yes_n$values[new_estimator_yes_yes_n$values < 1] <- 1
    new_estimator_yes_no_n <- mat2vec_rfast(est$estimator_yes_no)
    new_estimator_yes_no_n$values[new_estimator_yes_no_n$values < 1] <- 1
    new_estimator_no_no_n <- mat2vec_rfast(est$estimator_no_no)
    new_estimator_no_no_n$values[new_estimator_no_no_n$values < 1] <- 1
    new_estimator_no_yes_n <- mat2vec_rfast(est$estimator_no_yes)
    new_estimator_no_yes_n$values[new_estimator_no_yes_n$values < 1] <- 1
    print("coex estimation")
    yes_yes_n <- mat2vec_rfast(as.matrix(si_si))
    coex_n <- (yes_yes_n$values - mat2vec_rfast(as.matrix(est$estimator_yes_yes))$values)
    print("Cleaning RAM")
    rm(est)
    sum.for.div_n <- (1 / new_estimator_yes_yes_n$values + 1 / new_estimator_no_no_n$values + 1 /
      new_estimator_yes_no_n$values + 1 / new_estimator_no_yes_n$values)

    rm(new_estimator_yes_yes_n)
    gc()
    coex_n <- coex_n * sqrt(sum.for.div_n)
    coex_n <- coex_n / sqrt(object@n_cells)
    coex_n <- list("genes" = yes_yes_n$genes, "values" = coex_n)
    object@coex <- coex_n
    rm(coex_n)
    gc()
    return(object)
  }
)
