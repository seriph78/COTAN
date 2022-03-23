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
setMethod("get.coex","scCOTAN",
          function(object) {
              print("coex dataframe creation")
              hk <- object@hk
              ll <- obs_ct(object)

              object <- ll$object
              ll$object <- NA

              gc()
              ll$no_no <- ll$no_no[!rownames(ll$no_no) %in% hk,
                                   !colnames(ll$no_no) %in% hk]
              ll$yes_no <- ll$yes_no[!rownames(ll$yes_no) %in% hk,
                                     !colnames(ll$yes_no) %in% hk]
              ll$no_yes <- t(ll$yes_no)


              si_si <- object@yes_yes[!rownames(object@yes_yes) %in% hk,!colnames(object@yes_yes)
                                      %in% hk]

              est <- expected_ct(object)
              gc()


              new_estimator_yes_yes_n <- mat2vec_rfast(est$estimator_yes_yes)
              new_estimator_yes_yes_n$coex.values[new_estimator_yes_yes_n$coex.values < 1] <- 1
              new_estimator_yes_no_n <- mat2vec_rfast(est$estimator_yes_no)
              new_estimator_yes_no_n$coex.values[new_estimator_yes_no_n$coex.values < 1] <- 1
              new_estimator_no_no_n <- mat2vec_rfast(est$estimator_no_no)
              new_estimator_no_no_n$coex.values[new_estimator_no_no_n$coex.values < 1] <- 1
              new_estimator_no_yes_n <- mat2vec_rfast(est$estimator_no_yes)
              new_estimator_no_yes_n$coex.values[new_estimator_no_yes_n$coex.values < 1] <- 1


              print("coex estimation")

              yes_yes_n <- mat2vec_rfast(as.matrix(si_si))
              yes_no_n <- mat2vec_rfast(as.matrix(ll$yes_no))
              no_no_n <- mat2vec_rfast(as.matrix(ll$no_no))
              no_yes_n <- mat2vec_rfast(as.matrix(ll$no_yes))

              coex_n <- (yes_yes_n$coex.values - mat2vec_rfast(as.matrix(est$estimator_yes_yes))$coex.values)/new_estimator_yes_yes_n$coex.values +
                  (no_no_n$coex.values - mat2vec_rfast(as.matrix(est$estimator_no_no))$coex.values)/new_estimator_no_no_n$coex.values -
                  (no_yes_n$coex.values - mat2vec_rfast(as.matrix(est$estimator_no_yes))$coex.values)/new_estimator_no_yes_n$coex.values -
                  (yes_no_n$coex.values - mat2vec_rfast(as.matrix(est$estimator_yes_no))$coex.values)/new_estimator_yes_no_n$coex.values

              print("Cleaning RAM 1")
              rm(ll)
              rm(est)
              gc()
              sum.for.div_n <- (1/new_estimator_yes_yes_n$coex.values + 1/new_estimator_no_no_n$coex.values + 1/
                                  new_estimator_yes_no_n$coex.values + 1/new_estimator_no_yes_n$coex.values)

              print("Cleaning RAM 2")
              rm(new_estimator_yes_yes_n)
              rm(new_estimator_no_no_n)
              rm(new_estimator_yes_no_n)
              rm(new_estimator_no_yes_n)
              gc()

              coex_n <- coex_n / sqrt(sum.for.div_n)

              coex_n <- coex_n / sqrt(object@n_cells)

              coex_n <- list("genes"=yes_no_n$genes,"coex.values" =coex_n)

              object@coex <- coex_n#spMat(coex)
              rm(coex_n)
              gc()

              return(object)
          }
)
