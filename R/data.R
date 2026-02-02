# ---------- Dataset definitions ----------

#' @title Data-sets
#'
#' @description Simple data-sets included in the package
#'
#' @name Datasets
#'
NULL

# Raw sample dataset

#' @details `raw.dataset` is a sub-sample of a real *scRNA-seq* data-set
#'
#' @format `raw.dataset` is a data frame with \eqn{2000} genes and \eqn{815}
#'   cells
#'
#' @docType data
#'
#' @usage data(raw.dataset)
#'
#' @source `GEO` `GSM2861514`
#'
#' @rdname Datasets
#'
"raw.dataset"

# Raw `ERCC` dataset

#' @details `ERCCRaw` dataset
#'
#' @format `ERCCRaw` is a `data.frame`
#'
#' @docType data
#'
#' @usage data(ERCCraw)
#'
#' @source `ERCC`
#'
#' @rdname Datasets
#'
"ERCCraw"

#' @details `test.dataset` is an artificial data set obtained by sampling target
#'   negative binomial distributions on a set of \eqn{600} genes on \eqn{2} two
#'   cells *clusters* of \eqn{600} cells each. Each *clusters* has its own set
#'   of parameters for the distributions even, but a fraction of the genes has
#'   the same expression in both *clusters*.
#'
#' @format `test.dataset` is a `data.frame` with \eqn{600} genes and \eqn{1200}
#'   cells
#'
#' @docType data
#'
#' @usage data(test.dataset)
#'
#' @rdname Datasets
#'
"test.dataset"

#' @details `test.dataset.clusters1` is the *clusterization* obtained running
#'   `cellsUniformClustering()` on the `test.dataset`
#'
#' @format `test.dataset.clusters1` is a `character array`
#'
#' @docType data
#'
#' @usage data(test.dataset.clusters1)
#'
#' @rdname Datasets
#'
"test.dataset.clusters1"

#' @details `test.dataset.clusters2` is the *clusterization* obtained running
#'   `mergeUniformCellsClusters()` on the `test.dataset` using the previous
#'   *clusterization*
#'
#' @format `test.dataset.clusters2` is a `character array`
#'
#' @docType data
#'
#' @usage data(test.dataset.clusters2)
#'
#' @rdname Datasets
#'
"test.dataset.clusters2"


#' @details `vignette.split.clusters` is the clusterization obtained running
#'   `cellsUniformClustering()` on the vignette dataset (mouse cortex E17.5,
#'   `GEO`: `GSM2861514`)
#'
#' @format `vignette.split.clusters` is a `factor`
#'
#' @docType data
#'
#' @usage data(vignette.split.clusters)
#'
#' @rdname Datasets
#'
"vignette.split.clusters"


#' @details `vignette.merge.clusters` is the clusterization obtained running
#'   `mergeUniformCellsClusters()` on the vignette dataset (mouse cortex E17.5,
#'   `GEO`: `GSM2861514`) using the previous *clusterization*
#'
#' @format `vignette.merge.clusters` is a `factor`
#'
#' @docType data
#'
#' @usage data(vignette.merge.clusters)
#'
#' @rdname Datasets
#'
"vignette.merge.clusters"


#' @details `vignette.merge2.clusters` is the clusterization obtained re-running
#'   `mergeUniformCellsClusters()` on the vignette dataset (mouse cortex E17.5,
#'   `GEO`: `GSM2861514`) using the `vignette.split.clusters` *clusterization*,
#'   but with a sequence of progressively relaxed checks
#'
#' @format `vignette.merge2.clusters` is a `factor`
#'
#' @docType data
#'
#' @usage data(vignette.merge2.clusters)
#'
#' @rdname Datasets
#'
"vignette.merge2.clusters"

# ---------- Torch library section ----------

#' @title Installing `torch` `R` library (on `WSL-Linux`)
#'
#' @description A brief explanation of how to install the torch package on
#'   `WSL2` (`Windows Subsystem for Linux`), but it might work the same for
#'   other `Linux` systems. Naturally it makes a difference whether one wants to
#'   install support only for the `CPU` or also have the system `GPU` at the
#'   ready!
#'
#'   These are just *suggestions*, no guarantees are given:
#'   **try at your own peril!**
#'
#' @description The main resources to install `torch` is
#'   \url{https://torch.mlverse.org/docs/articles/installation.html} or
#'   \url{https://cran.r-project.org/web/packages/torch/vignettes/installation.html}
#'
#' @details For the `CPU`-only support one need to ensure that also numeric
#'   libraries are installed, like `BLAS` and `LAPACK` and/or `MKL` if your
#'   `CPU` is from *Intel*. Otherwise `torch` will be stuck at using a single
#'   core for all computations.
#'
#' @details For the `GPU`, currently only `cuda` devices are supported. Moreover
#'   only some specific versions of `cuda` (and corresponding `cudnn`) are
#'   effectively usable, so one needs to install them to actually use the `GPU`.
#'
#'   As of today only `cuda` 12.8 is supported, but check the `torch`
#'   documentation for more up-to-date information. Before downgrading your
#'   `cuda` version, please be aware that it is possible to maintain separate
#'   main versions of `cuda` at the same time on the system.
#'
#'   Below a link to install `cuda` 12.8 for `WSL2` given: use a local installer
#'   to be sure the wanted `cuda` version is being installed, and not the latest
#'   one: [`cuda` 12.8 for `WSL2`](https://developer.nvidia.com/cuda-12-8-0-download-archive?target_os=Linux&target_arch=x86_64&Distribution=WSL-Ubuntu&target_version=2.0&target_type=deb_local)
#'
#'   It can happen that after the installation of the new `torch` version,
#'   including the dependencies, `torch` actually fails to run claiming that
#'   `lantern` dependency was not installed.
#'   This can happen due to presence of old cache files, so one can simply
#'   delete the relevant cache folder as in the example below...
#'
#' @examplesIf FALSE
#'   # find the relevant cache folder and delete it
#'   folder <- tools::R_user_dir("torch","cache")
#'   print(folder)
#'   unlink(folder, recursive = TRUE, force = TRUE)
#'
#'   # reinstall explicitly for the toolchain you want
#'   torch::install_torch(cuda_version = "12.8", reinstall = TRUE)
#'
#'   # run this only **AFTER** restarting `R`
#'   COTAN:::canUseTorch(TRUE, "cuda")
#'
#' @name Installing_torch
#'
NULL
