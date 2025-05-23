% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/automaticCOTANObjectCreation.R
\name{COTAN_ObjectCreation}
\alias{COTAN_ObjectCreation}
\alias{COTAN}
\alias{proceedToCoex,COTAN-method}
\alias{proceedToCoex}
\alias{automaticCOTANObjectCreation}
\title{\code{COTAN} shortcuts}
\usage{
COTAN(raw = "ANY")

\S4method{proceedToCoex}{COTAN}(
  objCOTAN,
  calcCoex = TRUE,
  optimizeForSpeed = TRUE,
  deviceStr = "cuda",
  cores = 1L,
  cellsCutoff = 0.003,
  genesCutoff = 0.002,
  cellsThreshold = 0.99,
  genesThreshold = 0.99,
  saveObj = TRUE,
  outDir = "."
)

automaticCOTANObjectCreation(
  raw,
  GEO,
  sequencingMethod,
  sampleCondition,
  calcCoex = TRUE,
  optimizeForSpeed = TRUE,
  deviceStr = "cuda",
  cores = 1L,
  cellsCutoff = 0.003,
  genesCutoff = 0.002,
  cellsThreshold = 0.99,
  genesThreshold = 0.99,
  saveObj = TRUE,
  outDir = "."
)
}
\arguments{
\item{raw}{a matrix or dataframe with the raw counts}

\item{objCOTAN}{a newly created \code{COTAN} object}

\item{calcCoex}{a Boolean to determine whether to calculate the genes' \code{COEX}
or stop just after the \code{\link[=estimateDispersionBisection]{estimateDispersionBisection()}} step}

\item{optimizeForSpeed}{Boolean; when \code{TRUE} \code{COTAN} tries to use the \code{torch}
library to run the matrix calculations. Otherwise, or when the library is
not available will run the slower legacy code}

\item{deviceStr}{On the \code{torch} library enforces which device to use to run
the calculations. Possible values are \code{"cpu"} to us the system \emph{CPU},
\code{"cuda"} to use the system \emph{GPUs} or something like \code{"cuda:0"} to restrict
to a specific device}

\item{cores}{number of cores to use. Default is 1.}

\item{cellsCutoff}{\code{clean()} will delete from the \code{raw} data any gene that is
expressed in less cells than threshold times the total number of cells.
Default cutoff is \eqn{0.003 \; (0.3\%)}}

\item{genesCutoff}{\code{clean()} will delete from the \code{raw} data any cell that is
expressing less genes than threshold times the total number of genes.
Default cutoff is \eqn{0.002 \; (0.2\%)}}

\item{cellsThreshold}{any gene that is expressed in more cells than threshold
times the total number of cells will be marked as \strong{fully-expressed}.
Default threshold is \eqn{0.99 \; (99.0\%)}}

\item{genesThreshold}{any cell that is expressing more genes than threshold
times the total number of genes will be marked as \strong{fully-expressing}.
Default threshold is \eqn{0.99 \; (99.0\%)}}

\item{saveObj}{Boolean flag; when \code{TRUE} saves intermediate analyses and
plots to file}

\item{outDir}{an existing directory for the analysis output.}

\item{GEO}{a code reporting the GEO identification or other specific dataset
code}

\item{sequencingMethod}{a string reporting the method used for the sequencing}

\item{sampleCondition}{a string reporting the specific sample condition or
time point.}
}
\value{
a \code{COTAN} object

\code{proceedToCoex()} returns the updated \code{COTAN} object with genes'
\code{COEX} calculated. If asked to, it will also store the object, along all
relevant clean-plots, in the output directory.

\code{automaticCOTANObjectCreation()} returns the new \code{COTAN} object with
genes' \code{COEX} calculated. When asked, it will also store the object, along
all relevant clean-plots, in the output directory.
}
\description{
These functions create a \linkS4class{COTAN} object and/or also run
all the necessary steps until the genes' \code{COEX} matrix is calculated.
}
\details{
Constructor of the class \code{COTAN}

\code{proceedToCoex()} takes a newly created \code{COTAN} object (or the
result of a call to \code{dropGenesCells()}) and runs \code{\link[=calculateCoex]{calculateCoex()}}

\code{automaticCOTANObjectCreation()} takes a raw dataset, creates and
initializes a \code{COTAN} object and runs \code{\link[=proceedToCoex]{proceedToCoex()}}
}
\examples{
data("test.dataset")
obj <- COTAN(raw = test.dataset)

#
# In case one needs to run more steps to clean the datatset
# the following might apply
if (FALSE) {
  objCOTAN <- initializeMetaDataset(objCOTAN,
                                    GEO = "test",
                                    sequencingMethod = "artificial",
                                    sampleCondition = "test dataset")
#
# doing all the cleaning...
#
# in case the genes' `COEX` is not needed it can be skipped
# (e.g. when calling [cellsUniformClustering()])
  objCOTAN <- proceedToCoex(objCOTAN, calcCoex = FALSE,
                            cores = 6L, optimizeForSpeed = TRUE,
                            deviceStr = "cuda", saveObj = FALSE)
}


## Otherwise it is possible to run all at once.
objCOTAN <- automaticCOTANObjectCreation(
  raw = test.dataset,
  GEO = "code",
  sequencingMethod = "10X",
  sampleCondition = "mouse_dataset",
  calcCoex = TRUE,
  saveObj = FALSE,
  outDir = tempdir(),
  cores = 6L)

}
