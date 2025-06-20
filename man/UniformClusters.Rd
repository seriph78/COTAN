% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GDI-plot.R, R/cellsUniformClustering.R,
%   R/checkClusterUniformity.R, R/mergeUniformCellsClusters.R
\name{UniformClusters}
\alias{UniformClusters}
\alias{GDIPlot}
\alias{cellsUniformClustering}
\alias{checkClusterUniformity}
\alias{mergeUniformCellsClusters}
\title{Uniform Clusters}
\usage{
GDIPlot(
  objCOTAN,
  genes,
  condition = "",
  statType = "S",
  GDIThreshold = 1.43,
  GDIIn = NULL
)

cellsUniformClustering(
  objCOTAN,
  checker = NULL,
  GDIThreshold = NaN,
  initialResolution = 0.8,
  maxIterations = 25L,
  cores = 1L,
  optimizeForSpeed = TRUE,
  deviceStr = "cuda",
  useDEA = TRUE,
  distance = NULL,
  useCoexEigen = FALSE,
  dataMethod = "",
  genesSel = "HVG_Seurat",
  numGenes = 2000L,
  numReducedComp = 25L,
  hclustMethod = "ward.D2",
  initialClusters = NULL,
  initialIteration = 1L,
  saveObj = TRUE,
  outDir = "."
)

checkClusterUniformity(
  objCOTAN,
  clusterName,
  cells,
  checker,
  cores = 1L,
  optimizeForSpeed = TRUE,
  deviceStr = "cuda",
  saveObj = TRUE,
  outDir = "."
)

mergeUniformCellsClusters(
  objCOTAN,
  clusters = NULL,
  checkers = NULL,
  GDIThreshold = NaN,
  batchSize = 0L,
  cores = 1L,
  optimizeForSpeed = TRUE,
  deviceStr = "cuda",
  useDEA = TRUE,
  distance = NULL,
  hclustMethod = "ward.D2",
  allCheckResults = data.frame(),
  initialIteration = 1L,
  saveObj = TRUE,
  outDir = "."
)
}
\arguments{
\item{objCOTAN}{a \code{COTAN} object}

\item{genes}{a named \code{list} of genes to label. Each array will have different
color.}

\item{condition}{a string corresponding to the condition/sample (it is used
only for the title).}

\item{statType}{type of statistic to be used. Default is "S": Pearson's
chi-squared test statistics. "G" is G-test statistics}

\item{GDIThreshold}{legacy. The threshold level that is used in a
\linkS4class{SimpleGDIUniformityCheck}. It defaults to \eqn{1.43}}

\item{GDIIn}{when the \code{GDI} data frame was already calculated, it can be put
here to speed up the process (default is \code{NULL})}

\item{checker}{the object that defines the method and the threshold to
discriminate whether a \emph{cluster} is \emph{uniform transcript}. See
\link{UniformTranscriptCheckers} for more details}

\item{initialResolution}{a number indicating how refined are the clusters
before checking for \strong{uniformity}. It defaults to \eqn{0.8}, the same as
\code{\link[Seurat:FindClusters]{Seurat::FindClusters()}}}

\item{maxIterations}{max number of re-clustering iterations. It defaults to
\eqn{25}}

\item{cores}{number of cores to use. Default is 1.}

\item{optimizeForSpeed}{Boolean; when \code{TRUE} \code{COTAN} tries to use the \code{torch}
library to run the matrix calculations. Otherwise, or when the library is
not available will run the slower legacy code}

\item{deviceStr}{On the \code{torch} library enforces which device to use to run
the calculations. Possible values are \code{"cpu"} to us the system \emph{CPU},
\code{"cuda"} to use the system \emph{GPUs} or something like \code{"cuda:0"} to restrict
to a specific device}

\item{useDEA}{Boolean indicating whether to use the \emph{DEA} to define the
distance; alternatively it will use the average \emph{Zero-One} counts, that is
faster but less precise.}

\item{distance}{type of distance to use. Default is \code{"cosine"} for \emph{DEA} and
\code{"euclidean"} for \emph{Zero-One}. Can be chosen among those supported by
\code{\link[parallelDist:parDist]{parallelDist::parDist()}}}

\item{useCoexEigen}{Boolean to determine whether to project the data \code{matrix}
onto the first eigenvectors of the \strong{COEX} \code{matrix} or instead restrict
the data \code{matrix} to the selected genes before applying the \code{PCA} reduction}

\item{dataMethod}{selects the method to use to create the \code{data.frame} to
pass to the \code{\link[=UMAPPlot]{UMAPPlot()}}. See \code{\link[=getDataMatrix]{getDataMatrix()}} for more details.}

\item{genesSel}{Decides whether and how to perform the gene-selection
(defaults to \code{"HVG_Seurat"}). See \code{\link[=getSelectedGenes]{getSelectedGenes()}} for more details.}

\item{numGenes}{the number of genes to select using the above method. Will be
ignored when an explicit list of genes has been passed in}

\item{numReducedComp}{the number of calculated \strong{RDM} components}

\item{hclustMethod}{It defaults is \code{"ward.D2"} but can be any of the methods
defined by the \code{\link[stats:hclust]{stats::hclust()}} function.}

\item{initialClusters}{an existing \emph{clusterization} to use as starting point:
the \emph{clusters} deemed \strong{uniform} will be kept and the remaining cells will
be processed as normal}

\item{initialIteration}{the number associated tot he first iteration; it
defaults to 1. Useful in case of restart of the procedure to avoid
intermediate data override}

\item{saveObj}{Boolean flag; when \code{TRUE} saves intermediate analyses and
plots to file}

\item{outDir}{an existing directory for the analysis output. The effective
output will be paced in a sub-folder.}

\item{clusterName}{the tag of the \emph{cluster}}

\item{cells}{the cells belonging to the \emph{cluster}}

\item{clusters}{The \emph{clusterization} to merge. If not given the last
available \emph{clusterization} will be used, as it is probably the most
significant!}

\item{checkers}{a \code{list} of objects that defines the method and the
\emph{increasing} thresholds to discriminate whether to merge two \emph{clusters} if
deemed \emph{uniform transcript}. See \link{UniformTranscriptCheckers} for more
details}

\item{batchSize}{Number pairs to test in a single round. If none of them
succeeds the merge stops. Defaults to \eqn{2 (\#cl)^{2/3}}}

\item{allCheckResults}{An optional \code{data.frame} with the results of previous
checks about the merging of clusters. Useful to restart the \emph{merging}
process after an interruption.}
}
\value{
\code{GDIPlot()} returns a \code{ggplot2} object with a point got each gene,
where on the ordinates are the \code{GDI} levels and on the abscissa are the
average gene expression (log scaled). Also marked are the given \emph{threshold}
(in red) and the \eqn{50\%} and \eqn{75\%} quantiles (in blue).

\code{cellsUniformClustering()} returns a \code{list} with 2 elements:
\itemize{
\item \code{"clusters"} the newly found cluster labels array
\item \code{"coex"} the associated \code{COEX} \code{data.frame}
}

\code{checkClusterUniformity} returns a checker object of the same type
as the input one, that contains both threshold and results of the check:
see \link{UniformTranscriptCheckers} for more details

a \code{list} with:
\itemize{
\item \code{"clusters"} the merged cluster labels array
\item \code{"coex"} the associated \code{COEX} \code{data.frame}
}
}
\description{
This group of functions takes in input a \code{COTAN} object and
handle the task of dividing the dataset into \strong{Uniform Clusters}, that is
\emph{clusters} that have an homogeneous genes' expression. This condition is
checked by calculating the \code{GDI} of the \emph{cluster} and verifying that no
more than a small fraction of the genes have their \code{GDI} level above the
given \code{GDIThreshold}
}
\details{
\code{GDIPlot()} directly evaluates and plots the \code{GDI} for a sample.

\code{cellsUniformClustering()} finds a \strong{Uniform} \emph{clusterizations} by
means of the \code{GDI}. Once a preliminary \emph{clusterization} is obtained from
the \code{Seurat-package} methods, each \emph{cluster} is checked for \strong{uniformity}
via the function \code{\link[=checkClusterUniformity]{checkClusterUniformity()}}. Once all \emph{clusters} are
checked, all cells from the \strong{non-uniform} clusters are pooled together
for another iteration of the entire process, until all \emph{clusters} are
deemed \strong{uniform}. In the case only a few cells are left out (\eqn{\leq
  50}), those are flagged as \code{"-1"} and the process is stopped.

\code{checkClusterUniformity()} takes a \code{COTAN} object and a cells'
\emph{cluster} and checks whether the latter is \strong{uniform} by looking at the
genes' \code{GDI} distribution. The function runs \code{\link[=checkObjIsUniform]{checkObjIsUniform()}} on the
given input checker

\code{mergeUniformCellsClusters()} takes in a \strong{uniform}
\emph{clusterization} and progressively checks whether \emph{merging} two \emph{near
clusters} would form a \strong{uniform} \emph{cluster} still. Multiple thresholds
will be used from \eqn{1.37} up to the given one in order to prioritize
merge of the best fitting pairs.

This function uses the \emph{cosine distance} to establish the \emph{nearest clusters
pairs}. It will use the \code{\link[=checkClusterUniformity]{checkClusterUniformity()}} function to check
whether the merged \emph{clusters} are \strong{uniform}. The function will stop once
no \emph{tested pairs} of clusters can be \emph{merged} after testing all pairs in a
single batch
}
\examples{
data("test.dataset")

objCOTAN <- automaticCOTANObjectCreation(raw = test.dataset,
                                         GEO = "S",
                                         sequencingMethod = "10X",
                                         sampleCondition = "Test",
                                         cores = 6L,
                                         saveObj = FALSE)

groupMarkers <- list(G1 = c("g-000010", "g-000020", "g-000030"),
                     G2 = c("g-000300", "g-000330"),
                     G3 = c("g-000510", "g-000530", "g-000550",
                            "g-000570", "g-000590"))

gdiPlot <- GDIPlot(objCOTAN, genes = groupMarkers, cond = "test")
plot(gdiPlot)

## Here we override the default checker as a way to reduce the number of
## clusters as higher thresholds imply less stringent uniformity checks
##
## In real applications it might be appropriate to do so in the cases when
## the wanted resolution is lower such as in the early stages of the analysis
##

checker <- new("AdvancedGDIUniformityCheck")
stopifnot(identical(checker@firstCheck@GDIThreshold, 1.297))

checker2 <- shiftCheckerThresholds(checker, 0.1)
stopifnot(identical(checker2@firstCheck@GDIThreshold, 1.397))

splitList <- cellsUniformClustering(objCOTAN, cores = 6L,
                                    optimizeForSpeed = TRUE,
                                    deviceStr = "cuda",
                                    initialResolution = 0.8,
                                    checker = checker2,
                                    saveObj = FALSE)

clusters <- splitList[["clusters"]]

firstCluster <- getCells(objCOTAN)[clusters \%in\% clusters[[1L]]]

checkerRes <-
  checkClusterUniformity(objCOTAN, checker = checker2,
                         clusterName = clusters[[1L]], cells = firstCluster,
                         cores = 6L, optimizeForSpeed = TRUE,
                         deviceStr = "cuda", saveObj = FALSE)

objCOTAN <- addClusterization(objCOTAN,
                              clName = "split",
                              clusters = clusters,
                              coexDF = splitList[["coex"]],
                              override = FALSE)

stopifnot(identical(reorderClusterization(objCOTAN)[["clusters"]], clusters))

## It is possible to pass a list of checkers tot the merge function that will
## be applied each to the *resulting* merged *clusterization* obtained using
## the previous checker. This ensures that the most similar clusters are
## merged first improving the overall performance

mergedList <- mergeUniformCellsClusters(objCOTAN,
                                        checkers = c(checker, checker2),
                                        batchSize = 2L,
                                        clusters = clusters,
                                        cores = 6L,
                                        optimizeForSpeed = TRUE,
                                        deviceStr = "cpu",
                                        distance = "cosine",
                                        hclustMethod = "ward.D2",
                                        saveObj = FALSE)

objCOTAN <- addClusterization(objCOTAN,
                              clName = "merged",
                              clusters = mergedList[["clusters"]],
                              coexDF = mergedList[["coex"]],
                              override = TRUE)

stopifnot(identical(reorderClusterization(objCOTAN)[["clusters"]],
                    mergedList[["clusters"]]))

}
