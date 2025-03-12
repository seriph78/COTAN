## 2.6.4

Solved issue with use of `parallelDist::parDist()` by allowing a fall-back to
`stats::dist()`. This was needed to address failures running `parDist()` 
on Linux-aarch64 machines

## 2.6.3

Solved bug causing errors while using `torch` with a `CPU` device

Ensure the drop out cluster from `cellsUnifromClustersing()` [`-1`]
keeps its name if it has not been merged at the end of 
the function `mergeUniformCellsClusters()`

Stopped using broken BioConductor PCAtools::pca:
using BioSingular::runPCA() instead

Added new utility function `asClusterization()` that takes any input
representing a *clusterization* (`factor`, `vector` or `data.frame`)
and makes it into a `factor`. This function is now used in all functions taking
in a *clusterization* to standardize the given input

Added `initialIteration` as input parameter to *clusterization* functions so
to avoid override of partial data when the functions are being restarted

## 2.6.2

Added possibility to specify genes' selection method used in the
`cellsUniformClustering()` function

Improved `"AdvancedGDIUniformityCheck"` by adding a new check about 99% quantile

Solved issue in `clustersMarkersHeatmapPlot()`: now passing in a genes not in
the data-set will result in a corresponding gray column instead of an error

Added possibility to specify `clean()` thresholds in the functions
`proceedToCoex()` and `automaticCOTANObjectCreation()`

## 2.6.1

Improved function `clustersMarkersHeatmapPlot()`: now its shows a column for
each marker gene and the shown content is more expressive

Marked the function `clustersDeltaExpression()` as legacy: it has been replaced
with the function `DEAOnClusters()` in the package

Fixed minor bug in class `AdvancedGDIUniformityCheck` regarding third check:
was testing third highest GDI value instead of second

## 2.5.11

Fixed bug in the function `cellsUMAPPlot()`: restored possibility of passing a
genes `vector` as `genesSel` parameter. Also updated the documentation about the
available genes selection methods

## 2.5.10

Fixed typo in error message

Fixed bug in function `genesCoexSpace()`: now `primaryMarkers` can have only a
single gene

## 2.5.9

Exported utility functions about names arrays: `conditionsFromNames()`,
`niceFactorLevels()`, `factorToVector()`

## 2.5.8

feature/add_condition_arguments_to_plot_functions
Now the following plot functions take in conditions explicitly, instead of just
instructions to determine them from the cells' names. The changes involved:
`cellSizePlot()`,  `ECDPlot()`, `genesSizePlot()`,
`mitochondrialPercentagePlot()`, `scatterPlot()`

Added possibility to convert `COTAN` objects to/from `SingleCellExperiment`
objects. `SCE` objects created by the `Seurat` package are supported

Hardened arguments' checks for function `UMAPPlot()`

Solved issue with the function `establishGenesClusters()`: it was throwing an
error when one of the sub-lists in the `groupMarkers` argument did contain
only one element

## 2.5.7

Introduced new way to check for the **Uniform-Transcript** property of the 
*clusters* based on multiple thresholds calibrated so that the new method is
more effective at describing really **statistically uniform** *clusters*

Functions `cellsUniformClustering()` and `mergeUniformCellsClusters()` have been
re-factored so to support new class hierarchy for **UT** checkers.
This allows user to select which method to use for the checks;
as of now the following methods are supported:
* `"SimpleGDIUniformityCheck"`
* `"AdvancedGDIUniformityCheck"`

Avoided issue with pdf file creation: file handle was not closed
in case of errors

Added possibility of choosing number of features in `seuratHVG()`

Solved minor issue with with clusterization functions in cases when only one
cluster was created

## 2.5.6

Made function `heatmapPlot()` more easy to use and in line with the rest of
the `COTAN` package

Now the method `storeGDI()` can take in the output `data.frame` from
the function `calculateGDI()`

Solved few minor issues with the vignette and changed a few default parameters
in `cellsUMAPPlot()`, `pValueFromDEA()` and `findClustersMarkers()`

## 2.5.5

Stopped function `cellsUniformClustering()` from saving the internally created
`Seurat` object due to possibly long saving times

Split the now deprecated function `getNormalizedData()` into two separated
functions: `getNuNormData()` and `getLogNormData()`

Re-factored function `mergeUniformCellsClusters()` to be more precise:
now it merges clusters starting from the most similar in latest batch and
also runs the merging in multiple steps adjusting gradually the *GDI* threshold
ranging from a very strict up to the user given ones.

Fixed minor bugs in functions `GDIPlot()` and `clustersMarkersHeatmapPlot()`

## 2.5.4

Added possibility to display UMAP plots of cells clusters, using the function
`cellsUMAPPlot()`

## 2.5.3

Updated the vignette to the most recent changes

Allowed user to set the ratio of genes above the threshold allowed
in a Uniform Transcript cluster

## 2.5.2

Solved issue with usage checks about the `torch` library

Allowed user to explicitly **opt-out** from the `torch` library usage:
COTAN will avoid `torch` commands when the option `"COTAN.UseTorch"` is
set to `FALSE`

## 2.5.1

Added support for the `torch` library to help with the heavy lifting
calculations of the *genes' COEX* matrix, with consequent substantial speed-up,
especially when a **GPU** is available on the system

## COTAN 2.5.0

First release in Bioconductor 3.20

## COTAN 2.3.6

Refactored `DEAOnCluster()` to make it run faster.

Now clustering functions dump the `GDI` check results for all clusters

Changed default `GDI` threshold to 1.43

Added new input to `mergeUniformCellsClusters()` to allow proper resume of
interrupted merges

Added possibility to query whether the `COEX` matrix is available in a `COTAN`
object

## COTAN 2.3.5

Made checks more strict when adding a *clusterization* or *condition*

Increased reliability of clustering functions by improved error handling and 
by allowing retry runs on estimators functions

## COTAN 2.3.4

Speed-up of GDI calculation via Rfast package

Added possibility of using distance between clusters based on *Zero-One* matrix
instead of *DEA*

Added average floor to `logFoldChangeOnClusters()` to dampen extreme results
when genes are essentially absent from a cluster.

## COTAN 2.3.3

Added method to handle expression levels' change via log-normalized data:
`logFoldChangeOnClusters()`

Minor fix in the import of operators to align to new version of `roxigen2`

Restored default adjustment method of `pValueFromDEA()` to `"none"`
for backward compatibility reasons

## COTAN 2.3.2

Solved issue with `cleanPlots()` when the number of cells exceeded *65536*

Added methods to calculate the `COEX` matrix only on a subset of the columns

Now the function `pValueFromDEA()` returns the p-value adjusted for multi-test

## COTAN 2.3.1

Stopped using explicit PCA via irlba package:
using BioConductor PCAtools::pca instead

## COTAN 2.3.0

First release in Bioconductor 3.19

## COTAN 2.1.8

Made passing *clusterizations* to `COTAN` functions more easy:
now all functions that take a `COTAN` object and a *clusterization* as input
parameters can also take a *clusterization* name 

Added time-stamps to log entries when written on a log file

Fixed bug in the `clustersMarkersHeatmapPlot` function when given a
*clusterization* not matching the latest added to the `COTAN` object

Fixed issue with the highest possible resolution in `seuratClustering()`
function, needed when large datasets must be split in many clusters

## COTAN 2.1.7

Added new flag to the function `cleanPlots()` to suppress evaluation of the
*PCA* on the normalized data. In particular, this allows to reduce significantly
time spent within the function `checkClusterUniformity()`

Added `initialResolution` parameter to `cellsUniformClustering()`: it allows
users to specify the initial resolution used in the calls to
`Seurat::FindClusters()` method. It now uses the same default as Seurat

Added new method `estimateNuLinearByCluster()` that calculates `nu` ensuring
that its average is 1.0 in each given cluster

## COTAN 2.1.6

Added function `reorderClusterization()`: it reorders the given *clusterization*
so that *near clusters* have also *near labels*

The functions `cellsUniformClustering()` and `mergeUniformCellsClusters()` now
return the result of this new function

Separated p-value calculations from `DEAOnClusters()` into the new function
`pValueFromDEA()`. Those `data.frames` are no longer part of the list returned
by the functions `DEAOnClusters()` and `mergeUniformCellsClusters()`

Added function `getClusters()` to retrieve the wanted clusterization from the
cells' meta-dataset

Added function `calculateGenesCE()`: it returns the cross-entropy between
the expected absence of a gene reading against the observed state

Fixed minor issue with `logThis()` to file: it was always appending a new line
even when `appendLF` was set to `FALSE`

Now `checkClusterUniformity()` returns more GDI stats like the percentage of
genes above threshold or the last percentile of the GDI values

Revamped `mergeUniformCellsClusters()` to select in order all the the most
likely candidates pairs of clusters to merge. Provided new user parameter to
balance the merging of most possible candidates versus the time spent doing so

Improved `dropGenesCells()` method: it now retains all meta-data information
that is not related to the results of the other methods

Added zoomed UDE plot to `cleanPlots()` return. It suggests a possible
cut level for low UDE cells

## COTAN 2.1.5

Improved `mergeUniformCellsClusters()`: now it attempts to merge more
clusters pairs

Now errors in the `seuratClustering()` function are interpreted as remaining
cells not-clustered "-1". This applies mostly to cases when `Seurat` finds only
*singlets*

Added flag `calcCoex` to `proceedToCoex()` and `automaticCOTANObjectCreation()`
functions to allow user not to spend time calculating the *genes' COEX* when not
needed

Solved potential issue in the `clustersMarkersHeatmapPlot()` regarding clusters'
labels

Added new internal function `niceFactorLevels()` that ensures all the factors'
levels will have labels with the same length, via padding  the integers values
with '0' and string values with '_'

Relaxed tolerance on tests comparing against saved data

## COTAN 2.1.4

Speed-up by use of `parallelDist::parDist()` to calculate distances instead of
`stats::dist()`

Fixed regression tests failing on non-Linux architectures

## COTAN 2.1.3

Completed function `clustersMarkersHeatmapPlot()`

Added new utility function `normalizeNameAndLabels()`

Added `mergeClusters()` and `multiMergeClusters()` functions

Added support to `conditions` in cells' meta-data

Now *clusterizations* are stored as `factors`

Fixed `COTAN::validity` method in `AllClasses.R`

## COTAN 2.1.2

Fixed bug in `proceedToCoex()` in cases when `saveObj == TRUE`

## COTAN 2.1.1

Updated `README.md` and `NEWS.md`

Renamed methods dealing with housekeeping genes and fully-expressed cells to use
the more proper names fully-expressed genes and fully-expressing cells

Added possibility to users to set the cutoff and thresholds used by the `clean`
and related methods

## COTAN 2.1.0

First release in Bioconductor 3.18

## COTAN 1.99.4

Solved remaining documentation warnings

## COTAN 1.99.3

Updated the vignette, `README.md` and `NEWS.md`

## COTAN 1.99.2

Dropped second vignette: will be merged in the other one...

## COTAN 1.99.1

Minor bug fixes and new function `clustersMarkersHeatmapPlot()`

## COTAN 1.99.0

Included new functionalities for Bioc 2.17 release:

* created a new `COTAN` class to replace the old `scCOTAN`: this class provides
  internal invariants along a wide host of accessors that allows users to avoid
  peeking inside the class
  
* made a new multi-core implementation of the model parameters estimations and
  `COEX` calculations that achieves much higher speeds.
  
* added new functionality about **gene clusters** starting from given markers
  lists

* added new functionality about **uniform cell clustering** based on the maximum
  *GDI* level in the *cluster*
  
* added function to get a *differential expression* estimation for each cluster
  against background
  
* added function to get an *enrichment score* for each cluster given a list of
  markers specific for the cells' population
  
* added plots to asses data-set information at cleaning stage

## COTAN 0.99.13

After official release. PCA function changed to avoid basilisk and Python.
 
## COTAN 0.99.12

Release before the official release of Bioc 3.15. Main changes:
 - The way in which the COEX matrix is estimated and stored is changed to occupy
   less space and run faster.
 
## COTAN 0.99.10

Initial Bioconductor release
