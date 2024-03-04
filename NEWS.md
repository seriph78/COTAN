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
