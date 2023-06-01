## COTAN 2.1.5

Improved `mergeUniformCellsClusters()`: now it attempts to merge more
clusters pairs

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
