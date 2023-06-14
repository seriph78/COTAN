# COTAN v2

A Comprehensive and Versatile Framework for Single-Cell
Gene Co-Expression Studies and Cell Type Identification

## About

The estimation of gene co-expression in single-cell RNA sequencing (scRNA-seq)
is a critical step in the analysis of scRNA-seq data. The low efficiency of
scRNA-seq methodologies makes sensitive computational approaches crucial to
accurately infer transcription profiles in a cell population.

COTAN is a statistical and computational method that analyzes the co-expression
of gene pairs at the single-cell level. It employs an innovative mathematical
model that leads to a generalized contingency table framework. COTAN relies on
the zero unique molecular identifier (UMI) counts distribution instead of
focusing on positive counts to evaluate or extract different scores and
information for gene correlation studies and gene or cell clustering.

COTAN assesses whether gene pairs are correlated or anti-correlated, providing a
new correlation index with an approximate p-value for the associated test of
independence. It also checks whether single genes are differentially expressed,
scoring them with a newly defined global differentiation index (GDI). COTAN
plots and clusters genes according to their co-expression pattern with other
genes to study gene interactions and identify cell-identity markers.

COTAN v2 introduces a novel feature that uses gene GDI values to assess the
biological uniformity of a cell cluster. This feature allows researchers to
apply an iterative cell clustering pipeline and achieve a finer resolution of
uniform clusters. COTAN shows high sensitivity in extracting information from
small clusters and lowly expressed genes. Furthermore, COTAN leverages its
contingency table framework to directly identify genes that are over-represented
or under-represented in the cluster with respect to the rest of the data-set.

COTAN computes an enrichment score for a given list of marker genes, which can
be used to identify and merge small uniform clusters and to check a final
cluster identification.

In version 2.0.0  new functions and plots to check and clean the data-set were
included along several visualization tools to help users explore and interpret
their data. COTAN has a user-friendly interface that is easy to use and does not
require extensive programming skills.
The strength of COTAN is its ability to help researchers better understand
scRNA-seq data. By identifying gene modules, cell types, and new marker genes,
researchers gain insights into the underlying biology of their samples. This
helps disease diagnosis, drug discovery, and other applications. In summary,
COTAN is a powerful and versatile tool for the analysis of scRNA-seq data, with
the potential to facilitate the discovery of new cell types and biological
insights.

## Examples

Main source of examples for the COTAN v2 is the vignette: Guided_tutorial_v2.
There it is illustrated the preparatory cleaning steps, various analysis results
and plots done on the data-set "Mouse Cortex E17.5, GEO:GSM2861514"

Further more it is possible to look at some other examples on real data-sets
at <https://seriph78.github.io/Cotan_paper/index.html> and
<https://seriph78.github.io/COTAN_Datasets_analysis/>.

The first link shows how to handle the genes' clustering while the second shows
how to use the new cells' clustering functions to obtain
**uniform clusterizations**
[Please note: the first link has not been upgraded to the version 2.0.x,
so, while it is possible to reproduce all the steps there described,
they need to be manually adapted to the new interface to be executed]

## Installation

This current version can be installed as R package using devtools. Currently the
installation was tested on Linux, Windows and Mac but please note that due to
lack of multi-core support under Windows it might run slower.

devtools::install_github("seriph78/COTAN")
