# COTAN v2

This package provides a *comprehensive* and *versatile* framework for single-cell **gene Co-Expression** studies and cell type identification.

## About

The estimation of gene *co-expression* in single-cell RNA sequencing (scRNA-seq) is a critical step in the analysis of scRNA-seq data. The low efficiency of scRNA-seq methodologies makes sensitive computational approaches crucial to accurately infer transcription profiles in a cell population.

`COTAN` is a statistical and computational method that analyzes the co-expression of gene pairs at the single-cell level. It employs an innovative mathematical model that leads to a generalized contingency table framework. `COTAN` relies on the zero *Unique Molecular Identifier* (`UMI`) counts distribution instead of focusing on positive counts to evaluate or extract different scores and information for gene correlation studies and gene or cell clustering.

`COTAN` assesses whether gene pairs are *correlated* or *anti-correlated*, providing a new correlation index with an approximate `p-value` for the associated test of independence. It also checks whether single genes are differentially expressed, scoring them with a newly defined Global Differentiation Index (`GDI`). `COTAN` plots and clusters genes according to their co-expression pattern with other genes to study gene interactions and identify cell-identity markers.

`COTAN v2` introduces a novel feature that uses gene `GDI` values to assess the biological *uniformity* of a cell cluster. This feature allows researchers to apply an iterative cell clustering pipeline and achieve a finer resolution of uniform clusters. `COTAN` shows high sensitivity in extracting information from small clusters and lowly expressed genes. Furthermore, `COTAN` leverages its contingency table framework to directly identify genes that are over-represented or under-represented in the cluster with respect to the rest of the data-set. `COTAN` computes an enrichment score for a given list of marker genes, which can be used to identify and merge small uniform clusters and to check a final cluster identification.

From version 2.0.0 new functions and plots to check and clean the data-set were included along several visualization tools to help users explore and interpret their data. `COTAN` has a user-friendly interface that is easy to use and does not require extensive programming skills. The strength of `COTAN` is its ability to help researchers better understand scRNA-seq data. By identifying gene modules, cell types, and new marker genes, researchers gain insights into the underlying biology of their samples. This helps disease diagnosis, drug discovery, and other applications. In summary, `COTAN` is a powerful and versatile tool for the analysis of scRNA-seq data, with the potential to facilitate the discovery of new cell types and biological insights.

## Examples

Main source of examples for the `COTAN v2` is the *vignette*: [Guided_tutorial_v2](https://github.com/seriph78/COTAN/blob/devel/vignettes/Guided_tutorial_v2.Rmd). There it is illustrated the preparatory cleaning steps, various analysis results and plots done on the data-set `"Mouse Cortex E17.5, <GEO:GSM2861514>"`

Further more it is possible to look at some other examples on real data-sets at [COTAN paper](https://seriph78.github.io/Cotan_paper/index.html) and more extensively at [COTAN Datasets analysis](https://seriph78.github.io/COTAN_Datasets_analysis/).

The first link shows how to handle the genes' clustering while the second shows how to use the new cells' clustering functions to obtain **uniform clusterizations** [Please note: the first link has not been upgraded to the version 2, so, while it is possible to reproduce all the steps there described, they need to be manually adapted to the new interface to be executed]

## Installation

| Build                   | Status                                        |
|-------------------------|-----------------------------------------------|
| BioConductor-release    | [![BioConductor-release](https://images.weserv.nl/?url=bioconductor.org/shields/build/release/bioc/COTAN.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/COTAN) |
| BioConductor-devel      | [![BioConductor-devel](https://images.weserv.nl/?url=bioconductor.org/shields/build/devel/bioc/COTAN.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/COTAN/) |

Check out the user guide on the [Bioconductor landing page - `release`](https://bioconductor.org/packages/release/bioc/html/COTAN.html) (or [`devel`](https://bioconductor.org/packages/devel/bioc/html/COTAN.html)) for more details.

The latest snapshot can be installed directly as `R` package using `devtools`. The installation was tested on *Linux*, *Windows* and *Mac* but please note that due to lack of multi-core support under *Windows* it might run slower. `devtools::install_github("seriph78/COTAN")`

From version 2.5.0 `COTAN` can optionally use the `torch` library and thus, in case, use the `GPU` for some of its calculations, with substantial speed-ups. However this implies a possibly more complicated installation process: see the specific help page [Installing torch](https://github.com/seriph78/COTAN/blob/devel/man/Installing_torch.Rd) for some pointers.
