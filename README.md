# COTAN
Estimating co-expression of cell identity factors in single-cell transcriptomes is crucial to decode new mechanisms of cell state transition. Due to the intrinsic low efficiency of single-cell mRNA profiling, novel computational approaches are required to accurately infer gene co-expression in a cell population. We introduce COTAN, a statistical and computational method to analyze the co-expression of gene pairs at single cell level, providing the foundation for single-cell gene interactome analysis.

Some examples on real datasets can be found at https://seriph78.github.io/Cotan_paper/index.html

This current version can be installed as R package using devtools. Currently the installation was tested on Linux, Windows and Mac but there is one multicore function (mclapply) that is not supported under Windows so there can be some problems.

## Installation

devtools::install_github("seriph78/COTAN")

