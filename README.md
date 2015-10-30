# genmat

An R package for organizing, analyzing, and visualizing quantitative genomic data sets in two-dimensional matrices.

---

## Overview

genmat provides a function `matMake` to organize genomic data such as signals (wiggle, bedGraph) and intervals (e.g. bed, narrowPeak) into two-dimensional matrices by aligning such data at defined set of regions in the genome. Each row in a matrix corresponds to a particular region in the genome aligned at some feature, and each column is a base or window within that region. For example, one could align ChIP-seq read densities aligned at transcription start sites. Genomic data organized in this way is convenient for performaing calculations and generating visualizations such as those included in this package..

`matMake` contains many options to define how data is organized into matrices, such as:
- taking into account strandedness of features
- aligning at particular parts of intervals
- setting the genomic size each row represents
- windowing of data, take into account missing data
- metafeature matrices, where genomic intervals of varying sizes are scaled to the same size and a matrix is generated of data within and surrounding these scaled features.
- fragment size matrices, where matrices of bed intervals also store information on the size of the bed interval. These are useful for creating so called "vplots" (Henikoff et. al, 2011 PNAS).

genmat also provides functions for analyzing, visualizing, and performing calculations on matrices generated with `matMake`:
- `matCor` calculates pairwise correlations among a set of matrices
- `matHeamap` generates heatmaps of matrices
- `matHist` plots a histogram of scores in matrices
- `matIqrNorm` normalizes matrices to have a specified interquartile range
- `matLoess` applies loess smoothing to matrices
- `matOps` performs a wide variety of calculations on single matrices (e.g. transformations ) and sets of matrices (e.g. averages, differences, variance, ratios)
- `matPlotAverages` generates aggregate line plots by plotting column means
- `matPlotRows` generates line plots of specific intervals (rows)
- `matQuantileNorm` normalizes matrices with quantile normalization
- `matRead` reads matrices into R as a numeric matrix
- `matTileGrid` generates average line plots of sets of rows depending on their patterns in two other matrices
- `matWindow` smooths matrices using sliding window averages
- `matWrite` save matrix R objects as files
 
## Installation

Install devtools if not installed already:
```R
install.packages("devtools")
```

Then install github-hosted R dependencies and genmat:
```R
# install conifur (convenience functions for R)
devtools::install_github("dvera/conifur")

# install gyro (genomic wrapper scripts in R)
devtools::install_github("dvera/gyro")

# install rubber (R Utilities for Bed and BEdgRaphs
devtools::install_github("dvera/rubber")

# install genmat
devtools::install_github("dvera/genmat")
```

The following R packages are recommended:
```R
devtools::install_github("dvera/converge")
```

