# Single Cell RNA Sequencing Normalization, Dimensionality Reduction, and Integration

In this section we will describe strategies for filtering out "low-quality" cell profiles from scRNAseq data. We will perform quality control our two-sample combined and formatted PBMC scRNAseq data set [`Seurat`] object generated in the previous section.nPortions of this section have been adapted from a previous Tufts HPC [workshop](https://hbctraining.github.io/scRNA-seq/lessons/03_SC_quality_control-setup.html).

## Setting up R environment

We begin by setting up our R environment similar to the previous section.

### R library source

We will be reading in and writing files relative to our `intro_to_scrnaseq`. For simplicity, we will create an R object that is simply a character string that gives this path, and use it as a prefix for reading and writing files.

```R
LIB='/cluster/tufts/hpc/tools/R/4.0.0/'
.libPaths(c("",LIB))
```

### Read in R packages

For this section, requires three R packages:

1. `Seurat` : A package for working with and analysis of scRNAseq data. Comprehensive tutorials for available analyses with the `Seurat` R package are available on the project [website](https://satijalab.org/seurat/).
2. `ggplot2`: The standard package for creating plots in R. Much of the plotting functions are wrappers for `ggplot2` functionality.
3. `cowplot` : A nice package for combining plots into a single figure. Specifically we will make use of of the `plot_grid()` function.

```R
library(Seurat)
library(ggplot2)
library(cowplot)
```

### Set base directory

```R
baseDir <- "~/intro_to_scrnaseq/"
```
