# Single Cell RNA Sequencing Clustering

In this section we will describe procedures for clustering scRNAseq.

We will perform these procedures on our two-sample combined, formatted, QC'd, and integrated PBMC scRNAseq data set [`Seurat`] object generated in the previous sections. Portions of this section have been adapted from previous Tufts HPC workshops describing [clustering](https://hbctraining.github.io/scRNA-seq_online/lessons/07_SC_clustering_cells_SCT.html) and [clustering quality control](https://hbctraining.github.io/scRNA-seq_online/lessons/07_SC_clustering_cells_SCT.html).

## Setting up R environment

We begin by setting up our R environment similar to the previous sections.

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

## Read in integrated `Seurat` object

```R
integ_seurat <- readRDS(file.path(baseDir, "data/integrated_seurat.rds"))
```

## Data Clustering

The data clustering workflow from the `Seurat` package is carried out in three main steps

1. Principal component analysis, performed by `RunPCA()`.
2. K-nearest neighbor analysis, performed by `FindNeighbors()`.
3. Cluster estimation, performed by `FindClusters()`

### Principal Component Analysis

We've already performed PCA a number of times. We'll run it just one more time to illustrate the steps in the workflow. Moreover, we'll demonstrate a strategy for choosing important PCs to be used in cluster analysis.

```R
integ_seurat <- RunPCA(integ_seurat)
```

#### Choosing the number of principal components for clustering

Clustering analysis can feel like an analytically "squishy" process. There are many "knobs" we can turn in terms of parameters we can set that can affect our clustering output, and consequent downstream analysis. One such knob is the choice of PCs to use as input in our clustering procedure. If we choose to few PCs we risk failing to include biologically relevant information in our data. If we choose too many PCs we risk including spurious variability in our data. Optimizing this process is a difficult task and we generally rely on heuristics for choosing our parameters.

One such heuristic is plot the relationship between PC ranking and their respective standard deviation, where the standard deviation gives us a sense of the amount of variability each PC captures in our data. The general strategy is to use this "elbow plot", to make a decision about at what PC our standard deviation level out at subsequent PCs.

We can create this plot using the `ElbowPlot()` `Seurat` function.

```R
ElbowPlot(integ_seurat, ndims = 50)
```

![](images/elbowPlot.png)

A "squishy" part here, is where we decide the "elbow" is. Clearly, there is a leveling out after the first 7 PCs. However, we also see a drop after the 17th PC. Generally, it is advisable to heir on the side of more PCs, as choosing too few PCs we run the risk of lacking the resolution to identify rare cell types. Accordingly, we'll move on with 17 as our choice of PC dimension.

Next, we'll use a run UMAP to make an initial assessment as to whether our samples appear well integrated at 17 PCs.

```R
integ_seurat <- RunUMAP(integ_seurat, dims = seq(17))
```

```R
uPlot_dim17 <- UMAPPlot(integ_seurat, group.by = "sample") +
  ggtitle("UMAP",
          subtitle = "# PCs = 17") +
  theme(
    legend.position = "bottom"
  )

uPlot_dim17
```

![](images/dimRed5.png)

It appears that at 17 PCs our data relatively homogeneous in terms of the mixture of cell profiles in lower dimension, so we'll move on to clustering.

### Determine the K-nearest neighbor graph

The `Seurat` clustering workflow is a "graph" based method, which means that it takes as input a graph in which nodes are individual cell profiles and edges are connections between cell profiles, based on some similarity measure. Thus, before clustering we will construct a K-nearest neighbor graph (KNN). "Near" refers to the similarity of cell profiles
in our reduced 17 dimension data set. "K" refers to the number of the "nearest" neighbors to choose for each cell profile.

In `Seurat`, KNN is performed using the `FindNeighbors()` function. If you look at the documentation you will see that there are many different parameters you can set for this process. For example, K is set to 20, with the argument `k.param`. In practice, researchers generally leave the defaults alone, so we won't get into the weeds about optimizing this step.

```{r}
integ_seurat <- FindNeighbors(object = integ_seurat, 
                                   dims = seq(17))
```

### Perform clustering with the Louvain algorithm

By default, `Seurat` performs clustering on the KNN graph, using the Louvain algorithm. The Louvain algorithm works by recursively merging the most similar cell profiles in clusters, based on an optimization metric, called "modularity". The algorithm will stop after a certain modularity value has been reached, yeilding the final cluster estimates. 

In `Seurat` the Louvain algorithm is performed by the `FindClusters()` function. Again, there are knobs we can turn in this process, particularly the `resolution` parameter, which controls the stopping value. This is set to 0.8, by default, and lower values will result in earlier stoppage, i.e. fewer clusters, and the opposition for higher values. here, we can take a heuristical approach by clustering under several resolution parameters, and then view our results in lower dimension with UMAP. Again, this is a "squishy" process, and prior knowledge about the general types of cells we expect in our data is very helpful.

Here we will run `FindClusters()` with three levels of resolution by setting `resolution = c(0.4, 0.8, 1.4)`. This will create three new columns in our metadata, "integrated_snn_res.0.4", "integrated_snn_res.0.8", "integrated_snn_res.1.4", each will the cluster assignments under each parameter.

```{r}
integ_seurat <- FindClusters(object = integ_seurat,
                                  resolution = c(0.4, 0.8, 1.4))
```

Next, we plot the UMAPs, color and labeling the cell profiles each different resolution.

```R
uPlotCl0p4 <- UMAPPlot(integ_seurat, group.by = "integrated_snn_res.0.4", label = TRUE) +
  ggtitle("UMAP", 
          subtitle = "Res. = 0.4, # PCs = 17") +
  theme(
    legend.position = "none"
  )

uPlotCl0p8 <- UMAPPlot(integ_seurat, group.by = "integrated_snn_res.0.8", label = TRUE) +
  ggtitle("UMAP", 
          subtitle = "Res. = 0.8, # PCs = 17") +
  theme(
    legend.position = "none"
  )

uPlotCl1p4 <- UMAPPlot(integ_seurat, group.by = "integrated_snn_res.1.4", label = TRUE) +
  ggtitle("UMAP", 
          subtitle = "Res. = 1.4, # PCs = 17") +
  theme(
    legend.position = "none"
  )

plot_grid(uPlotCl0p4, uPlotCl0p8, uPlotCl1p4, nrow = 1)
```

![](images/UMAPclust_pc17.png)

Here, we can see that at resolution = 0.4. We appear to capture clusters of cell profiles that appear to be most similar in terms of their "bunching". At higher resolutions, the data appears like it might be "over clustered", i.e. some of the cluster assignments appear arbitrary.

## Inspecting cluster results

### Visualize distributions of QC variables in clusters

```R
VlnPlot(integ_seurat, 
        group.by = "integrated_snn_res.0.4", 
        features = c("nUMI", "nGene", "log10GenesPerUMI", "percMitoUMI"), 
        ncol = 1)
```

![](images/clCheck17.png)

### Check cell cluster assignments accross samples

```R
table(integ_seurat[[c("integrated_snn_res.0.4", "sample")]])
```
![](images/clVsamp.png)

## Save seurat object

```R
saveRDS(integ_seurat, file.path(baseDir, "data/clustered_seurat.rds"))
```
