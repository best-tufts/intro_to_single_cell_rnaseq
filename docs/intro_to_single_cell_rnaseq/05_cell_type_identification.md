## Cell type identification

As before, we load our libraries and set our working directory:
```R
LIB='/cluster/tufts/hpc/tools/R/4.0.0/'
.libPaths(c("",LIB))
library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
setwd("/cluster/tufts/patralab/rbator01/workshops/2023_06_intro_scrnaseq/")
```

We begin by loading our integrated samples.
```R
seurat_integrated = readRDS("results/integrated_seurat.rds")
```

Set our identities to be the clusters found at the resolution 0.4
```R
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
```

View a UMAP plot of the clusters.
```R
DimPlot(seurat_integrated, label=T)
```
![](images/umap_res0.4.png)


Next, we'll use the [SingleR](https://github.com/LTLA/SingleR) tool with a reference database of expression profiles of known cell types in order to identify our cells and clusters. As mentioned in the lecture, this method measures the correlation of overall gene expression between cells in a reference database with cells in the query dataset in order to label cells  

![](images/singler.png)

To start, we'll use a general database of Human pure cell-types called the Human Primary Cell Type Atlas.  This dataset along with several others is available through the [celldex](https://rdrr.io/github/LTLA/celldex/man/HumanPrimaryCellAtlasData.html) R library. To load:
```R
hpca = HumanPrimaryCellAtlasData()
```

The HPCA object is of the data type called a `Summarized Experiment` which allows one to store count data matrices in assays along with metadata which annotate each cell/sample in the count data.

```R
head(hpca)
```

!!! info "output"
```R
class: SummarizedExperiment 
dim: 6 713 
metadata(0):
assays(1): logcounts
rownames(6): A1BG A1BG-AS1 ... A2M-AS1 A2ML1
rowData names(0):
colnames(713): GSM112490 GSM112491 ... GSM92233 GSM92234
colData names(3): label.main label.fine label.ont
```

Summarized Experiments have the following form:

![](images/summarized_experiment.png)

Well use in particular the label.main column of the metadata, which has the following cell-types:

```R
unique(hpca$label.main)
```

!!! info "output"
```R
 [1] "DC"                   "Smooth_muscle_cells"  "Epithelial_cells"     "B_cell"              
 [5] "Neutrophils"          "T_cells"              "Monocyte"             "Erythroblast"        
 [9] "BM & Prog."           "Endothelial_cells"    "Gametocytes"          "Neurons"             
[13] "Keratinocytes"        "HSC_-G-CSF"           "Macrophage"           "NK_cell"             
[17] "Embryonic_stem_cells" "Tissue_stem_cells"    "Chondrocytes"         "Osteoblasts"         
[21] "BM"                   "Platelets"            "Fibroblasts"          "iPS_cells"           
[25] "Hepatocytes"          "MSC"                  "Neuroepithelial_cell" "Astrocyte"           
[29] "HSC_CD34+"            "CMP"                  "GMP"                  "MEP"                 
[33] "Myelocyte"            "Pre-B_cell_CD34-"     "Pro-B_cell_CD34+"     "Pro-Myelocyte" 
```

Our data to be labeled is input into SingleR as a normalized count matrix, which we can extract from the `RNA` assay our `seurat_integrated` object:
```R
query_counts = seurat_integrated@assays$RNA@data
```

SingleR can be run both on the cluster level and the individual cell level. For cluster-level annotation, the average expression profile of each cluster is used and a single label is generated. This is much faster to run, so we'll start here.

```R
query_clusters = seurat_integrated@meta.data$integrated_snn_res.0.4
```

Next, we'll run SingleR on the cluster level, which should take only a few seconds, and save it to a file.
```R
pred_cluster <- SingleR(test = query_counts,
                        ref = hpca,
                        assay.type.test="logcounts",
                        clusters = query_clusters,
                        labels = hpca$label.main)

saveRDS(pred_cluster, "results/singler_hpca_cluster_res0.4.rds")
```






