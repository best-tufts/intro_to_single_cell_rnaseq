## Differential Expression 
Differential Expression in scRNAseq has multiple meanings. In this section we will focus on two types:
1) Genes that are overexpressed in on cell-type compared to all other cell-types, known as "Marker Genes".
2) Genes that are statistically different between groups of cells with different phenotypes or conditions, known as "Differentially Expressed Genes". 

To start, we set our library path:
```R
LIB='/cluster/tufts/hpc/tools/R/4.0.0/'
.libPaths(c("",LIB))
```

We require new packages:
1) openxlsx
2) metap
3) clusterProfiler
4) org.Hs.eg.db

```R
library(tidyverse)
library(Seurat)
library(openxlsx)
library(metap)
library(clusterProfiler)
library(org.Hs.eg.db)
```

Set the base dir:
```R
baseDir <- "~/intro_to_scrnaseq/"
```

Load data and select resolution
```R
seurat_integrated = readRDS(file.path(baseDir, "results/clustered_seurat.rds"))
```

Set our identities to be the clusters found at the resolution 0.4 and set the RNA assay to be the default assay:
```R
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
DefaultAssay(seurat_integrated) = "RNA"
```

Visualizations 
```R
FeaturePlot(seurat_integrated, features=c("CD8A"), label=T)

VlnPlot(seurat_integrated, 
        features=c("CD8A"),
        split.by="orig.ident",
        split.plot = TRUE, 
        pt.size=0)

VlnPlot(seurat_integrated, 
        features=c("PYURF"),
        split.by="orig.ident",
        split.plot = TRUE, 
        pt.size=0)
```
look more closely at t-cell subsets

```R
tcell = subset(seurat_integrated, idents = c(0,1,7,11))
DefaultAssay(tcell) = "RNA"
DimPlot(tcell)
```

We need to reintegrate, here we load pre-processed integrated t-cell subset:

```R
tcell_int = readRDS(file.path(baseDir, "data/integrated_tcell.rds"))
```

```R
Idents(object = tcell_int) <- "integrated_snn_res.0.1"
DefaultAssay(tcell_int) = "RNA"
DimPlot(tcell_int, label=T)
```

Let's calculate conserved markers:
```R
all_conserved_markers = data.frame()
clusters = unique(tcell_int$integrated_snn_res.0.1)

for(cl in clusters){
  markers_conserved = FindConservedMarkers(tcell_int,
                                           ident.1 = cl,
                                           grouping.var = "orig.ident",
                                           only.pos = TRUE,
                                           logfc.threshold = 0.25)
  
  markers_conserved$cluster = cl
  markers_conserved$gene= rownames(markers_conserved)
  all_conserved_markers = rbind(all_conserved_markers, markers_conserved)
}
```

View markers
```R 
table(all_conserved_markers$cluster)
view(all_conserved_markers)
```

Write markers:
```R
write.xlsx(all_conserved_markers, 
            file.path(baseDir,"results/findconservedmarkers_tcell_res_0.1.xlsx")
```

Filter the markers
```R
markers_filter = all_conserved_markers  %>%
  dplyr::filter(max_pval<0.05)
```

How many in each cluster:
```R
table(markers_filter$cluster)
```

Let's find the top 3:
```R
markers_top3 = markers_filter %>%
  group_by(cluster) %>%
  slice_max(order_by=ctrl_avg_log2FC, n=3)
```

Plot:
```R
DotPlot(tcell, features=unique(markers_top3$gene)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

Find the top 100:
```R
markers_top100 = markers_filter %>%
  group_by(cluster) %>%
  slice_max(order_by=ctrl_avg_log2FC, n=100)
```

GO Functional Enrichment of top 100:
```R
ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = markers_top100,
                     OrgDb = org.Hs.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP",
                     universe=rownames(tcell_int@assays$RNA@counts))
```

Plot:
```R
dotplot(ck, show=3)
```

Save:
```R
saveRDS(ck, 
        file.path(baseDir,"results/go_tcell_res_0.1_top_100.rds"))
write.xlsx(ck, 
           file.path(baseDir,"results/go_tcell_res_0.1_top_100.xlsx"))
```
