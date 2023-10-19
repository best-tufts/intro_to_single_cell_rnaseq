## Differential Expression 

Differential Expression in scRNAseq has multiple meanings. In this section we will focus on two types:
1) Genes that are overexpressed in on cell-type compared to all other cell-types, known as "Marker Genes".
2) Genes that are statistically different between groups of cells with different phenotypes or conditions, known as "Differentially Expressed Genes". 


```R
LIB='/cluster/tufts/hpc/tools/R/4.0.0/'
.libPaths(c("",LIB))
library(tidyverse)
library(Seurat)
library(openxlsx)
library(metap)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("~/intro_to_scrnaseq/")
```

Load data and select resolution
```R
seurat_integrated = readRDS("results/integrated_seurat.rds")
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
DimPlot(seurat_integrated, label=T)
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

We need to reintegrate, here we load pre-processed:

```R
tcell_int = readRDS("results/integrated_tcell_0_1_7_11.rds")

Idents(object = tcell_int) <- "integrated_snn_res.0.1"
DefaultAssay(tcell_int) = "RNA"
DimPlot(tcell_int, label=T)

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

view(all_conserved_markers)

table(all_conserved_markers$cluster)
write.xlsx(all_conserved_markers, "results/findconservedmarkers_tcell_res_0.1_0_1_7_11.xlsx")
```

Filter the markers
```R
markers_filter = all_conserved_markers  %>%
  dplyr::filter(max_pval<0.05)

table(markers_filter$cluster)

markers_top3 = markers_filter %>%
  group_by(cluster) %>%
  slice_max(order_by=ctrl_avg_log2FC, n=3)

DefaultAssay(tcell) = "RNA"
DotPlot(tcell, features=unique(markers_top3$gene)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

markers_top100 = markers_filter %>%
  group_by(cluster) %>%
  slice_max(order_by=ctrl_avg_log2FC, n=100)
```

GO Functional Enrichment
```R
ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = markers_top100,
                     OrgDb = org.Hs.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP",
                     universe=rownames(tcell_int@assays$RNA@counts))

saveRDS(ck, "results/go_tcell_res_0.1_0_1_7_11_top_100.rds")
write.xlsx(ck, "results/go_tcell_res_0.1_0_1_7_11_top_100.xlsx")

dotplot(ck, show=3)
```
