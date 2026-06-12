LIB='/cluster/tufts/hpc/tools/R/4.0.0/'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
})

baseDir <- "~/intro_to_scrnaseq/"

integ_seurat = readRDS(file.path(baseDir, "data/clustered_seurat.rds"))

hpca = HumanPrimaryCellAtlasData()

query_counts = integ_seurat@assays$RNA@data

pred_cell <- SingleR(test = query_counts,
                     ref = hpca,
                     assay.type.test="logcounts",
                     labels = hpca$label.main)

saveRDS(pred_cell, file.path(baseDir, "results/singler_hpca_cell.rds"))
