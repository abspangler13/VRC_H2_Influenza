library(Seurat)
library(tidyverse)
library(sessioninfo)

##### Do DE and return all genes for fgsea #####
seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering", "H2_Mem_G_clusters.rds"))
seurat <- SetIdent(seurat, value = "Cluster.res0.4")

DefaultAssay(seurat) <- "RNA"
seurat <- ScaleData(seurat, features <- rownames(seurat))

de.genes <- FindAllMarkers(object = seurat, logfc.threshold = 0, return.thresh = 1)

saveRDS(de.genes, file = here::here("analysis","data_objects","06_DE","de_rna_all_genes.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
