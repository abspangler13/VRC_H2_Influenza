library(here)
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(sessioninfo)
library(harmony)

A316 <- readRDS(file = here::here("analysis","data_objects","03_vdj","A316_final_vdj_all.rds"))

## remove mito genes here?

#RNA data
DefaultAssay(A316) <- "RNA"
A316 <- NormalizeData(A316, normalization.method = "LogNormalize")
A316 <- FindVariableFeatures(A316, selection.method = "vst", nfeatures = 2000)
VariableFeatures(A316) <- grep("IG[HKL]V|TRBV", VariableFeatures(A316), invert = TRUE, value = TRUE)
A316 <- ScaleData(A316)
A316 <- RunPCA(A316, features = VariableFeatures(A316))
A316 <- RunHarmony(A316, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")

#make umaps of both PCs and harmony corrected PCs
A316 <- RunUMAP(A316, reduction = "pca", reduction.name = "rna.umap", dims = 1:30, assay = "RNA")
A316 <- RunUMAP(A316, reduction = "harmony", reduction.name = "harmony.rna.umap", dims = 1:30, assay = "RNA")

pdf(file = here::here("analysis","plots","035_integration","rna_DimRed.pdf"))
DimPlot(A316, reduction = "rna.umap",group.by = "Subject")
DimPlot(A316, reduction = "harmony.rna.umap",group.by = "Subject")
DimPlot(A316, reduction = "rna.umap",group.by = "run")
DimPlot(A316, reduction = "harmony.rna.umap",group.by = "run")
dev.off()

#MS4A1 (high for MemB), PRDM1 (high for PB), CD3E (T cells), CD14 (Monocytes) Justification for removing junk cells aka non-bcells
pdf(file = here::here("analysis","plots","035_integration","broad_celltype_markers.pdf"))
FeaturePlot(A316, features = "MS4A1",reduction = "harmony.rna.umap")
FeaturePlot(A316, features = "PRDM1",reduction = "harmony.rna.umap")
FeaturePlot(A316, features = "CD3E",reduction = "harmony.rna.umap")
FeaturePlot(A316, features = "CD14",reduction = "harmony.rna.umap")
FeaturePlot(A316, features = "CD38",reduction = "harmony.rna.umap")
dev.off()

DefaultAssay(A316) <- "Prot"
A316 <- ScaleData(A316, features = rownames(A316))
A316 <- RunPCA(A316, assay = "Prot", slot = "data", features = rownames(A316), reduction.name = "apca")
A316 <- RunHarmony(A316, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")
#make umaps of both PCs and harmony corrected PCs
A316 <- RunUMAP(A316, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
A316 <- RunUMAP(A316, reduction = "harmony.prot", dims = 1:18, assay = "Prot", reduction.name = "harmony.prot.umap", reduction.key = "har_protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("analysis","plots","035_integration","prot_DimRed.pdf"))
DimPlot(A316, reduction = "prot.umap",group.by = "Subject")
DimPlot(A316, reduction = "harmony.prot.umap",group.by = "Subject")
DimPlot(A316, reduction = "prot.umap",group.by = "run")
DimPlot(A316, reduction = "harmony.prot.umap",group.by = "run")
dev.off()

# do clusering on harmony corrected PCs only
A316 <- FindMultiModalNeighbors(A316, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:30, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn") #makes graphs
A316 <- FindClusters(A316, graph.name = "harmony.snn", algorithm = 3, resolution = 0.2, verbose = FALSE)
A316 <- FindClusters(A316, graph.name = "harmony.snn", algorithm = 3, resolution = 0.3, verbose = FALSE)
A316 <- RunUMAP(A316, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("analysis","plots","035_integration","harmony_DimPlot_UMAP_all_0.2.pdf"))
DimPlot(A316, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.2")
DimPlot(A316, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.2", split.by = "orig.ident", ncol = 5)
DimPlot(A316, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.2", split.by = "Subject", ncol = 4)
A316 %>% filter(vac_grp %in% c("4A","4B")) %>% DimPlot(reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.2", split.by = "orig.ident", ncol = 4)
A316 %>% filter(vac_grp %in% c("4A","4B")) %>% DimPlot(reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.2", split.by = "Subject", ncol = 4)
dev.off()

pdf(file = here::here("analysis","plots","035_integration","harmony_DimPlot_UMAP_all_0.3.pdf"))
DimPlot(A316, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.3")
DimPlot(A316, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.3", split.by = "orig.ident", ncol = 5)
DimPlot(A316, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.3", split.by = "Subject", ncol = 4)
A316 %>% filter(vac_grp %in% c("4A","4B")) %>% DimPlot(reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.3", split.by = "orig.ident", ncol = 4)
A316 %>% filter(vac_grp %in% c("4A","4B")) %>% DimPlot(reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.3", split.by = "Subject", ncol = 4)
dev.off()

pdf(file = here::here("analysis","plots","035_integration","feature_counts.pdf"))
VlnPlot(A316, features = "nCount_RNA", pt.size = 0.0)
VlnPlot(A316, features = "nCount_Prot", pt.size = 0.0)

VlnPlot(A316, features = "nCount_RNA",group.by = "Subject", pt.size = 0.0)
VlnPlot(A316, features = "nCount_Prot", group.by = "Subject",pt.size = 0.0)

VlnPlot(A316, features = "nCount_RNA",group.by = "orig.ident", pt.size = 0.0)
VlnPlot(A316, features = "nCount_Prot", group.by = "orig.ident",pt.size = 0.0)
dev.off()

pdf(file = here::here("analysis","plots","035_integration","other_celltype_markers.pdf"))
VlnPlot(A316, features = c("P-CD14","P-CD56","P-CD3","P-CD4"),group.by = "harmony.snn_res.0.2", pt.size = 0,ncol=2)
VlnPlot(A316, features = c("CD14","CD56","CD3","CD4"),group.by = "harmony.snn_res.0.2", pt.size = 0,ncol=2,assay = "RNA",log = TRUE)
dev.off()

saveRDS(A316, file = here::here("analysis","data_objects","035_integration","integrated_A316.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
