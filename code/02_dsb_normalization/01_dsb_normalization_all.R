library(Seurat)
library(tidyverse)
library(sessioninfo)
library(dsb)
library(limma)

# Read in negative and positive objects as determined by hashtag demultiplexting. Load in good object which is after QC.
A316.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_n_all.rds"))
A316.good <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_p_all.rds"))

dim(A316.good)
# [1]  33538 156370

######QC protein libraries########
prot.mult = (3*mad(A316.good$log_nCount_Prot))
prot.upper = median(A316.good$log_nCount_Prot) + prot.mult
A316.good <- A316.good %>% subset(subset = log_nCount_Prot < prot.upper)

dim(A316.good)
# [1]  33538 156147


############################################Do DSB normalization ########################################
# https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html 

isotypes = rownames(A316.good@assays$Prot@counts)[grep("iso",rownames(A316.good@assays$Prot@counts))]

normalized_dsb_matrix_sm = DSBNormalizeProtein(cell_protein_matrix = as.matrix(A316.good@assays$Prot@counts),
                                               empty_drop_matrix = as.matrix(A316.n@assays$Prot@counts), use.isotype.control = TRUE, isotype.control.name.vec = isotypes)

# #####Run second linear correction using Limma to reduce batch effect######
# normalized_dsb_matrix_sm = removeBatchEffect(x = normalized_dsb_matrix_sm, batch = A316.good$run, batch2 = A316.good$Subject)

#Save assay data into normalized data slot 
A316.good <- SetAssayData(A316.good, slot = "data", new.data = normalized_dsb_matrix_sm, assay = "Prot")

saveRDS(A316.good, file = here::here("analysis","data_objects","02_dsb_normalization","A316_dsb_all.rds"))

# A316.VDJ <-A316.good
# DefaultAssay(A316.VDJ) <- "Prot"

# #do we need to set variable features here? VariableFeatures(bm) <- rownames(bm[["ADT"]])
# A316.VDJ <- ScaleData(A316.VDJ, features = rownames(A316.VDJ))
# A316.VDJ <- RunPCA(A316.VDJ, assay = "Prot", slot = "data", features = rownames(A316.VDJ), reduction.name = "apca")
# pdf(file = here::here("analysis","plots","05_clustering","Elbow_VDJ_protein_all.pdf"))
# ElbowPlot(A316.VDJ)
# dev.off()

# A316.VDJ <- RunUMAP(A316.VDJ, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
# pdf(file = here::here("analysis","plots","05_clustering","DimPlot_prot_UMAP_all.pdf"))
# DimPlot(A316.VDJ, reduction = "prot.umap", label = TRUE)
# DimPlot(A316.VDJ, reduction = "prot.umap", label = TRUE, group.by = "run")
# DimPlot(A316.VDJ, reduction = "prot.umap", label = TRUE, group.by = "orig.ident")
# DimPlot(A316.VDJ, reduction = "prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
# dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
