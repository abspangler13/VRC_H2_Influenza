#try azimuth https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html

# devtools::install_github("satijalab/seurat-data")
# devtools::install_github("satijalab/azimuth")
# need more than 40G

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyseurat)
library(tidyverse)
A316 <- readRDS(file = here::here("analysis","data_objects","05_clustering","A316_final.rds"))

DefaultAssay(A316) <- "RNA"
A316 <- RunAzimuth(A316, reference = "pbmcref", do.adt = TRUE)
saveRDS(A316,file = here::here("analysis","data_objects","05_clustering","A316_final_azimuth.rds"))

p1 <- DimPlot(A316, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, reduction = "wnn.umap") 
p2 <- DimPlot(A316, group.by = "broad.celltype", reduction = "wnn.umap")

pdf(file = here::here("analysis","plots","04_HA","azimuth_pbmc_predicted_celltype.pdf"), width = 36)
p1 + p2 
dev.off()

table(A316$predicted.celltype.l2)
#    B intermediate          B memory           B naive CD4 Proliferating 
#             22752             13270              8984                86 
#  NK Proliferating       Plasmablast          Platelet 
#                 1             19775                 7 
A316.sub <- A316 %>% filter(predicted.celltype.l2 %in% c("B intermediate","B memory","B naive","Plasmablast"))

p3 <- DimPlot(A316.sub, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
p4 <- DimPlot(A316.sub, group.by = "broad.celltype")
pdf(file = here::here("analysis","plots","05_clustering","azimuth_pbmc_predicted_celltype_subset.pdf"), width = 36)
p3 + p4
dev.off()

table(A316.sub$broad.celltype)

A316.remove <- A316 %>% filter(predicted.celltype.l2 %in% c("CD4 Proliferating","NK Proliferating","Platelet"))
p5 <- DimPlot(A316.remove, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
p6 <- DimPlot(A316.remove, group.by = "broad.celltype")
pdf(file = here::here("analysis","plots","04_HA","azimuth_pbmc_predicted_celltype_remove.pdf"), width = 36)
p5 + p6
dev.off()

#create table of removed cells 
removed.tab <- A316.remove %>% join_features(features = rownames(A316.remove@assays$Probes)) %>% pivot_wider(names_from = .feature, values_from = .abundance_Probes) %>% select(c(1:95,200:205))
write.csv(removed.tab,file = here::here("analysis","data_objects","04_HA","azimuth_removed_cells.csv"))

table(A316.remove$broad.celltype)

#    MemBCell PlasmaBlast 
#          13          81 

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()