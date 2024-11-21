library(Seurat)
library(sessioninfo)
library(harmony)
library(tidyseurat)
library(tidyverse)
library(ggridges)

###################################### Clustering based on Prot and RNA transcripts #################################
A316 <- readRDS(file = here::here("analysis","data_objects","04_HA","A316_specificity.rds"))

#have to add in Specificities that Sarah adjusted 
spec.sfa <- read.csv(here::here("analysis","230914_files_SFA","230914_Final","specificity.SFA.csv"),row.names = 1)
A316 <- A316 %>% left_join(spec.sfa,by = "CELL")

## Add Broad specificity column 
A316 <- A316 %>% mutate(Spec.Broad = case_when(specificity.SFA == "H2H1 Head" ~ "H2_Cross", 
                                                                specificity.SFA == "H2H3st" ~"H2_Cross",
                                                                specificity.SFA == "H2head" ~ "H2_Only",
                                                                specificity.SFA == "H2st" ~ "H2_Cross",
                                                                specificity.SFA == "H2H5 Head" ~ "H2_Cross",
                                                                specificity.SFA == "H2neg" ~ "Neg",
                                                                TRUE ~ "unclear"))
table(A316$Spec.Broad)
# H2_Cross  H2_Only      Neg  unclear 
#    17704    37222    18707       32 

## Cluster on Just H2 specific MemBcells first
H2.M <- A316 %>% filter(specificity.SFA %in% c("H2H1 Head", "H2H3st", "H2head", "H2st", "H2H5 Head") & Timepoint != "HA-" & broad.celltype == "MemBCell")

table(H2.M$broad.celltype)

# MemBCell 
#    35252 

## Did clustering in integration script, but will re-do it here because we dropped junk cells in 04_HA. Also want to drop HA- cells before clustering
DefaultAssay(H2.M) <- "RNA"
H2.M <- NormalizeData(H2.M, normalization.method = "LogNormalize")
H2.M <- FindVariableFeatures(H2.M, selection.method = "vst", nfeatures = 600)
VariableFeatures(H2.M) <- grep("IG[HKL]V|TRBV", VariableFeatures(H2.M), invert = TRUE, value = TRUE)
H2.M <- ScaleData(H2.M)
H2.M <- RunPCA(H2.M, features = VariableFeatures(H2.M))
H2.M <- RunHarmony(H2.M, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")

#make umaps of both PCs and harmony corrected PCs
H2.M <- RunUMAP(H2.M, reduction = "pca", reduction.name = "rna.umap", dims = 1:30, assay = "RNA")
H2.M <- RunUMAP(H2.M, reduction = "harmony", reduction.name = "harmony.rna.umap", dims = 1:30, assay = "RNA")

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_rna_UMAP_MemB.pdf"))
VariableFeaturePlot(H2.M)
ElbowPlot(H2.M)
DimPlot(H2.M, reduction = "rna.umap",group.by = "Subject")
DimPlot(H2.M, reduction = "harmony.rna.umap",group.by = "Subject")
DimPlot(H2.M, reduction = "rna.umap",group.by = "run")
DimPlot(H2.M, reduction = "harmony.rna.umap",group.by = "run")
DimPlot(H2.M, reduction = "rna.umap", label = TRUE)
DimPlot(H2.M, reduction = "rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(H2.M, reduction = "harmony.rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(H2.M, reduction = "rna.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

# Run on Prot individually #
DefaultAssay(H2.M) <- "Prot"
ProtFeatures <- rownames(H2.M@assays$Prot@data)
VariableProtFeatures <- ProtFeatures[- c(47, 48, 53, 56, 58)]
H2.M <- ScaleData(H2.M, features = VariableProtFeatures)
H2.M <- RunPCA(H2.M, assay = "Prot", slot = "data", features = rownames(H2.M), reduction.name = "apca")
H2.M <- RunHarmony(H2.M, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")

#make umaps of both PCs and harmony corrected PCs
H2.M <- RunUMAP(H2.M, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
H2.M <- RunUMAP(H2.M, reduction = "harmony.prot", dims = 1:18, assay = "Prot", reduction.name = "harmony.prot.umap", reduction.key = "har_protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_prot_UMAP_MemB.pdf"))
ElbowPlot(H2.M)
DimPlot(H2.M, reduction = "prot.umap",group.by = "Subject")
DimPlot(H2.M, reduction = "harmony.prot.umap",group.by = "Subject")
DimPlot(H2.M, reduction = "prot.umap",group.by = "run")
DimPlot(H2.M, reduction = "harmony.prot.umap",group.by = "run")
DimPlot(H2.M, reduction = "prot.umap", group.by = "orig.ident")
DimPlot(H2.M, reduction = "harmony.prot.umap", group.by = "orig.ident")
DimPlot(H2.M, reduction = "harmony.prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
dev.off()

# Do multi-modal clustering
H2.M <- FindMultiModalNeighbors(H2.M, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:15, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn")
H2.M <- FindClusters(H2.M, graph.name = "harmony.snn", algorithm = 3, resolution = 0.4, verbose = FALSE)
H2.M <- RunUMAP(H2.M, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_UMAP_MemB.pdf"))
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE) + theme(aspect.ratio = 1)
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "broad.celltype", label = TRUE)
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "run", split.by = "run", ncol = 2)
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "orig.ident")
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "Subject")
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_UMAP_clustering_MemB.pdf"))
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE)
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Subject", ncol = 4)
DimPlot(H2.M, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "run", ncol = 2)
dev.off()

table(H2.M$harmony.snn_res.0.4)

#    0    1   10    2    3    4    5    6    7    8    9 
# 7693 6318  121 5280 4945 4619 2543 1538 1331  448  416 

saveRDS(H2.M, file = here::here("analysis","data_objects","05_clustering","H2_MemB_clusters.rds"))

# look at diff exp protein markers all clusters with res 0.4 #
Prot <- data.frame(t(as.matrix(H2.M@assays$Prot@data)))
Prot$CELL <- rownames(Prot)
head(Prot)
meta <- H2.M@meta.data[, c("CELL", "harmony.snn_res.0.4", "vac_grp", "specificity.SFA", "Timepoint", "Subject", "c_call", "broad.celltype")]
Prot.meta <- left_join(Prot, meta, by = "CELL")

############# DSb analysis PDF #############
dsb.analysis <- Prot.meta[,c(1:62)] %>% filter(!(harmony.snn_res.0.4 %in% c("10")))
Proteins <- colnames(dsb.analysis)[1:60]

dsb.plot.list <- list()

pdf(file = here::here("analysis","plots","05_clustering","Protein_ridge_plots.pdf"))
for ( i in 1:length(Proteins)) {
  
  dsb.plot.list[[i]] <- dsb.analysis %>%
    ggplot(aes_string(x = Proteins[i], y = "harmony.snn_res.0.4")) + 
    geom_density_ridges(alpha = 0.6) + 
    theme_classic()
  
  plot(dsb.plot.list[[i]])
  
}
dev.off()

##Now remove some bad clusters and cluster only on Memory B Cells
# get rid of ("6","7", "9", "10")
H2.M.G <- H2.M %>% filter(!(harmony.snn_res.0.4 %in% c("6","7", "9", "10")))

## Re-do clustering with just memory cells
DefaultAssay(H2.M.G) <- "RNA"
H2.M.G <- NormalizeData(H2.M.G, normalization.method = "LogNormalize")
H2.M.G <- FindVariableFeatures(H2.M.G, selection.method = "vst", nfeatures = 600)
VariableFeatures(H2.M.G) <- grep("IG[HKL]V|TRBV", VariableFeatures(H2.M.G), invert = TRUE, value = TRUE)
H2.M.G <- ScaleData(H2.M.G)
H2.M.G <- RunPCA(H2.M.G, features = VariableFeatures(H2.M.G))
H2.M.G <- RunHarmony(H2.M.G, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")

H2.M.G <- RunUMAP(H2.M.G, reduction = "pca", reduction.name = "rna.umap", dims = 1:30, assay = "RNA")
H2.M.G <- RunUMAP(H2.M.G, reduction = "harmony", reduction.name = "harmony.rna.umap", dims = 1:30, assay = "RNA")

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_rna_UMAP_memB_G.pdf"))
VariableFeaturePlot(H2.M.G)
ElbowPlot(H2.M.G)
DimPlot(H2.M.G, reduction = "rna.umap",group.by = "Subject")
DimPlot(H2.M.G, reduction = "harmony.rna.umap",group.by = "Subject")
DimPlot(H2.M.G, reduction = "rna.umap",group.by = "run")
DimPlot(H2.M.G, reduction = "harmony.rna.umap",group.by = "run")
DimPlot(H2.M.G, reduction = "rna.umap", label = TRUE)
DimPlot(H2.M.G, reduction = "rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(H2.M.G, reduction = "harmony.rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(H2.M.G, reduction = "rna.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
DimPlot(H2.M.G, reduction = "harmony.rna.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

DefaultAssay(H2.M.G) <- "Prot"
ProtFeatures <- rownames(H2.M.G@assays$Prot@data)
VariableProtFeatures <- ProtFeatures[- c(47, 48, 53, 56, 58)]
H2.M.G <- ScaleData(H2.M.G, features = VariableProtFeatures)
H2.M.G <- RunPCA(H2.M.G, assay = "Prot", slot = "data", features = rownames(H2.M.G), reduction.name = "apca")
H2.M.G <- RunHarmony(H2.M.G, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")

#make umaps of both PCs and harmony corrected PCs
H2.M.G <- RunUMAP(H2.M.G, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
H2.M.G <- RunUMAP(H2.M.G, reduction = "harmony.prot", dims = 1:18, assay = "Prot", reduction.name = "harmony.prot.umap", reduction.key = "har_protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_prot_UMAP_memB_G.pdf"))
ElbowPlot(H2.M.G)
DimPlot(H2.M.G, reduction = "prot.umap",group.by = "Subject")
DimPlot(H2.M.G, reduction = "harmony.prot.umap",group.by = "Subject")
DimPlot(H2.M.G, reduction = "prot.umap",group.by = "run")
DimPlot(H2.M.G, reduction = "harmony.prot.umap",group.by = "run")
DimPlot(H2.M.G, reduction = "prot.umap", group.by = "orig.ident")
DimPlot(H2.M.G, reduction = "harmony.prot.umap", group.by = "orig.ident")
DimPlot(H2.M.G, reduction = "harmony.prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
dev.off()

# Combine Prot and RNA clusters #
H2.M.G <- FindMultiModalNeighbors(H2.M.G, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:15, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn")
H2.M.G <- FindClusters(H2.M.G, graph.name = "harmony.snn", algorithm = 3, resolution = 0.4, verbose = FALSE)

H2.M.G <- RunUMAP(H2.M.G, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

table(H2.M.G$harmony.snn_res.0.4)
#    0    1    2    3    4    5    6    7    8    9 
# 7614 6010 5456 4951 3699 2479 1113  407  115    2 

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_UMAP_memB_G.pdf"))
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "broad.celltype", label = TRUE)
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "run", split.by = "run", ncol = 2)
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "orig.ident")
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "Subject")
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","DimPlot_UMAP_clustering_memB_G.pdf"))
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE) + theme(aspect.ratio = 1)
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE)
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Subject", ncol = 4)
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "run", ncol = 2)
dev.off()

saveRDS(H2.M.G, file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))
H2.M.G <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))
### put clusters 7, 8, 9 into cluster 0

H2.M.G$harmony.snn_res.0.4 <- H2.M.G$harmony.snn_res.0.4 %>% recode(`7` = "0", `8` = "0", `9` = "0")
# H2.M.G <- H2.M.G %>% mutate(adj.harmony.snn_res.0.4 = case_when(harmony.snn_res.0.4 == "7" ~ "0", 
#                                                                 harmony.snn_res.0.4 == "8" ~ "0",
#                                                                 harmony.snn_res.0.4 == "9" ~ "0",
#                                                                 TRUE ~ TRUE))
#case_match(char_vec, "a" ~ "Apple", "b" ~ "Banana", .default = char_vec)

H2.M.G@meta.data <- mutate(H2.M.G@meta.data, Cluster.res0.4 = case_when(harmony.snn_res.0.4 == "0" ~ "AM2_0", 
                                                                        harmony.snn_res.0.4 == "1" ~ "AM3_acute_1",
                                                                        harmony.snn_res.0.4 == "2" ~ "RM_2",
                                                                        harmony.snn_res.0.4 == "3" ~ "AM3_chronic_3",
                                                                        harmony.snn_res.0.4 == "4" ~ "AM1a_4",
                                                                        harmony.snn_res.0.4 == "6" ~ "AM1a_6",
                                                                        harmony.snn_res.0.4 == "5" ~ "IgM_5",
                                                                        harmony.snn_res.0.4 == "7" ~ "AM2_0",
                                                                        TRUE ~ "unclear"))


# plyr::count(H2.M.G@meta.data$Cluster.res0.4)

saveRDS(H2.M.G, file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))
#H2.M.G <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))

DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "Cluster.res0.4", split.by = "Timepoint", ncol = 3)
DimPlot(H2.M.G, reduction = "harmony.wnn.umap", group.by = "Cluster.res0.4")

H2.M.G@meta.data %>% ggplot(aes(x = Timepoint, y = specificity.SFA, fill = specificity.SFA))+
  geom_bar(stat="identity") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()