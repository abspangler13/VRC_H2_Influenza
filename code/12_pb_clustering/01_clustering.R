# adapted from 240319_AnalyzePB_for_Abby.R from Sarah Andrews

# Required Packages #

library(plyr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(shazam)
library(alakazam)
library(tidyr)
library(gplots)
library(tidyseurat)
library(ggridges)
library(harmony)


# Object with all H2+ Cells #
A316.H2 <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_clusters.rds"))


# All PB #
H2.P <- A316.H2 %>% filter(broad.celltype == "PlasmaBlast")
plyr::count(H2.P@meta.data$specificity.SFA)
#           x  freq
# 1 H2H1 Head  2047
# 2    H2H3st    55
# 3 H2H5 Head   287
# 4    H2head 13819
# 5      H2st  2885


##### Do only PB without preclustering #####

# start over on Prot and RNA cluster with H2 Plasmablast only #

DefaultAssay(H2.P) <- "RNA"

H2.P <- NormalizeData(H2.P, normalization.method = "LogNormalize")
pdf(file = here::here("analysis","plots","12_pb_clustering","varible_feature_plot.pdf"))
VariableFeaturePlot(H2.P)
dev.off()

H2.P <- FindVariableFeatures(H2.P, selection.method = "vst", nfeatures = 600)
VariableFeatures(H2.P) <- grep("IG[HKL]V|TRBV", VariableFeatures(H2.P), invert = TRUE, value = TRUE)

H2.P <- ScaleData(H2.P)

H2.P <- RunPCA(H2.P, features = VariableFeatures(H2.P))
H2.P <- RunHarmony(H2.P, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")

pdf(file = here::here("analysis","plots","12_pb_clustering","elbow_plot.pdf"))
ElbowPlot(H2.P)
dev.off()

DefaultAssay(H2.P) <- "Prot"
ProtFeatures <- rownames(H2.P@assays$Prot@data)
ProtFeatures
VariableProtFeatures <- ProtFeatures[- c(47, 48, 53, 56, 58)]
VariableProtFeatures
H2.P <- ScaleData(H2.P, features = VariableProtFeatures)


H2.P <- RunPCA(H2.P, assay = "Prot", slot = "data", features = rownames(H2.P), reduction.name = "apca")
pdf(file = here::here("analysis","plots","12_pb_clustering","prot_elbow_plot.pdf"))
ElbowPlot(H2.P)
dev.off()
H2.P <- RunHarmony(H2.P, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")

# Combine Prot and RNA clusters #
H2.P <- FindMultiModalNeighbors(H2.P, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:15, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn")

H2.P <- RunUMAP(H2.P, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
H2.P <- FindClusters(H2.P, graph.name = "harmony.snn", algorithm = 3, resolution = 0.2, verbose = FALSE)

pdf(file = here::here("analysis","plots","12_pb_clustering","multi_modal_clusters_dimplot.pdf"))
DimPlot(H2.P, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.2", label = TRUE) + theme(aspect.ratio = 1)
DimPlot(H2.P, reduction = "harmony.wnn.umap", split.by = "Timepoint", label = TRUE, ncol = 3) + theme(aspect.ratio = 1)
DimPlot(H2.P, reduction = "harmony.wnn.umap", split.by = "clone.cell.type", label = TRUE) + theme(aspect.ratio = 1)
dev.off()

##### All clusters but 6 differ at transcript, not protein level, 6 is primarily Subject 309 and pretty different ###

# Rename cluster meta data columne for merging with Memory Clusters #

H2.P@meta.data <- mutate(H2.P@meta.data, PB.Clusters = case_when(harmony.snn_res.0.2 == "0" ~ "PB_0", 
                                                                 harmony.snn_res.0.2 == "1" ~"PB_1",
                                                                 harmony.snn_res.0.2 == "2" ~ "PB_2",
                                                                 harmony.snn_res.0.2 == "3" ~ "PB_3",
                                                                 harmony.snn_res.0.2 == "4" ~ "PB_4",
                                                                 harmony.snn_res.0.2 == "5" ~ "PB_5M",
                                                                 harmony.snn_res.0.2 == "6" ~ "PB_6junk",
                                                                 harmony.snn_res.0.2 == "7" ~ "PB_7junk",
                                                                 TRUE ~ "unclear"))


plyr::count(H2.P@meta.data$PB.Clusters)
#          x freq
# 1     PB_0 7539
# 2     PB_1 5658
# 3     PB_2 2118
# 4     PB_3 1352
# 5     PB_4  992
# 6    PB_5M  920
# 7 PB_6junk  512
# 8 PB_7junk    2


##### Remove cluster 6 and 7 and recluster #####

H2P.G <- H2.P %>% filter(!harmony.snn_res.0.2 %in% c("6", "7"))
plyr::count(H2P.G@meta.data$harmony.snn_res.0.2)
#   x freq
# 1 0 7539
# 2 1 5658
# 3 2 2118
# 4 3 1352
# 5 4  992
# 6 5  920


DefaultAssay(H2P.G) <- "RNA"

H2P.G <- NormalizeData(H2P.G, normalization.method = "LogNormalize")
pdf(file = here::here("analysis","plots","12_pb_clustering","rna_varible_feature_plot2.pdf"))
VariableFeaturePlot(H2P.G)
dev.off()

H2P.G <- FindVariableFeatures(H2P.G, selection.method = "vst", nfeatures = 600)
VariableFeatures(H2P.G) <- grep("IG[HKL]V|TRBV", VariableFeatures(H2P.G), invert = TRUE, value = TRUE)

H2P.G <- ScaleData(H2P.G)

H2P.G <- RunPCA(H2P.G, features = VariableFeatures(H2P.G))
H2P.G <- RunHarmony(H2P.G, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")
pdf(file = here::here("analysis","plots","12_pb_clustering","rna_elbow_plot2.pdf"))
ElbowPlot(H2P.G)
dev.off()

DefaultAssay(H2P.G) <- "Prot"
ProtFeatures <- rownames(H2P.G@assays$Prot@data)
ProtFeatures
VariableProtFeatures <- ProtFeatures[- c(47, 48, 53, 56, 58)]
VariableProtFeatures
H2P.G <- ScaleData(H2P.G, features = VariableProtFeatures)


H2P.G <- RunPCA(H2P.G, assay = "Prot", slot = "data", features = rownames(H2P.G), reduction.name = "apca")
pdf(file = here::here("analysis","plots","12_pb_clustering","prot_elbow_plot2.pdf"))
ElbowPlot(H2P.G)
dev.off()
H2P.G <- RunHarmony(H2P.G, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")

# Combine Prot and RNA clusters #
H2P.G <- FindMultiModalNeighbors(H2P.G, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:15, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn")

H2P.G <- RunUMAP(H2P.G, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
H2P.G <- FindClusters(H2P.G, graph.name = "harmony.snn", algorithm = 3, resolution = 0.2, verbose = FALSE)
H2P.G <- FindClusters(H2P.G, graph.name = "harmony.snn", algorithm = 3, resolution = 0.1, verbose = FALSE)

pdf(file = here::here("analysis","plots","12_pb_clustering","multi_modal_clusters_dimplot2.pdf"))
DimPlot(H2P.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.1", label = TRUE) + theme(aspect.ratio = 1)
dev.off()


H2P.G@meta.data <- mutate(H2P.G@meta.data, PB.Clusters.res01 = case_when(harmony.snn_res.0.1 == "0" ~ "PB_0", 
                                                                 harmony.snn_res.0.1 == "1" ~"PB_1",
                                                                 harmony.snn_res.0.1 == "2" ~ "PB_2",
                                                                 harmony.snn_res.0.1 == "3" ~ "PB_3",
                                                                 TRUE ~ "unclear"))

pdf(file = here::here("analysis","plots","12_pb_clustering","multi_modal_clusters_dimplot2_ann.pdf"))
DimPlot(H2P.G, reduction = "harmony.wnn.umap", group.by = "PB.Clusters.res01", label = TRUE) + theme(aspect.ratio = 1)
dev.off()

saveRDS(H2P.G, file = here::here("analysis","data_objects","12_pb_clustering","H2P_G.rds"))
