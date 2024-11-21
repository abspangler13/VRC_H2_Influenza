# Load package
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(RColorBrewer)
library(sessioninfo)
library(tidyseurat)

# 0. load data --------------------------------
#seu.obj <- readRDS(file = here::here("analysis","data_objects","04_HA","A316_specificity.rds"))
seu.obj <- readRDS(file = here::here("analysis","data_objects","05_clustering", "H2_Mem_G_clusters.rds"))
Idents(seu.obj) <- "Cluster.res0.4"

seu.obj <- seu.obj %>% filter(vac_grp == "4B")

DefaultAssay(seu.obj) <- "RNA"
# 1. Wrap from Seurat -------------------------     

cds <- as.cell_data_set(seu.obj) 

# 2. Assign UMAP coordinate - cell embeddings -------
#move already made batch corrected umap into the UMAP slot so that cluster_cells() can use it. 
cds@int_colData@listData$reducedDims@listData$UMAP <- cds@int_colData@listData$reducedDims$HARMONY.WNN.UMAP

# 1. Cluster cells (using clustering info from seurat's UMAP) ---------------------------

## Assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

## Assign the cluster info 
list_cluster <- seu.obj@active.ident ## change object
#list_cluster
cds@clusters$UMAP$clusters <- list_cluster

raw_cds <- cds
saveRDS(raw_cds, file = here::here("analysis","data_objects","07_pseudotime","raw_cds_seurat_clusters_4B.rds"))

# Plot
pdf(file = here::here("analysis","plots","07_pseudotime","clusters_4B.pdf"))
plot_cells(cds,
  color_cells_by = 'Cluster.res0.4',
  label_groups_by_cluster = FALSE,
  group_label_size = 5) +
theme(legend.position = "right")
dev.off()

# 3. Learn trajectory graph ------------------------
cds <- learn_graph(cds)

saveRDS(cds, file = here::here("analysis","data_objects","07_pseudotime","learn_graph_cds_seurat_clusters_4B.rds"))
#cds <- readRDS(file = here::here("analysis","data_objects","07_pseudotime","learn_graph_cds_seurat_clusters.rds"))

pdf(file = here::here("analysis","plots","07_pseudotime","trajectory_graph_4B.pdf"))
plot_cells(cds,
           color_cells_by = 'Cluster.res0.4',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE,
           group_label_size = 5)
dev.off()

# 4. Order the cells in pseudotime -------------------
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == "RM_2"])) #Naive cluster selected as root


# cells ordered by monocle3 pseudotime
cds$pseudotime_RM_2 <- pseudotime(cds) #extract pseudo time from cds object
data.pseudo <- as.data.frame(colData(cds))

pdf(file = here::here("analysis","plots","07_pseudotime","pseudotime_RM_2_4B.pdf"))
plot_cells(cds[,which(cds$Cluster.res0.4 == "RM_2")],
           color_cells_by = 'Cluster.res0.4', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE) +
           theme(legend.position = "right") +
           ggtitle("RM_2")

plot_cells(cds,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           trajectory_graph_color = "black",
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE) +
           theme(legend.position = "right") 

ggplot(data.pseudo, aes(pseudotime_RM_2, reorder(Cluster.res0.4, pseudotime_RM_2,median), fill = Cluster.res0.4)) +
  geom_boxplot() +
  theme_bw()

dev.off()



# cells ordered by monocle3 pseudotime
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == "IgM_5"]))
cds$pseudotime_IgM_5 <- pseudotime(cds) #extract pseudo time from cds object
data.pseudo <- as.data.frame(colData(cds))

pdf(file = here::here("analysis","plots","07_pseudotime","pseudotime_IgM_5_4B.pdf"))
plot_cells(cds[,which(cds$Cluster.res0.4 == "IgM_5")],
           color_cells_by = 'Cluster.res0.4', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE) +
           theme(legend.position = "right") +
           ggtitle("IgM_5")

plot_cells(cds,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           trajectory_graph_color = "black",
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE) +
           theme(legend.position = "right") 
ggplot(data.pseudo, aes(pseudotime_IgM_5, reorder(Cluster.res0.4, pseudotime_IgM_5,median), fill = Cluster.res0.4)) +
  geom_boxplot() +
  theme_bw()

  dev.off()



cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == "AM1a_4"]))


# cells ordered by monocle3 pseudotime
cds$pseudotime_AM1a_4 <- pseudotime(cds) #extract pseudo time from cds object
data.pseudo <- as.data.frame(colData(cds))

pdf(file = here::here("analysis","plots","07_pseudotime","pseudotime_AM1a_4_4B.pdf"))
plot_cells(cds[,which(cds$Cluster.res0.4 == "AM1a_4")],
           color_cells_by = 'Cluster.res0.4', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE) +
           theme(legend.position = "right") +
           ggtitle("AM1a_4")

plot_cells(cds,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           trajectory_graph_color = "black",
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE) +
           theme(legend.position = "right")
ggplot(data.pseudo, aes(pseudotime_AM1a_4, reorder(Cluster.res0.4, pseudotime_AM1a_4 ,median), fill = Cluster.res0.4)) +
  geom_boxplot() +
  theme_bw()

dev.off()

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == "AM3_chronic_3"]))

# cells ordered by monocle3 pseudotime
cds$pseudotime_AM3_chronic_3 <- pseudotime(cds) #extract pseudo time from cds object
data.pseudo <- as.data.frame(colData(cds))

pdf(file = here::here("analysis","plots","07_pseudotime","pseudotime_AM3_chronic_3_4B.pdf"))
plot_cells(cds[,which(cds$Cluster.res0.4 == "AM3_chronic_3")],
           color_cells_by = 'Cluster.res0.4', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE) +
           theme(legend.position = "right") +
           ggtitle("AM3_chronic_3")

plot_cells(cds,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           trajectory_graph_color = "black",
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE) +
           theme(legend.position = "right")
ggplot(data.pseudo, aes(pseudotime_AM3_chronic_3, reorder(Cluster.res0.4, pseudotime_AM3_chronic_3 ,median), fill = Cluster.res0.4)) +
  geom_boxplot() +
  theme_bw()
dev.off()

saveRDS(cds, file = here::here("analysis","data_objects","07_pseudotime","cds_pseudotime_4B.rds"))
# cds <- readRDS(file = here::here("analysis","data_objects","07_pseudotime","cds_rna_seurat_clusters.rds"))


# 5. Finding genes that change as a function of pseudotime --------------------
# deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

# deg_bcells %>% 
#   arrange(q_value) %>% 
#   filter(status == 'OK') %>% 
#   head()

# FeaturePlot(seu.obj, features = c('E2F2', 'STMN1', 'CD52'))


# # visualizing pseudotime in seurat

# seu.obj$pseudotime <- pseudotime(cds)
# Idents(seu.obj) <- seu.obj$redefined_cluster
# FeaturePlot(seu.obj, features = "pseudotime", label = T)

#------------------- END

##### Filter object 

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()