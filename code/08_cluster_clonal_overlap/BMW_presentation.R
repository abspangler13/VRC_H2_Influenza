library(pheatmap)
library(tidyseurat)

seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))

cells.per.clone <- as.matrix(table(seurat@meta.data$clone_subject_id))
cells.per.clone <- as.data.frame(cells.per.clone)

cells.per.clone.small <- cells.per.clone %>% filter(V1 > 10) %>% select(V1)
hist(cells.per.clone.small$V1)
dim(cells.per.clone.small)
summary(cells.per.clone.small$V1)


dat <- seurat@meta.data %>% select(clone_subject_id,Cluster.res0.4)

dat.small <- dat %>% filter(clone_subject_id %in% rownames(cells.per.clone.small))

map <- table(dat.small$clone_subject_id,dat.small$Cluster.res0.4)

pheatmap(map,scale = "none")

seurat <- SetIdent(seurat, value = "Cluster.res0.4")
DimPlot(seurat, reduction = "harmony.prot.umap")

seurat.top.clones <- seurat %>% filter(clone_subject_id %in% topclones)

write.csv(table(seurat.top.clones$clone_subject_id,seurat.top.clones$Cluster.res0.4), file = here::here("analysis","data_objects","08_cluster_clonal_overlap","top_clones_cluster.csv"))
write.csv(table(seurat.top.clones$clone_subject_id,seurat.top.clones$vac_grp), file = here::here("analysis","data_objects","08_cluster_clonal_overlap","top_clones_vac_grp.csv"))

clone.clusters <- as.data.frame.matrix(table(seurat.top.clones$clone_subject_id,seurat.top.clones$Cluster.res0.4))
clone.clusters$top.cluster <- colnames(clone.clusters)[apply(clone.clusters,1,which.max)]
clone.clusters$clone_subject_id <- rownames(clone.clusters)

vac.grp <- seurat.top.clones %>% select(clone_subject_id,vac_grp) %>% filter(clone_subject_id %in% clone.clusters$clone_subject_id) %>% distinct(clone_subject_id, .keep_all = TRUE)

dat <- left_join(clone.clusters,vac.grp, by = "clone_subject_id")

summary <- as.data.frame.matrix(table(dat$top.cluster,dat$vac_grp))

pheatmap(summary)