library(Seurat)
library(tidyseurat)
library(pheatmap)
library(patchwork)
library(ggplotify)
library(ggplot2)
library(sessioninfo)

#setwd("/Users/spanglerab/Desktop/VRC316")
#seurat <- readRDS("230927_files_SFA/H2.M.G")
seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))

clusters <- unique(seurat$Cluster.res0.4)
mat <- matrix(, nrow = length(clusters), ncol = length(clusters))

for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    a <- seurat %>% filter(Cluster.res0.4 == clusters[i]) %>% select(clone_subject_id)
    b <- seurat %>% filter(Cluster.res0.4 == clusters[j]) %>% select(clone_subject_id)
    intr <- intersect(a$clone_subject_id,b$clone_subject_id)
    mat[i,j] <- dim(seurat %>% filter(Cluster.res0.4 == clusters[i] & clone_subject_id %in% intr))[2] / dim(seurat %>% filter(Cluster.res0.4 == clusters[i]))[2]
  }
}

dat <- as.data.frame(mat)
rownames(dat) <- clusters
colnames(dat) <- clusters
#dat <- read.csv(file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell.csv"),row.names = 1)
#pheatmap(dat,cluster_rows = FALSE, cluster_cols = FALSE)

reoder <- c("AM1a_4","AM1a_6","AM2_0","AM3_acute_1","AM3_chronic_3","IgM_5","RM_2")
dat <- dat[reoder ,reoder ]

write.csv(dat, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell.csv"))
#dat <- read.csv(file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell.csv"),row.names = 1)

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell.pdf"))
pheatmap(dat,display_numbers = round(dat,2), main = "All H2",cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(dat,display_numbers = round(dat,2))
dev.off()

p.all <- as.ggplot(pheatmap(dat,display_numbers = round(dat,2), main = "All",cluster_rows = FALSE, cluster_cols = FALSE))
p.all.clustered <- as.ggplot(pheatmap(dat,display_numbers = round(dat,2), cluster_cols = FALSE))


####### H2 Only #########

unique(seurat$Spec.Broad)
# [1] "H2_Cross" "H2_Only" 

h2.only <- seurat %>% filter(Spec.Broad == "H2_Only")

clusters <- unique(h2.only$Cluster.res0.4)
mat <- matrix(, nrow = length(clusters), ncol = length(clusters))

for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    a <- h2.only %>% filter(Cluster.res0.4 == clusters[i]) %>% select(clone_subject_id)
    b <- h2.only %>% filter(Cluster.res0.4 == clusters[j]) %>% select(clone_subject_id)
    intr <- intersect(a$clone_subject_id,b$clone_subject_id)
    mat[i,j] <- dim(h2.only %>% filter(Cluster.res0.4 == clusters[i] & clone_subject_id %in% intr))[2] / dim(h2.only %>% filter(Cluster.res0.4 == clusters[i]))[2]
  }
}

dat <- as.data.frame(mat)
rownames(dat) <- clusters
colnames(dat) <- clusters

#dat <- read.csv(file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_h2only.csv"),row.names = 1)
#pheatmap(dat,cluster_rows = FALSE, cluster_cols = FALSE)

reoder <- c("AM1a_4","AM1a_6","AM2_0","AM3_acute_1","AM3_chronic_3","IgM_5","RM_2")
dat <- dat[reoder ,reoder ]

write.csv(dat, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_h2only.csv"))
#dat <- read.csv(file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_h2only.csv"),row.names = 1)

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_h2only.pdf"))
pheatmap(dat,display_numbers = round(dat,2), main = "H2 Only",cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(dat,display_numbers = round(dat,2), main = "H2 Only")
dev.off()

p.h2.only <- as.ggplot(pheatmap(dat,display_numbers = round(dat,2), main = "H2 Only",cluster_rows = FALSE, cluster_cols = FALSE))
p.h2.only.clustered <- as.ggplot(pheatmap(dat,display_numbers = round(dat,2), main = "H2 Only"), cluster_cols = FALSE)

####### H2 Cross #########
h2.cross <- seurat %>% filter(Spec.Broad == "H2_Cross")

clusters <- unique(h2.cross$Cluster.res0.4)
mat <- matrix(, nrow = length(clusters), ncol = length(clusters))
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    a <- h2.cross %>% filter(Cluster.res0.4 == clusters[i]) %>% select(clone_subject_id)
    b <- h2.cross %>% filter(Cluster.res0.4 == clusters[j]) %>% select(clone_subject_id)
    intr <- intersect(a$clone_subject_id,b$clone_subject_id)
    mat[i,j] <- dim(h2.cross %>% filter(Cluster.res0.4 == clusters[i] & clone_subject_id %in% intr))[2] / dim(h2.cross %>% filter(Cluster.res0.4 == clusters[i]))[2]
  }
}

dat <- as.data.frame(mat)
rownames(dat) <- clusters
colnames(dat) <- clusters

#dat <- read.csv(file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_h2only.csv"),row.names = 1)
#pheatmap(dat,cluster_rows = FALSE, cluster_cols = FALSE)

reoder <- c("AM1a_4","AM1a_6","AM2_0","AM3_acute_1","AM3_chronic_3","IgM_5","RM_2")
dat <- dat[reoder ,reoder]

write.csv(dat, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_h2cross.csv"))

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_h2cross.pdf"))
pheatmap(dat,display_numbers = round(dat,2), main = "H2 Cross",cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(dat,display_numbers = round(dat,2), main = "H2 Cross", cluster_cols = FALSE)
dev.off()

p.h2.cross <- as.ggplot(pheatmap(dat,display_numbers = round(dat,2), main = "H2 Cross",cluster_rows = FALSE, cluster_cols = FALSE))
p.h2.cross.clustered <- as.ggplot(pheatmap(dat,display_numbers = round(dat,2), main = "H2 Cross"), cluster_cols = FALSE)

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_all.pdf"), width = 24)
p.all + p.h2.cross + p.h2.only
dev.off()

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","clonal_cluster_overlap_by_cell_clustered.pdf"), width = 24)
p.all.clustered + p.h2.cross.clustered + p.h2.only.clustered
dev.off()


####### H2_only by vac_grp #############
seurat.list <- SplitObject(seurat, split.by = "vac_grp")
vac.grps <- unique(seurat$vac_grp)
heatmaps <- list()
for(x in 1:length(vac.grps)){
  seurat <- seurat.list[[x]] %>% filter(Spec.Broad == "H2_Only",vac_grp == vac.grps[x])
  clusters <- unique(seurat$Cluster.res0.4)
  mat <- matrix(, nrow = length(clusters), ncol = length(clusters))

  for(i in 1:length(clusters)){
    for(j in 1:length(clusters)){
      a <- seurat %>% filter(Cluster.res0.4 == clusters[i]) %>% select(clone_subject_id)
      b <- seurat %>% filter(Cluster.res0.4 == clusters[j]) %>% select(clone_subject_id)
      intr <- intersect(a$clone_subject_id,b$clone_subject_id)
      mat[i,j] <- dim(seurat %>% filter(Cluster.res0.4 == clusters[i] & clone_subject_id %in% intr))[2] / dim(seurat %>% filter(Cluster.res0.4 == clusters[i]))[2]
    }
  }

  dat <- as.data.frame(mat)
  rownames(dat) <- clusters
  colnames(dat) <- clusters
  #reorder <- c("AM1a_4","AM1a_6","AM2_0","AM3_acute_1","AM3_chronic_3","IgM_5","RM_2")
  reorder <- sort(clusters)
  dat <- dat[reorder,reorder]
  write.csv(dat, file = here::here("analysis","data_objects","08_cluster_clonal_overlap",paste0("clonal_cluster_overlap_h2_only_",vac.grps[x],".csv")))
  heatmaps[[x]] <- as.ggplot(pheatmap(dat,display_numbers = round(dat,2), main = paste0("H2 Only ",vac.grps[x]),cluster_rows = FALSE, cluster_cols = FALSE))
}

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","clonal_cluster_overlap_h2_only_vacgrp.pdf"), width = 16)
patchwork::wrap_plots(heatmaps, 
                      nrow = 2, ncol = 2)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()