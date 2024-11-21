library(Seurat)
library(tidyverse)
library(tidyseurat) 
library(sessioninfo)
library(readxl)
library(pheatmap)

seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering", "H2_Mem_G_clusters.rds"))
seurat <- SetIdent(seurat, value = "Cluster.res0.4")

seurat %>% group_by(clone_subject_id) %>% summarise()

dat <- as.data.frame(table(seurat$clone_subject_id,seurat$Cluster.res0.4))

dat <- dat %>% pivot_wider(names_from = Var2, values_from = Freq)

AM1a_4 <- dat %>% filter(AM1a_4 > 0) %>% select(Var1)
AM1a_6 <- dat %>% filter(AM1a_6 > 0) %>% select(Var1)
AM2_0 <- dat %>% filter(AM2_0 > 0) %>% select(Var1)
AM3_acute_1 <- dat %>% filter(AM3_acute_1 > 0) %>% select(Var1)
AM3_chronic_3 <- dat %>% filter(AM3_chronic_3 > 0) %>% select(Var1)
IgM_5 <- dat %>% filter(IgM_5 > 0) %>% select(Var1)
RM_2 <- dat %>% filter(RM_2 > 0) %>% select(Var1)

nrow(intersect(RM_2,RM_2))

dat <- read_excel(path = here::here("analysis","data_objects","08_cluster_clonal_overlap","cluster_clonal_overlap.xlsx"), sheet = 2, col_names = TRUE)
dat <- as.data.frame(dat)
rownames(dat) <- dat$`...1`
dat <- dat[,-1]
pheatmap(dat,cluster_rows = FALSE, cluster_cols = FALSE)