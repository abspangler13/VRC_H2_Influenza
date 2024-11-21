library(Seurat)
library(tidyseurat)
library(pheatmap)
library(cowplot)

A316 <- readRDS(file = here::here("analysis","data_objects","05_clustering","A316_final_azimuth.rds"))
#A316 <- readRDS(file = here::here("analysis","data_objects","05_clustering","A316_final.rds"))


### RNA
pdf(file = here::here("analysis","plots","05_clustering","dimplot_rna_subject.pdf"), width = 14)
A316 %>% DimPlot(split.by = "Subject", reduction = "rna.umap", group.by = "Subject", ncol=4)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","dimplot_rna_run.pdf"), width = 14)
A316 %>% DimPlot(split.by = "run", reduction = "rna.umap", group.by = "run", ncol=2)
A316 %>% DimPlot(split.by = "run", reduction = "rna.umap", group.by = "Subject", ncol=2)
dev.off()

A316 %>% DimPlot(split.by = "vac_grp", reduction = "rna.umap", ncol = 2, group.by = "vac_grp") 
A316 %>% DimPlot(split.by = "orig.ident", reduction = "rna.umap", ncol = 5, group.by = "orig.ident") 

A316 %>% filter(vac_grp %in% c("4A","4B")) %>% DimPlot(split.by = "Subject", reduction = "rna.umap", ncol = 4, group.by = "Subject")

####### Prot
pdf(file = here::here("analysis","plots","05_clustering","dimplot_prot_subject.pdf"), width = 14)
A316 %>% DimPlot(split.by = "Subject", reduction = "prot.umap", group.by = "Subject", ncol=4)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","dimplot_prot_run.pdf"), width = 14)
A316 %>% DimPlot(split.by = "run", reduction = "prot.umap", group.by = "run", ncol=2)
A316 %>% DimPlot(split.by = "run", reduction = "prot.umap", group.by = "Subject", ncol=2)
dev.off()
A316 %>% filter(vac_grp %in% c("4A","4B")) %>% DimPlot(split.by = "Subject",reduction = "prot.umap",group.by="Subject",ncol = 4)

pdf(file = here::here("analysis","plots","05_clustering","dimplot_prot_run.pdf"), width = 14)
A316 %>% DimPlot(split.by = "vac_grp", reduction = "prot.umap", group.by = "Subject", ncol=2)
dev.off()

##### WNN
pdf(file = here::here("analysis","plots","05_clustering","dimplot_wnn_subject.pdf"), width = 14)
A316 %>% DimPlot(split.by = "Subject", reduction = "wnn.umap", group.by = "Subject", ncol=4)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","dimplot_wnn_subject_cluster.pdf"), width = 14)
A316 %>% DimPlot(split.by = "Subject", reduction = "wnn.umap", ncol=4)
dev.off()

A316 %>% filter(vac_grp %in% c("4A","4B")) %>% DimPlot(split.by = "Subject",reduction = "wnn.umap",group.by="Subject",ncol = 4)


df <- table(A316$Subject,A316$wsnn_res.0.5)
write.csv(df,file = here::here("analysis","data_objects","05_clustering","subject_by_cluster_tab.csv"))

# A316 <- FindMultiModalNeighbors(A316, reduction.list = list("pca", "apca"), dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight")
A316 <- FindClusters(A316, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = FALSE)
A316 <- FindClusters(A316, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)

breaksList = seq(-4, 4, by = 1)

subject.cluster <- as.matrix(table(A316$Subject,A316$wsnn_res.0.2))
pheatmap(subject.cluster,cluster_rows=F, cluster_cols=F)
pheatmap(subject.cluster,scale = "row",cluster_rows=F, cluster_cols=F)
pheatmap(subject.cluster,scale = "column",cluster_rows=F, cluster_cols=F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdYlBu")))(8),breaks = seq(-4, 4, by = 1))

har.subject.cluster <- as.matrix(table(A316$Subject,A316$harmony.snn_res.0.2))
pheatmap(har.subject.cluster,cluster_rows=F, cluster_cols=F)
pheatmap(har.subject.cluster,scale = "row",cluster_rows=F, cluster_cols=F) 
pheatmap(har.subject.cluster,scale = "column",cluster_rows=F, cluster_cols=F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdYlBu")))(8), breaks = seq(-4, 4, by = 1))

vac_grp.cluster <- as.matrix(table(A316$vac_grp,A316$wsnn_res.0.2))
pheatmap(vac_grp.cluster,cluster_rows=F, cluster_cols=F)
pheatmap(vac_grp.cluster,scale = "row",cluster_rows=F, cluster_cols=F)
pheatmap(vac_grp.cluster,scale = "column",cluster_rows=F, cluster_cols=F)


run.cluster <- as.matrix(table(A316$run,A316$wsnn_res.0.2))
pheatmap(run.cluster,cluster_rows=F, cluster_cols=F)
pheatmap(run.cluster,scale = "row",cluster_rows=F, cluster_cols=F)
pheatmap(run.cluster,scale = "column",cluster_rows=F, cluster_cols=F)

A316 %>% DimPlot(split.by = "Subject", group.by = "wsnn_res.0.2", reduction = "wnn.umap", ncol=4)

saveRDS(A316,file = here::here("analysis","data_objects","05_clustering","A316_final_azimuth.rds"))


#### look at UMI counts for both RNA and Protein 
A316 %>% DimPlot(split.by = "Subject", reduction = "wnn.umap", ncol = 4, group.by = "wsnn_res.0.2")


DefaultAssay(A316) <- "Prot"
pdf(file = here::here("analysis","plots","05_clustering","UMI_plots.pdf"),width = 15)
p <- A316 %>% FeaturePlot(features = "log_nCount_Prot",reduction = "wnn.umap",split.by = "Subject",combine = FALSE)
plot_grid(plotlist = p, ncol = 5)
p2 <- A316 %>% FeaturePlot(features = "log_nCount_RNA",reduction = "wnn.umap",split.by = "Subject",combine = FALSE) 
plot_grid(plotlist = p2, ncol = 5)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","UMI_vlnplots.pdf"))
A316 %>% VlnPlot(features = "log_nCount_Prot",group.by = "Subject",pt.size = 0)
A316 %>% VlnPlot(features = "log_nCount_Prot",group.by = "run",pt.size = 0)
A316 %>% VlnPlot(features = "log_nCount_RNA",group.by = "Subject", pt.size = 0)
A316 %>% VlnPlot(features = "log_nCount_RNA",group.by = "run",pt.size = 0)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","feature_vlnplots.pdf"))
A316 %>% VlnPlot(features = "nFeature_Prot",group.by = "Subject",pt.size = 0)
A316 %>% VlnPlot(features = "nFeature_Prot",group.by = "run",pt.size = 0)
A316 %>% VlnPlot(features = "nFeature_RNA",group.by = "Subject", pt.size = 0)
A316 %>% VlnPlot(features = "nFeature_RNA",group.by = "run", pt.size = 0)
dev.off()

prot_umi_mean <- colMeans(as.matrix(A316@assays$Prot@data))
A316$prot_umi_mean <- prot_umi_mean
A316.drop %>% VlnPlot(features = "prot_umi_mean",group.by = "Subject",pt.size = 0)

RNA_umi_mean <- colMeans(as.matrix(A316@assays$RNA@data))
A316$RNA_umi_mean <- RNA_umi_mean
A316 %>% VlnPlot(features = "RNA_umi_mean",group.by = "Subject",pt.size = 0)


A316 %>% DimPlot(reduction = "wnn.umap",group.by = "Subject",split.by = "Subject",ncol =4 )


##### check batch corrected object 
A316.cor <- readRDS(file = "../Experiment_316_2/analysis/data_objects/05_clustering/A316_final.rds")
prot_umi_mean <- colMeans(as.matrix(A316.cor@assays$Prot@data))
A316.cor$prot_umi_mean <- prot_umi_mean
A316.cor %>% VlnPlot(features = "prot_umi_mean",group.by = "Subject",pt.size = 0)

RNA_umi_mean <- colMeans(as.matrix(A316.cor@assays$RNA@data))
A316.cor$RNA_umi_mean <- RNA_umi_mean
A316.cor %>% VlnPlot(features = "RNA_umi_mean",group.by = "Subject",pt.size = 0)

A316.cor %>% DimPlot(reduction = "harmony.rna.umap",group.by = "Subject",split.by = "Subject",ncol =4)

pdf(file = here::here("analysis","plots","05_clustering","UMI_vlnplots_cor.pdf"))
A316.cor %>% VlnPlot(features = "log_nCount_Prot",group.by = "Subject",pt.size = 0)
A316.cor %>% VlnPlot(features = "log_nCount_Prot",group.by = "run",pt.size = 0)
A316.cor %>% VlnPlot(features = "log_nCount_RNA",group.by = "Subject", pt.size = 0)
A316.cor %>% VlnPlot(features = "log_nCount_RNA",group.by = "run",pt.size = 0)
dev.off()

pdf(file = here::here("analysis","plots","05_clustering","feature_vlnplots_cor.pdf"))
A316.cor %>% VlnPlot(features = "nFeature_Prot",group.by = "Subject",pt.size = 0)
A316.cor %>% VlnPlot(features = "nFeature_Prot",group.by = "run",pt.size = 0)
A316.cor %>% VlnPlot(features = "nFeature_RNA",group.by = "Subject", pt.size = 0)
A316.cor %>% VlnPlot(features = "nFeature_RNA",group.by = "run", pt.size = 0)
dev.off()


summary(A316@meta.data$nCount_RNA)

A316 %>% DimPlot(reduction = "rna.umap",group.by = "Subject",split.by = "Subject",ncol =4 )

# We have clusters that expression non bcell markers. Cd14, 56,3,4. See how multimodal cluster looks after that
A316 %>% FeaturePlot(features = c("P-CD14","P-CD56","P-CD3","P-CD4"),reduction = "prot.umap")
A316 %>% VlnPlot(features = c("P-CD14","P-CD56","P-CD3","P-CD4"),group.by = "Subject", pt.size = 0,ncol=2)
A316.drop %>% VlnPlot(features = c("P-CD14","P-CD56","P-CD3","P-CD4"),group.by = "wsnn_res.0.2", pt.size = 0,ncol=2)
A316 %>% DimPlot(reduction = "prot.umap",group.by = "wsnn_res.0.2")
A316 %>% DimPlot(reduction = "prot.umap",group.by = "wsnn_res.0.2",split.by = "Subject",ncol = 4)
table(A316$Subject,A316$wsnn_res.0.2)

A316.drop %>% VlnPlot(features = c("CD14","CD56","CD3","CD4"),group.by = "wsnn_res.0.2", pt.size = 0,ncol=2,assay = "RNA",log = TRUE)

A316.cor %>% FeaturePlot(features = c("P-CD14","P-CD56","P-CD3","P-CD4"),reduction = "prot.umap")
A316.cor %>% VlnPlot(features = c("P-CD14","P-CD56","P-CD3","P-CD4"),group.by = "Subject", pt.size = 0,ncol=2)
A316.cor %>% VlnPlot(features = c("P-CD14","P-CD56","P-CD3","P-CD4"),group.by = "wsnn_res.0.2", pt.size = 0,ncol=2)
A316.cor %>% DimPlot(reduction = "prot.umap",group.by = "wsnn_res.0.2")
A316.drop %>% DimPlot(reduction = "prot.umap",group.by = "wsnn_res.0.2",split.by = "Subject",ncol = 4)
table(A316.cor$Subject,A316.cor$wsnn_res.0.2)

A316.drop %>% VlnPlot(features = "nCount_Prot",group.by = "wsnn_res.0.2", pt.size = 0)
A316.drop %>% VlnPlot(features = "log_nCount_Prot",group.by = "wsnn_res.0.2", pt.size = 0)

A316 %>% VlnPlot(features = "nCount_Prot",group.by = "wsnn_res.0.2", pt.size = 0)
A316 %>% VlnPlot(features = "log_nCount_Prot",group.by = "wsnn_res.0.2", pt.size = 0)
A316 %>% VlnPlot(features = "prot_umi_mean",group.by = "wsnn_res.0.2",pt.size = 0) 
A316.cor %>% VlnPlot(features = "prot_umi_mean",group.by = "wsnn_res.0.2",pt.size = 0)

A316 %>% VlnPlot(features = "nCount_Prot",group.by = "orig.ident", pt.size = 0)
A316 %>% VlnPlot(features = "log_nCount_Prot",group.by = "orig.ident", pt.size = 0)
A316 %>% VlnPlot(features = "prot_umi_mean",group.by = "orig.ident",pt.size = 0) 
A316.cor %>% VlnPlot(features = "prot_umi_mean",group.by = "orig.ident",pt.size = 0)



#### LISI https://github.com/immunogenomics/lisi
devtools::install_github("immunogenomics/lisi")
library(devtools)
library(lisi)

pcs <- A316@reductions

#### ADT QC and normaliziation
BiocManager::install("DropletUtils")
library(DropletUtils)

A316.good %>% VlnPlot(features = "nCount_Prot",group.by = "orig.ident",pt.size=0)