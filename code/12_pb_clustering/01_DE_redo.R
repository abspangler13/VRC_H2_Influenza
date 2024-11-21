library(Seurat)
library(tidyverse)
library(pheatmap)
library(tidyseurat) 
library(sessioninfo)
library(scales)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 

seurat <- readRDS(file = here::here("analysis","data_objects","12_pb_clustering", "H2P_G.rds"))
seurat <- SetIdent(seurat, value = "PB.Clusters.res01")

DefaultAssay(seurat) <- "RNA"
seurat <- ScaleData(seurat, features <- rownames(seurat))

#de.genes <- FindAllMarkers(object = seurat, logfc.threshold = 0, return.thresh = 1)
## from Chris and Chaim 
## https://teams.microsoft.com/l/message/19:d67347066e214b548342b7908ef9982a@thread.v2/1699378747854?context=%7B%22contextType%22%3A%22chat%22%7D
de.genes <- FindAllMarkers(object = seurat, min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1) 
saveRDS(de.genes,file = here::here("analysis","data_objects","12_pb_clustering","pb_de_rna_all_genes.rds"))
#de.genes <- readRDS(file = here::here("analysis","data_objects","12_pb_clustering","pb_de_rna_all_genes.rds"))

pdf(file = here::here("analysis","plots","12_pb_clustering","de_rna_QC.pdf"))
hist(de.genes$p_val)
hist(de.genes$p_val_adj)
hist(de.genes$avg_log2FC)
dev.off()

## pull out enriched 
enriched_genes <- de.genes %>% filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
features <- enriched_genes$gene

## look at number of enriched genes per cluster
table(enriched_genes$cluster)
# PB_2 PB_0 PB_1 PB_3 
#  363  147 1792 1281 

## look at number of cells per cluster
table(seurat$PB.Clusters.res01)
#  PB_0  PB_1  PB_2  PB_3 
# 13318  3011  1402   848

## see if any genes are enriched in multiple clusters
which(table(enriched_genes$gene, enriched_genes$cluster) > 1)
# integer(0) none of the genes are differentially expressed in multiple clusters

### remove ribo and histone genes before taking the top 30
ribo <- rownames(de.genes)[grep("^RP",de.genes$gene)]
histo <- rownames(de.genes)[grep("^HIST",de.genes$gene)]
de.genes <- de.genes %>% filter(!gene %in% ribo)
de.genes <- de.genes %>% filter(!gene %in% histo)

###save enriched genes without histo and ribo
write.csv(de.genes, file = here::here("analysis","data_objects","12_pb_clustering", "enriched_rna_de.csv"))

##too many differentially expressed genes so find top 10 for each cluster
markers.seurat <- de.genes %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 30) %>% slice_max(order_by = avg_log2FC, n = 30)
table(markers.seurat$cluster)
# PB_2 PB_0 PB_1 PB_3 
#   30   30   30   30 

features <- markers.seurat$gene

pdf(file = here::here("analysis","plots","12_pb_clustering","enriched_genes.pdf"), height = 15)
DoHeatmap(subset(seurat,downsample = 1000), features = features)
dev.off()


#### Protein #####
DefaultAssay(seurat) <- "Prot"
seurat <- ScaleData(seurat, features <- rownames(seurat))
de.prot <- FindAllMarkers(object = seurat, min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1) 
saveRDS(de.prot,file = here::here("analysis","data_objects","12_pb_clustering","pb_de_prot_all_genes.rds"))

summary(de.prot$avg_log2FC)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -36.7722  -4.8811  -0.9339  -2.3187   1.0714  17.9641 

pdf(file = here::here("analysis","plots","12_pb_clustering","de_prot_QC.pdf"))
hist(de.prot$p_val)
hist(de.prot$p_val_adj)
hist(de.prot$avg_log2FC)
dev.off()

enriched_prot <- de.prot %>% filter(avg_log2FC > 0.5 & p_val_adj < 0.05)

# How many enriched proteins per cluster
table(enriched_prot$cluster)
# PB_2 PB_0 PB_1 PB_3 
#    2   16   11   19 

which(table(enriched_prot$gene, enriched_prot$cluster) > 1)
# integer(0)

features <- enriched_prot$gene
pdf(file = here::here("analysis","plots","12_pb_clustering","enriched_prots.pdf"))
DoHeatmap(subset(seurat,downsample = 1000), features = features)
dev.off()


# ### annotated dotplot for RNA
# more_enriched_genes <- enriched_genes %>% filter(avg_log2FC > 0.8)
# pdf(file = here::here("analysis","plots","06_DE","rna_dotplot.pdf"),width = 16)
# DotPlot(seurat, assay = "RNA", features = unique(more_enriched_genes$gene)) + theme(axis.text.x = element_text(angle = 90))
# dev.off()


# #### Top fc violing plots for RNA
# top.fc <- enriched_genes %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 2)
# pdf(file = here::here("analysis","plots","06_DE","topfc_vlnplot.pdf"),width = 16)
# VlnPlot(seurat, assay = "RNA", features = top.fc$gene, ncol = 5, pt.size = 0)
# dev.off()

pdf(file = here::here("analysis","plots","12_pb_clustering","VlnPlot_nCount_RNA.pdf"))
VlnPlot(seurat, features = "nCount_RNA")
dev.off()

pdf(file = here::here("analysis","plots","12_pb_clustering","VlnPlot_nCount_RNA.pdf"))