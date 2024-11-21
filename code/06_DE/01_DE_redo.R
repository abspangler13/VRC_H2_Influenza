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

seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering", "H2_Mem_G_clusters.rds"))
seurat <- SetIdent(seurat, value = "Cluster.res0.4")

DefaultAssay(seurat) <- "RNA"
seurat <- ScaleData(seurat, features <- rownames(seurat))

all.markers <- readRDS(file = here::here("analysis","data_objects","06_DE",paste0("rna_markers_0.4_seurat.rds")))

enriched_genes <- all.markers %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 20) %>% slice_max(order_by = avg_log2FC, n = 20)
features <- enriched_genes$gene

table(enriched_genes$cluster)
        #  RM_2        AM1a_4 AM3_chronic_3   AM3_acute_1         AM2_0 
        #    15            80           111            98             4 
        # IgM_5        AM1a_6 
        #    76            74 
table(seurat$Cluster.res0.4)
    #    AM1a_4        AM1a_6         AM2_0   AM3_acute_1 AM3_chronic_3 
    #      3699          1113          8138          6010          4951 
    #     IgM_5          RM_2 
    #      2479          5456

table(enriched_genes$gene, enriched_genes$cluster)

pdf(file = here::here("analysis","plots","06_DE","enriched_genes.pdf"), height = 35)
DoHeatmap(subset(seurat,downsample = 1000), features = features)
dev.off()

deg.overlap <- as.data.frame.matrix(table(enriched_genes$gene, enriched_genes$cluster))

pdf(file = here::here("analysis","plots","06_DE","deg_overlap.pdf"), height = 35)
pheatmap(deg.overlap)
dev.off()

#### Protein #####
DefaultAssay(seurat) <- "Prot"
seurat <- ScaleData(seurat, features <- rownames(seurat))
prot.all.markers <- readRDS(file = here::here("analysis","data_objects","06_DE","prot_all_markers_0.4_seurat.rds"))

hist(prot.all.markers$avg_log2FC)

table(prot.all.markers$cluster)

enriched_prot <- prot.all.markers %>% filter(avg_log2FC > 0)

table(enriched_prot$cluster)

features <- enriched_prot$gene
pdf(file = here::here("analysis","plots","06_DE","enriched_prots.pdf"))
DoHeatmap(subset(seurat,downsample = 1000), features = features)
dev.off()

dep.overlap <- as.data.frame.matrix(table(enriched_prot$gene, enriched_prot$cluster))

pdf(file = here::here("analysis","plots","06_DE","dep_overlap.pdf"))
pheatmap(dep.overlap)
dev.off()


#### RNA ####

seurat_raw <-seurat
genes <- unique(enriched_genes$gene)
rna <- t(as.matrix(seurat@assays$RNA@scale.data[genes,]))
annotation <- seurat@meta.data %>% select(Cluster.res0.4)
rna <- cbind(rna,annotation$Cluster.res0.4)
#colnames(rna)[358] <- "Cluster"
colnames(rna)[103] <- "Cluster"
rna <- as.data.frame(rna)
#rna <- rna %>% mutate_at(c(1:357),as.numeric)
rna <- rna %>% mutate_at(c(1:102),as.numeric)
rna <- rna %>% arrange(Cluster)

df <- rna %>% group_by(Cluster) %>% summarise_at(genes,mean)

df.transposed <- df %>% pivot_longer(cols= -1) %>% pivot_wider(names_from = "Cluster",values_from = "value") %>% rename(gene = name)
df.transposed <- as.data.frame(df.transposed)
rownames(df.transposed) <- df.transposed$gene
df.transposed <- df.transposed[,-1]
df.transposed <- df.transposed[,c("RM_2","AM1a_4","AM3_chronic_3", "AM3_acute_1", "AM2_0","IgM_5","AM1a_6")]

pdf(file = here::here("analysis","plots","06_DE","RNA_pheatmap_small.pdf"),height = 20)
pheatmap(df.transposed,scale = "row", cluster_rows = FALSE, cluster_cols = FALSE) 
dev.off()


#### gene annotation to heatmap
#write.csv(table(enriched_genes$gene, enriched_genes$cluster), file = here::here("analysis","data_objects","06_DE","gene_ann_20.csv"))

#manually designate annotations in excel and then read it back in
gene_ann <- read.csv(file = here::here("analysis","data_objects","06_DE","gene_ann_20.csv"),row.names = 1)
gene_ann <- gene_ann[,c("ann1","ann2")] #,"ann3"

my_color_palette <- hue_pal()(7)
annoCol<-list(ann1=c(AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"),
ann2 = c(`NA` = "#FFFFFF" , AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"),
#ann3 = c(`NA` = "#FFFFFF" , AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"),
unique.dat.Cluster. = c(AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"))

# [1] "AM1a_4"        "AM1a_6"        "AM2_0"         "AM3_acute_1"  
# [5] "AM3_chronic_3" "IgM_5"         "RM_2"  
# [1] "#F8766D" "#C49A00" "#53B400" "#00C094" "#00B6EB" "#A58AFF" "#FB61D7"
col_ann <- data.frame(unique(rna$Cluster))
rownames(col_ann) <- unique(rna$Cluster) 

order_genes <- gene_ann %>% arrange(match(ann1, c("AM3_chronic_3","AM3_acute_1","AM1a_4","AM1a_6","IgM_5","RM_2","AM2_0")))
order_genes <- rownames(order_genes)
df.transposed <- df.transposed[order_genes,] ### problem line

df.transposed1 <- df.transposed[1:170,c("AM3_chronic_3","AM3_acute_1","AM1a_4","AM1a_6","IgM_5","RM_2","AM2_0")]
df.transposed2 <- df.transposed[171:357,c("AM3_chronic_3","AM3_acute_1","AM1a_4","AM1a_6","IgM_5","RM_2","AM2_0")]
pdf(file = here::here("analysis","plots","06_DE","RNA_pheatmap_ann_small.pdf"),height = 23)
pheatmap(df.transposed,scale = "row",annotation_row = gene_ann,annotation_col = col_ann, annotation_colors = annoCol, cluster_rows = FALSE) 
dev.off()

pdf(file = here::here("analysis","plots","06_DE","RNA_pheatmap_ann1_small.pdf"),height = 23)
pheatmap(df.transposed1,scale = "row",annotation_row = gene_ann,annotation_col = col_ann, annotation_colors = annoCol, cluster_rows = FALSE,cluster_cols=FALSE) 
dev.off()
pdf(file = here::here("analysis","plots","06_DE","RNA_pheatmap_ann2_small.pdf"),height = 23)
pheatmap(df.transposed2,scale = "row",annotation_row = gene_ann,annotation_col = col_ann, annotation_colors = annoCol, cluster_rows = FALSE,cluster_cols=FALSE) 
dev.off()
#### Protein ####
write.csv(table(enriched_prot$gene, enriched_prot$cluster), file = here::here("analysis","data_objects","06_DE","prot_ann.csv"))

prot_ann <- read.csv(file = here::here("analysis","data_objects","06_DE","prot_ann.csv"),row.names = 1)
prot_ann <- prot_ann[,c("ann1","ann2","ann3")]

prot <- unique(enriched_prot$gene)

dat <- t(as.matrix(seurat@assays$Prot@scale.data[prot,]))
cell_ann <- seurat@meta.data %>% select(Cluster.res0.4)
dat <- cbind(dat,cell_ann$Cluster.res0.4)
colnames(dat)[26] <- "Cluster"
dat <- as.data.frame(dat)
dat <- dat %>% mutate_at(c(1:25),as.numeric)
dat <- dat %>% arrange(Cluster)

df <- dat %>% group_by(Cluster) %>% summarise_at(prot,mean)

df.transposed <- df %>% pivot_longer(cols= -1) %>% pivot_wider(names_from = "Cluster",values_from = "value") %>% rename(gene = name)
df.transposed <- as.data.frame(df.transposed)
rownames(df.transposed) <- df.transposed$gene
df.transposed <- df.transposed[,-1]

my_color_palette <- hue_pal()(7)
annoCol<-list(ann1=c(AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"),
ann2 = c(`NA` = "#FFFFFF" , AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"),
ann3 = c(`NA` = "#FFFFFF" , AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"),
unique.dat.Cluster. = c(AM1a_4 = "#F8766D",AM1a_6 = "#C49A00",AM2_0 = "#53B400",AM3_acute_1 = "#00C094",AM3_chronic_3 = "#00B6EB",IgM_5 = "#A58AFF",RM_2 = "#FB61D7"))

# [1] "AM1a_4"        "AM1a_6"        "AM2_0"         "AM3_acute_1"  
# [5] "AM3_chronic_3" "IgM_5"         "RM_2"  
# [1] "#F8766D" "#C49A00" "#53B400" "#00C094" "#00B6EB" "#A58AFF" "#FB61D7"
col_ann <- data.frame(unique(dat$Cluster))
rownames(col_ann) <- unique(dat$Cluster) 

#re-order columns
df.transposed <- df.transposed[c("P-IgD","P-IgM","P-CD31","P-CD11c","P-CD72","P-FcRL4","P-CD95","P-CD120B","P-CCR6",
"P-ITGB7","P-CD21","P-CD35","P-CD69","P-CD25","P-CD24","P-CXCR5","P-CD45RB","P-CD1c","P-CD71","P-CD86","P-CD38","P-CD62L","P-CD27","P-FCRL5","P-CD10"),
c("RM_2","AM1a_4","AM3_chronic_3", "AM3_acute_1", "AM2_0","IgM_5","AM1a_6")]

pheatmap(df.transposed,scale = "row",annotation_row = prot_ann,annotation_col = col_ann, annotation_colors = annoCol, cluster_rows = FALSE) 

#### now reduce proteins down to list from sarah
prot_list <- read.csv(file = here::here("analysis","data_objects","06_DE","prot_list.csv"))
prot_list$mod <- paste0("P-",prot_list$protein)

prot_ann$mod <- rownames(prot_ann)
prot_ann <- inner_join(prot_list,prot_ann, by = "mod")
rownames(prot_ann) <- prot_ann$mod
prot_ann <- prot_ann[,c("ann1","ann2","ann3")]


df.transposed.small <- df.transposed[rownames(prot_ann),]
df.transposed.small <- df.transposed.small[c("P-CD11c","P-CD72","P-CD95","P-CCR6","P-CD27","P-CD38","P-CD62L",
"P-CD71","P-FCRL5","P-CD21","P-CD35","P-CXCR5","P-CD24","P-CD45RB"),]
pheatmap(df.transposed.small,scale = "row",annotation_row = prot_ann,annotation_col = col_ann, annotation_colors = annoCol, cluster_rows = FALSE) 

#### now use all proteins on list from Sarah
prot_list <- read.csv(file = here::here("analysis","data_objects","06_DE","prot_list.csv"))
prot_list$mod <- paste0("P-",prot_list$protein)

# "P-FCRL4" "P-IGM" "P-IGD"  "P-IGG" "P-CD1C"
# "P-FcRL4" "P-IgM" "P-IgD"  "P-IgG"  "P-CD1c" 

### set up expression 
dat <- t(as.matrix(seurat@assays$Prot@scale.data[prot_list$mod,]))
cell_ann <- seurat@meta.data %>% select(Cluster.res0.4)
dat <- cbind(dat,cell_ann$Cluster.res0.4)
colnames(dat)[25] <- "Cluster"
dat <- as.data.frame(dat)
dat <- dat %>% mutate_at(c(1:24),as.numeric)
dat <- dat %>% arrange(Cluster)

df <- dat %>% group_by(Cluster) %>% summarise_at(prot_list$mod,mean)
df.transposed <- df %>% pivot_longer(cols= -1) %>% pivot_wider(names_from = "Cluster",values_from = "value") %>% rename(gene = name)
df.transposed <- as.data.frame(df.transposed)
rownames(df.transposed) <- df.transposed$gene
df.transposed <- df.transposed[,-1]

col_ann <- data.frame(unique(dat$Cluster))
rownames(col_ann) <- unique(dat$Cluster) 

pheatmap(df.transposed,scale = "row",annotation_col = col_ann, annotation_colors = annoCol, cluster_rows = FALSE)

### annotated dotplot for RNA
more_enriched_genes <- enriched_genes %>% filter(avg_log2FC > 0.8)
pdf(file = here::here("analysis","plots","06_DE","rna_dotplot.pdf"),width = 16)
DotPlot(seurat, assay = "RNA", features = unique(more_enriched_genes$gene)) + theme(axis.text.x = element_text(angle = 90))
dev.off()


#### Top fc violing plots for RNA
top.fc <- enriched_genes %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 2)
pdf(file = here::here("analysis","plots","06_DE","topfc_vlnplot.pdf"),width = 16)
VlnPlot(seurat, assay = "RNA", features = top.fc$gene, ncol = 5, pt.size = 0)
dev.off()


### pairwise analysis between acute and chronic
acute.chronic <- FindMarkers(seurat, assay = "RNA", ident.1 = "AM3_chronic_3", ident.2 = "AM3_acute_1",logfc.threshold = 0.5)
saveRDS(acute.chronic, file = here::here("analysis","data_objects","06_DE","acute_chronic_DE.rds"))
enriched.ac <- acute.chronic %>% filter(avg_log2FC > 0)

DotPlot(seurat %>% filter(Cluster.res0.4 %in% c("AM3_chronic_3","AM3_acute_1")), assay = "RNA", features = rownames(enriched.ac)) + theme(axis.text.x = element_text(angle = 90))
VlnPlot(seurat %>% filter(Cluster.res0.4 %in% c("AM3_chronic_3","AM3_acute_1")), assay = "RNA", features = rownames(enriched.ac), ncol = 7, pt.size = 0)

### pairwise analysis between AM1's 
am1.pairwise <- FindMarkers(seurat, assay = "RNA", ident.1 = "AM1a_4", ident.2 = "AM1a_6",logfc.threshold = 0.5)
saveRDS(am1.pairwise, file = here::here("analysis","data_objects","06_DE","am1_pairwise_DE.rds"))
enriched.am1 <- am1.pairwise %>% filter(avg_log2FC > 0)
DotPlot(seurat %>% filter(Cluster.res0.4 %in% c("AM1a_4","AM1a_6")), assay = "RNA", features = rownames(am1.pairwise)) + theme(axis.text.x = element_text(angle = 90))

### top fc Prot
top.fc.prot <- enriched_prot %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 2)
VlnPlot(seurat, assay = "Prot", features = top.fc.prot$gene, ncol = 5, pt.size = 0)

##### Do DE and return all genes for fgsea #####
seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering", "H2_Mem_G_clusters.rds"))
seurat <- SetIdent(seurat, value = "Cluster.res0.4")

DefaultAssay(seurat) <- "RNA"
seurat <- ScaleData(seurat, features <- rownames(seurat))

de.genes <- FindAllMarkers(object = seurat, min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)