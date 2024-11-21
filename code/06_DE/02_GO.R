# enrichGO https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/11_gene_ontology/deprecated/gene_ontology.R

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(readr)
library(ggplot2)
library(sessioninfo)

tonsil_in <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))

DefaultAssay(tonsil_in) <- "RNA"
tonsil_in <- SetIdent(tonsil_in, value = "Cluster.res0.4")

all.markers <- readRDS(file = here::here("analysis","data_objects","06_DE","rna_markers_0.4_seurat.rds"))
enriched.markers <- all.markers %>% filter(avg_log2FC > 0.0,p_val_adj < 0.05)
un.enriched.markers <- all.markers %>% filter(avg_log2FC < -0.0, p_val_adj < 0.05)

saveRDS(enriched.markers,file = here::here("analysis","data_objects","06_DE","rna_enriched_markers.rds"))
saveRDS(un.enriched.markers,file = here::here("analysis","data_objects","06_DE","rna_unenriched_markers.rds"))

# make gene lists for enrichR online tool
clusters <- unique(enriched.markers$cluster)
gene.lists <- list()
deg.stat <- data.frame(tot.deg=rep(NA,12),
                        mito=rep(NA,12),
                        ig=rep(NA,12),
                        ribo=rep(NA,12))
for(i in 1:length(clusters)){
    gene.list <- enriched.markers %>% filter(cluster == clusters[i])
    deg.stat$tot.deg[i] <- dim(gene.list)[1]
    #remove mito genes
    gene.list <- gene.list %>% filter(!grepl("^MT-",gene))
    deg.stat$mito[i] <- deg.stat$tot.deg[i] - dim(gene.list)[1]
    #remove variable IG genes
    gene.list <- gene.list %>% filter(!grepl("IG[HKL]V|TRBV",gene))
    deg.stat$ig[i] <- deg.stat$tot.deg[i] - dim(gene.list)[1]
    #remove ribosomal genes
    gene.list <- gene.list %>% filter(!grepl("^RP[LS]",gene))
    deg.stat$ribo[i] <- deg.stat$tot.deg[i] - dim(gene.list)[1]

    write.csv(gene.list, file = here::here("analysis","data_objects","06_DE",paste0(clusters[i],"_enriched_markers.csv")))
    write.csv(deg.stat, file = here::here("analysis","data_objects","06_DE","deg_stats.csv"))
    x <- as.vector(gene.list %>% dplyr::select(gene))
    gene.lists[[i]] <- bitr(x$gene,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}


names(gene.lists) <- clusters

clust_CC <- list()
clust_MF <- list()
clust_BP <- list()

for(i in 1:length(clusters)){
    clust_CC[[i]] <- enrichGO(
    gene = gene.lists[[i]]$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

clust_MF[[i]] <- enrichGO(
    gene = gene.lists[[i]]$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

clust_BP[[i]] <- enrichGO(
    gene = gene.lists[[i]]$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

}

clust_compare <- list()
for(i in 1:length(clusters)){
    clust_compare[[i]] <- gene.lists[[i]]$ENTREZID
}
names(clust_compare) <- clusters

comp_CC <- compareCluster(clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "CC",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

comp_MF <- compareCluster(clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "MF",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

comp_BP <- compareCluster(clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "BP",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

save(gene.lists,clust_CC,clust_BP,clust_MF,comp_BP,comp_CC,comp_MF, file = here::here("analysis","data_objects","06_DE","enrichGO_results.rds"))
# load(file = here::here("analysis","data_objects","06_DE","enrichGO_results.rds"))

## comparing clusters
pdf(file = here::here("analysis","plots","06_DE","all_comparative_GO.pdf"), width = 18, height = 14)

dotplot(comp_CC, showCategory = 9, label_format = 50) +
    ggtitle("Top 9 Comparative GO Cellular Compartment")
dotplot(comp_MF, showCategory = 9, label_format = 50) +
    ggtitle("Top 9 Comparative GO Molecular Function")
dotplot(comp_BP, showCategory = 9, label_format = 50) +
    ggtitle("Top 9 Comparative GO Biological Process")

comp_CC <- pairwise_termsim(comp_CC)
emapplot(comp_CC,
    showCategory = 9, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 9 Comparative GO Cellular Compartment Modules")
comp_MF <- pairwise_termsim(comp_MF)
emapplot(comp_MF,
    showCategory = 9, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 9 Comparative GO Molecular Function Modules")
comp_BP <- pairwise_termsim(comp_BP)
emapplot(comp_BP,
    showCategory = 9, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 9 Comparative GO Biological Process Modules")

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()