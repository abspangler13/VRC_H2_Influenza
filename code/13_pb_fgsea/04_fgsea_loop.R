library(fgsea)
library(data.table)
library(tidyverse)
library(tidyseurat)
library(Seurat)
library(sessioninfo)

#code adapted from this tutorial: https://biostatsquid.com/fgsea-tutorial-gsea/ 

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# Analysis ====================================================

## 1. Read in data -----------------------------------------------------------
df_raw <- readRDS(file = here::here("analysis","data_objects","12_pb_clustering","pb_de_rna_all_genes.rds")) #give it all the genes

seurat_raw <- readRDS(file = here::here("analysis","data_objects","12_pb_clustering", "H2P_G.rds"))
seurat_raw <- SetIdent(seurat_raw, value = "PB.Clusters.res01")

clusters <- unique(seurat_raw$PB.Clusters.res01)

## 2. Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
bg_path <- here::here("analysis","data_objects","13_pb_fgsea","fgsea_dbs")  #Path to the downloaded fgsea databases/genesets


for(i in 4:length(clusters)) {
  print(paste0("running cluster ",clusters[i]))

  #subset seurat object and DE object to contain only cells from cluster of interest
  df <- df_raw %>% filter(cluster == clusters[i])
  seurat <- seurat_raw %>% filter(PB.Clusters.res01 == clusters[i])

  # For GSEA
    # Filter out the gmt files for KEGG, Reactome and GOBP

    #Identify the genes I have in my data
    my_genes <- as.vector(rownames(seurat@assays$RNA))

    #Rank the genes by fold change
    rankings <- df$avg_log2FC 
    names(rankings) <- df$gene # genes as names#
    head(rankings)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

    # Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
    max_ranking <- max(rankings[is.finite(rankings)])
    min_ranking <- min(rankings[is.finite(rankings)])
    rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
    rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
    
    saveRDS(rankings, file = here::here("analysis","data_objects","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_gene_rank_stats.rds")))

    pdf(file = here::here("analysis","plots","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_gene_rankings.pdf")))
    plot(rankings)
    ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
        geom_point() +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    dev.off()

    list.files(bg_path,pattern = '.gmt')
    gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
    gmt_files
    for(j in 1:length(gmt_files)){
        print(paste0("running cluster ",clusters[i]," and gene set ",list.files(bg_path,pattern = '.gmt')[j]))

        ##Prepare the gene sets we downloaded by subsetting it with genes from my data
        bg_genes <- prepare_gmt(gmt_files[j], my_genes, savefile = TRUE)  ##saves to same folder where original .gmt file is located

        ## 4. Run GSEA ---------------------------------------------------------------
        # Easy peasy! Run fgsea with the pathways 
        GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                        stats = rankings,
                        scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                        minSize = 10,
                        maxSize = 500,
                        nproc = 1,
                        nPermSimple = 1000000) # for parallelisation
        
        print(head(GSEAres[order(padj), ]))
        print(paste0("significantly up regulated pathways ",sum(GSEAres[, padj < 0.05 & ES > 0]))) ##print the number of significantly upregulated genesets
        print(paste0("significantly down regulated pathways ",sum(GSEAres[, padj < 0.05 & ES < 0]))) ##print the number of significantly downregulated genesets
        sum(is.na(GSEAres$NES)) 

        saveRDS(GSEAres, file = here::here("analysis","data_objects","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_",list.files(bg_path,pattern = '.gmt')[j],'_GSEAres.rds')))
        ## 6. Check results ------------------------------------------------------
        # Top 6 enriched pathways (ordered by p-val)
        #GSEAres %>% filter(padj < 0.05) %>% arrange(NES)

        #how many significant pathways
        if(sum(GSEAres[, padj < 0.05 & ES > 0]) > 0){
            topPathwaysUp <- GSEAres[ES > 0 & padj < 0.05][head(order(padj)), pathway]
            
            # pdf(file = here::here("analysis","plots","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_",list.files(bg_path,pattern = '.gmt')[j],'_gsea_top_pathways_up.pdf')), width = 20, height = 15)
            # plot(plotGseaTable(bg_genes[topPathwaysUp], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5))
            # dev.off()

            # plot the significantly enriched pathway
            pdf(file = here::here("analysis","plots","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_",list.files(bg_path,pattern = '.gmt')[j],'_gsea_top_pathways_up_plot.pdf')))
            for(k in 1:length(topPathwaysUp)){
                myplot <- plotEnrichment(bg_genes[[topPathwaysUp[k]]],
                        rankings) + labs(title = topPathwaysUp[k])
                plot(myplot)
            }
            dev.off()
        } 
        if(sum(GSEAres[, padj < 0.05 & ES < 0]) > 0){
            topPathwaysDown <- GSEAres[ES < 0 & padj < 0.05][head(order(padj)), pathway]

            # pdf(file = here::here("analysis","plots","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_",list.files(bg_path,pattern = '.gmt')[j],'_gsea_top_pathways_down.pdf')), width = 20, height = 15)
            # plot(plotGseaTable(bg_genes[topPathwaysDown], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5))
            # dev.off()

             # plot the significantly enriched pathway
            pdf(file = here::here("analysis","plots","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_",list.files(bg_path,pattern = '.gmt')[j],'_gsea_top_pathways_down_plot.pdf')))
            for(k in 1:length(topPathwaysDown)){
                myplot <- plotEnrichment(bg_genes[[topPathwaysDown[k]]],
                        rankings) + labs(title = topPathwaysDown[k])
                plot(myplot)
            }
            dev.off()
        }    
        
        #Collapse Pathways that can be collapsed 
        collapsedPathways <- collapsePathways(
          fgseaRes = GSEAres[order(pval)][padj < 0.05],
          pathways = bg_genes,
          stats = rankings
          )
          saveRDS(collapsedPathways, file = here::here("analysis","data_objects","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_",list.files(bg_path,pattern = '.gmt')[j],"_CollapsePathways.rds")))
          GSEAresSig <- GSEAres[order(pval)][padj < 0.05]
          GSEAresSig$parentPathway <- collapsedPathways$parentPathways
          saveRDS(GSEAresSig,file = here::here("analysis","data_objects","13_pb_fgsea","new_de_collapse",paste0(clusters[i],"_",list.files(bg_path,pattern = '.gmt')[j],"_GSEAresCollapsed.rds")))
    }

}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


