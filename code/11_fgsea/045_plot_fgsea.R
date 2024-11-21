library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(sessioninfo)
library(pheatmap)
#library(stingr)
#library(tidyseurat)
#library(Seurat)

gs_patterns <- c("h.all.v2023","c2.cp.reactome","c7.vax.v2023","c7.immunesigdb.v2023","c5.go.mf.v2023","c5.go.cc.v2023","c5.go.bp.v2023")#
gs_names <- c("Hallmark","Reactome","Vaccine Response","ImmuneSigDB","GO Molecular Function","GO Cellular Component","GO Biological Process") #
hist_plots <- list()
num.sig.paths <- data.frame(gene_set=character(),
                 num=integer(),
                 num.up=integer(),
                 num.down=integer(),
                 stringsAsFactors=FALSE)

for(i in 1:length(gs_patterns)){

    file <- list.files(path = here::here("analysis","data_objects","11_fgsea","new_de_collapse"),pattern = gs_patterns[i],full.names = TRUE)
    file <- file[grep("GSEAresCollapsed.rds", file)]

    #x <- readRDS(file = "/hpcdata/vrc_vip/Abby/Experiment_316/analysis/data_objects/11_fgsea/new_de_collapse/AM2_0_h.all.v2023.2.Hs.symbols.gmt_GSEAresCollapsed.rds" )


    df <- file %>%
    map(readRDS) %>% 
    bind_rows(.id = "cluster") %>% 
    mutate(cluster = recode(cluster,`1` = "AM1a_4", `2` = "AM1a_6", `3` = "AM2_0",`4` = "AM3_acute_1",`5` = "AM3_chronic_3",`6` = "IgM_5", `7` = "RM_2")) 

    sig.path <- df %>% filter(padj < 0.05) %>% select(pathway) %>% distinct()
    df.sig <- df %>% filter(pathway %in% sig.path$pathway)

    ### remove pathways that can be collapsed (aka have a parent pathway)
    df.sig <- df.sig %>% filter(is.na(parentPathway))
    write.csv(table(df.sig$parentPathway,df.sig$cluster), file = here::here("analysis","data_objects","11_fgsea","new_de_collapse",paste0(gs_patterns[i],"_table_parentPathway.csv")))

    save(sig.path,df.sig, file = here::here("analysis","data_objects","11_fgsea","new_de_collapse",paste0(gs_names[i],"_sig_pathways_all_clusters.Rdata")))

    num.sig.paths[i,1] <- gs_names[i]
    num.sig.paths[i,2] <- nrow(sig.path)
    num.sig.paths[i,3] <- sum(df.sig$NES > 0)
    num.sig.paths[i,4] <- sum(df.sig$NES < 0)

    hist_plots[[i]] <- ggplot(df.sig, aes(x=NES)) + geom_histogram() + ggtitle(gs_names[i]) + theme_bw()

    pdf(file = here::here("analysis","plots","11_fgsea","new_de_collapse",paste0(gs_patterns[i],"dotplot.pdf")),width = 12,height = 35)
    ggplot(df.sig, aes(x= cluster, y=pathway, size=NES, color=padj)) + geom_point(alpha = 0.8) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    dev.off()


    ##prepare data for pheatmap
    mat <- df.sig %>% select(cluster,NES,pathway) %>% pivot_wider(names_from = cluster, values_from = NES) #%>% na.omit()
    mat <- as.data.frame(mat)
    rownames(mat) <- mat$pathway
    mat <- mat[,-1]
    mat[is.na(mat)] <- 0

    pdf(file = here::here("analysis","plots","11_fgsea","new_de_collapse",paste0(gs_patterns[i],"heatmap.pdf")), height = 30, width = 16)
    pheatmap(mat)
    dev.off()

}

    pdf(file = here::here("analysis","plots","11_fgsea","new_de_collapse","NES_histogram.pdf"), height = 12, width = 16)
    patchwork::wrap_plots(hist_plots,ncol = 2)
    dev.off()

    write.csv(num.sig.paths, file = here::here("analysis","data_objects","11_fgsea","new_de_collapse","num_sig_paths.csv"))

    x <- num.sig.paths %>%
        select(gene_set,num.up,num.down) %>%
        pivot_longer(!gene_set, names_to = "status", values_to = "count") %>% as.data.frame()

    pdf(file = here::here("analysis","plots","11_fgsea","new_de_collapse","num_sig_pathways.pdf"), height = 12, width = 16)
        ggplot(data=x, aes(x=gene_set, y=count, fill=status)) +
            geom_bar(stat="identity", position=position_dodge())+
            scale_fill_brewer(palette="Paired")+
            theme_minimal() +
            theme(axis.text.y=element_text(size=16),axis.text.x=element_text(size=12),
                axis.title=element_text(size=14,face="bold"))
    dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
