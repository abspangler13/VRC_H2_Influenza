library(Seurat) 
library(tidyverse)


####KEGG. Take enriched marker gene csv files for each cluster and paste the gene names into https://maayanlab.cloud/Enrichr/. Go to pathways, KEGG and export the table. 
clusters <- unique(seurat$Cluster.res0.4)

# [1] "RM_2"          "AM1a_4"        "AM3_chronic_3" "AM3_acute_1"  
# [5] "AM2_0"         "IgM_5"         "AM1a_6" 

kegg.list <- list()
for(i in 1:length(clusters)){
    kegg.list[[i]] <- read_tsv(file = here::here("analysis","data_objects","06_DE","KEGG",paste0(clusters[i],"_KEGG_2021_Human_table (14).txt")))
}
names(kegg.list) <- clusters

terms <- c()
for(i in 1:length(kegg.list)){
    x <- kegg.list[[i]]$Term
    terms <- c(terms,x)
}

unique.terms <- unique(terms)

dat <- data.frame()
for(i in 1:length(clusters)){
    x <- kegg.list[[i]] %>% dplyr::select(Term,`Adjusted P-value`,`Combined Score`) %>% mutate(cluster = names(kegg.list)[i])
    dat <- rbind(dat,x)
}

pdf(file = here::here("analysis","plots","06_DE","KEGG_compare_clust_filtered.pdf"),height = 32)
ggplot(dat, aes(x= cluster, y=Term, size=`Combined Score`, color=`Adjusted P-value`)) + geom_point(alpha = 0.8) + 
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()