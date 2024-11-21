# module load immcantation/4.4.0
# singularity shell $IMMCAN_SING
# R

### Make plots of trees we just built

# Load required packages
# library(alakazam)
library(dowser)
library(dplyr)
library(ggtree)
# library(readr)
# library(ggplot2)
# library(RColorBrewer)


trees <- readRDS(file = "analysis/data_objects/10_dowser/HL_trees.rds")
plots <- plotTrees(trees, tipsize=3)
treesToPDF(plots, here::here("analysis","plots","10_dowser","HL_tree_plots_raxml.pdf"), nrow = 1, ncol = 2)


###### 1 ##### Cluster.res0.4

# AM1a_4: “darkorange” 	"#FF8C00"
# AM1a_6: “gold1” "#FFD700"
# AM2: “mediumpurple2” "#9F79EE"
# RM: “royalblue4” "#27408B"
# AM3_acute: “brown1” "#FF4040"
# AM3_chronic: “brown4” "#8B2323"
# IgM: “gray” "#BEBEBE"


custom_palette = c("#FF8C00","#FFD700","#9F79EE","#27408B","#FF4040","#8B2323","#BEBEBE","#000000")
names(custom_palette) <- c("AM1a_4","AM1a_6","AM2_0","RM_2","AM3_acute_1","AM3_chronic_3","IgM_5","Germline")

trees <- readRDS(file = "analysis/data_objects/10_dowser/HL_trees.rds")
plots <- plotTrees(trees, tips="Cluster.res0.4",tip_palette = custom_palette,tipsize=3)
titles <- trees %>% select(clone_id,vac_grp,Specificity) %>% distinct(clone_id,.keep_all = TRUE)
for(i in 1:nrow(titles)){
        plots[[i]] <- plots[[i]] + ggtitle(paste0(titles$clone_id[i]," ",titles$vac_grp[i]," ",titles$Specificity[i]))
}
treesToPDF(plots, here::here("analysis","plots","10_dowser","HL_tree_plots_raxml_cluster.pdf"), nrow = 1, ncol = 2)


##### 2 ##### Timepoint.num
# trees <- readRDS(file = "analysis/data_objects/10_dowser/Trees_igphyml_Time_meta.rds")
plots <- plotTrees(trees, tips="Timepoint.num",tip_palette = "RdYlBu", tipsize = 3)
plots <- lapply(plots, function(x){
    x + geom_tiplab(aes(label = Timepoint.wk),offset = 0.03)
})

titles <- trees %>% select(clone_id,vac_grp,Specificity) %>% distinct(clone_id,.keep_all = TRUE)
for(i in 1:nrow(titles)){
        plots[[i]] <- plots[[i]] + ggtitle(paste0(titles$clone_id[i]," ",titles$vac_grp[i]," ",titles$Specificity[i]))
}
treesToPDF(plots, here::here("analysis","plots","10_dowser","tree_plots_phyigml_timepoint.pdf"), nrow = 1, ncol = 2, height = 20)


##### 3 ##### c_call
trees <- readRDS(file = "analysis/data_objects/10_dowser/Trees_igphyml_c_call_meta.rds")
plots <- plotTrees(trees, tips="c_call",tipsize=3,tip_palette = "Set1")

titles <- trees %>% select(clone_id,vac_grp,Specificity) %>% distinct(clone_id,.keep_all = TRUE)
for(i in 1:nrow(titles)){
        plots[[i]] <- plots[[i]] + ggtitle(paste0(titles$clone_id[i]," ",titles$vac_grp[i]," ",titles$Specificity[i]))
}
treesToPDF(plots, here::here("analysis","plots","10_dowser","tree_plots_phyigml_c_call.pdf"), nrow = 1, ncol = 2)


###### 4 ##### "Timepoint.num","Cluster.res0.4"
# plots <- plotTrees(trees)
# plots <- lapply(plots, function(x){
#     x + geom_tippoint(mapping = aes(shape = Cluster.res0.4, color = Timepoint.num,size = 4)) +
#       scale_color_gradient2(low="blue",mid = "yellow",high = "red", midpoint = 7) 
# })
# titles <- trees %>% select(clone_id,vac_grp,Specificity) %>% distinct(clone_id,.keep_all = TRUE)
# for(i in 1:nrow(titles)){
#         plots[[i]] <- plots[[i]] + ggtitle(paste0(titles$clone_id[i]," ",titles$vac_grp[i]," ",titles$Specificity[i]))
# }
# treesToPDF(plots, here::here("analysis","plots","10_dowser",paste0("tree_plots_phyigml_timepoint_cluster.pdf")), nrow = 1, ncol = 2)

##### 5 ##### "Timepoint.num",'c_call'

# plots <- plotTrees(trees)
# plots <- lapply(plots, function(x){
#     x + geom_tippoint(mapping = aes(shape = c_call, color = Timepoint.num,size = 4)) +
#       scale_color_gradient2(low="blue",mid = "yellow",high = "red", midpoint = 7) 
# })
# titles <- trees %>% select(clone_id,vac_grp,Specificity) %>% distinct(clone_id,.keep_all = TRUE)
# for(i in 1:nrow(titles)){
#         plots[[i]] <- plots[[i]] + ggtitle(paste0(titles$clone_id[i]," ",titles$vac_grp[i]," ",titles$Specificity[i]))
# }
# treesToPDF(plots, here::here("analysis","plots","10_dowser",paste0("tree_plots_phyigml_",names[k],".pdf")), nrow = 1, ncol = 2)

##### 6 ##### "Timepoint.num",'c_call', cluster

# plots <- plotTrees(trees)
# plots <- lapply(plots, function(x){
#     x + geom_tippoint(mapping = aes(shape = c_call, color = Timepoint.num,size = 4)) +
#       scale_color_gradient2(low="blue",mid = "yellow",high = "red", midpoint = 7) 
# })
# titles <- trees %>% select(clone_id,vac_grp,Specificity) %>% distinct(clone_id,.keep_all = TRUE)
# for(i in 1:nrow(titles)){
#         plots[[i]] <- plots[[i]] + ggtitle(paste0(titles$clone_id[i]," ",titles$vac_grp[i]," ",titles$Specificity[i]))
# }
# treesToPDF(plots, here::here("analysis","plots","10_dowser",paste0("tree_plots_phyigml_",names[k],".pdf")), nrow = 1, ncol = 2)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()