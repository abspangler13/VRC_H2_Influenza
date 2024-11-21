### make some plots for manuscript
### protein expression of CD21, CD27, CD95, CXCR3, CD38, CD71, Cd1c, CD62L, CD11c, FCRL5, CD22, CD72, CD32, CD85j, CD45RB, CCR6, FCRL4
### vertical for figure 5C. genes as rows, clusters as columns
### cluster order and naming: 
# C1 – Naïve
# C2 – RM
# C3 – AM2
# C4 – AM1a_4
# C5 – AM1b_6
# C6 – AM3_acute
# C7 – AM3_chronic

library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

#https://divingintogeneticsandgenomics.com/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-2, 0, -2, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(gsub("P-","",feature)) + ggtitle("") + 
    labs(title = NULL) +
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = 24, angle = 90), 
          axis.text.y = element_text(size = 20), 
          plot.margin = plot.margin) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-2, 0, -2, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust = 1,size=18), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))
features<- c("P-CD21", "P-CD27", "P-CD45RB", "P-CD38",  "P-CD62L", "P-CD71",  "P-FCRL5", "P-CD11c","P-CD22","P-CD72","P-CD32","P-CD95", "P-CD1c", "P-CXCR3", "P-CCR6", "P-FcRL4","P-CD73")

custom_palette = c("#FF8C00","#FFD700","#9F79EE","#27408B","#FF4040","#8B2323","#BEBEBE","#000000")
names(custom_palette) <- c("AM1a_4","AM1b_6","AM2","RM","AM3_acute","AM3_chronic","Naïve","Germline")

seurat@meta.data <-seurat@meta.data %>% mutate(Annotation = case_when(Cluster.res0.4 == "RM_2" ~ "RM",
                                                        Cluster.res0.4 == "AM1a_4" ~ "AM1a_4",
                                                        Cluster.res0.4 == "AM3_chronic_3" ~ "AM3_chronic",
                                                        Cluster.res0.4 == "AM3_acute_1" ~ "AM3_acute",
                                                        Cluster.res0.4 == "AM2_0" ~ "AM2",
                                                        Cluster.res0.4 == "IgM_5" ~ "Naïve",
                                                        Cluster.res0.4 == "AM1a_6" ~ "AM1b_6"))

seurat$Annotation <- factor(x = seurat$Annotation, levels = c('Naïve', 'RM','AM2','AM1b_6','AM1a_4','AM3_acute','AM3_chronic'))
Idents(seurat) <- "Annotation"

pdf(file = here::here("analysis","plots","05_clustering","prot_violin_5C.pdf"),height = 60)
VlnPlot(seurat, features = features, pt.size=0, cols=custom_palette, flip = TRUE, stack = TRUE,fill.by = "ident") #
#StackedVlnPlot(obj = seurat, features = features)
dev.off()

## need to re-order clusters. convert ann to a factor and specify levels. export as vector file and reduce white space in between plots 
## make y axis scale label smaller. remove "P-" from y-axis label
## plot.margin order = top, right, bottom, left

pdf(file = here::here("analysis","plots","05_clustering","stacked_prot_violin_5C.pdf"),height = 25)
StackedVlnPlot(obj = seurat, features = features,ncol = 1, pt.size=0, cols=custom_palette)
dev.off()

svg(file = here::here("analysis","plots","05_clustering","stacked_prot_violin_5C.svg"),height = 25)
StackedVlnPlot(obj = seurat, features = features,ncol = 1, pt.size=0, cols=custom_palette)
dev.off()

