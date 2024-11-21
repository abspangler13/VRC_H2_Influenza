# https://htmlpreview.github.io/?https://github.com/leeprichman/ClonoCluster/blob/main/Tutorial.html

library(data.table)
library(magrittr)
library(Seurat)
library(ClonoCluster)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(tidyseurat)

seurat <- readRDS("230927_files_SFA/H2.M.G")

##randomly downsample seurat just because my computer can't handle it
seurat <- seurat %>% subset(downsample = 500)

# read in barcodes (clone id's for us)
my.bt <- as.data.table(seurat %>% select(CELL,clone_subject_id))
colnames(my.bt) <- c("rn","Barcode")

# what does it look like
head(my.bt)

# read in count matrix
my.cm <-as.data.frame(seurat@assays$RNA@data)
my.cm$rn <- rownames(my.cm)
my.cm <- my.cm[,c(3893,1:3892)]
my.cm <- as.data.table(my.cm)
my.cm <- tdt(my.cm)
my.cm <- dt2m(my.cm)

# what does it look like
my.cm[1:5,1:5]

# pca <- seurat@reductions$harmony
# need to figure out why rownames are genes for seurat PCA instead of cell barcodes

my.pca <- irlba_wrap(my.cm, npc = 25)
my.pca[1:5, 1:5]

# lets get our range of alpha values

my.als <- seq(0, 1, by = 0.1)

# return the cluster assignments for range of alphas
clust <- clonocluster(my.pca, my.bt, alpha = my.als, beta = 0.1, res = 1.5)

saveRDS(clust, file = "ClonoCluster_clust.rds")
head(clust)

bt <- my.bt
cm <- my.cm
pca <- my.pca
als <- my.als

# These are your hybrid clusters. Let’s compute some confusion matrix statistics and plot them just so we can see the trends. 
# Removing singlets makes the visualization much better so we will do that as well.

# loop over alphas
confusion <- lapply(clust[, alpha %>% unique], function(a){
  
  # function to compute confusion statistics
  ct <- cast_confusion(clust[alpha == a, .(rn, Group)], # subsetted, 2 column table is the input
                       bt) # barcode assignments
  
  # append alpha value                      
  ct[, alpha := a]
  
  return(ct)
  
}) %>% data.table::rbindlist() # merge to one big table

# get barcode ids of singlets
singlets <- bt[, .N, by = Barcode] %>% .[N == 1, Barcode]

# plot the cohen's kappa
ggplot(confusion[!barcode %chin% singlets], # ignore singlets
       aes(x = as.factor(alpha), y = cohens_k)) +
  geom_boxplot(fill = "dodgerblue") +
  theme_bw() +
  ttheme +
  theme(plot.title = element_text(size = 12)) +
  ggtitle("Concordance plot") +
  xlab("\u03B1") +
  ylab("Cohen's \u03BA")

# returns a table with column V1L number of groups at that value of alpha
nt <- clust[, Group %>% unique %>% length, by = "alpha"]

# plot it
ggplot(nt, aes(x = alpha, y = V1)) +
  geom_point(color = "dodgerblue", size = 3) +
  geom_line(color = "grey", linetype = "dashed") +
  theme_bw() +
  ttheme +
  xlab("\u03B1") +
  ylab("# of clusters") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))

notable_alphas <- c(0, 0.2, 0.4, 0.9)

# function to plot two sankeys colored by both ending nodes
p <- Plot_alluvia(clust[alpha %in% notable_alphas], # subset table on these alpha values
                  bt, # barcode table
                  title = "Sample Sankey",
                  xlab = "\u03B1", # unicode symbol for alpha
                  ylab = "# of cells",
                  border_size = 1, # border around the nodes
                  label_nodes = FALSE, # labels nodes but hard to viz here
                  cols = ClonoCluster::c25 # colors of nodes and ribbons
)

# change the subtitle of the second plot because we aren't coloring by barcodes
p[[2]] <- p[[2]] + ggtitle("Sample Sankey", subtitle = "Colored by \u03B1 = 0.9 clusters")

# plot colored by initial transcriptome
p[[1]]

# plot colored by highest alpha value in table
p[[2]]

# to skip this step:
# auc_table <- data.table::fread(file.path(dir, "YG1_markers.txt"))
auc_all <- Find_Markers_ROC(clust[alpha %in% notable_alphas], cm, n_threads = 1)

# you can save this table with:
# auc_table %>% data.table::fwrite("mytable.txt")

auc_all <- FindAllMarkers_Seurat(so, # seurat object
                                 clust, # ClonoCluster function output
                                 test.use = "roc") # if you let this default you will get wilcox.tests

# what does it look like?
head(auc_all)

# If you want to take the best auc for any cluster at each alpha
auc_table <- auc_all[order(-auc)] %>% unique(by = c("rn", "alpha"))

head(auc_table)