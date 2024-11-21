# Load required packages
library(alakazam)
library(dowser)
library(dplyr)
library(ggtree)

# load example AIRR tsv data
my.data <- readRDS(file = here::here("analysis","data_objects","10_dowser","seurat_germline.rds"))
my.data <- my.data@meta.data

#"B007" "B014" "B028" "B090" "B000" "B180" "A000" "A028" "A007"
my.data <- my.data %>% mutate(Timepoint.num = case_when(Timepoint == "A000" ~ 1,
                                                        Timepoint == "A007" ~ 2,
                                                        Timepoint == "A028" ~ 3,
                                                        Timepoint == "B000" ~ 4,
                                                        Timepoint == "B007" ~ 5,
                                                        Timepoint == "B014" ~ 6,
                                                        Timepoint == "B028" ~ 7,
                                                        Timepoint == "B090" ~ 8,
                                                        Timepoint == "B180" ~ 9))


######
B180 <- my.data %>% filter(Timepoint == "B180")
B180 <- unique(B180$clone_subject_id)

## subset data to only include these clones
my.data <- my.data %>% filter(clone_subject_id %in% B180)
######

x <- as.data.frame(table(my.data$clone_subject_id))
x <- x[order(x$Freq),]

x.clones <- x %>% filter(Freq > 50) %>% select(Var1)

# subset data for this example
#ExampleAirr = ExampleAirr[ExampleAirr$clone_id %in% c("3170", "3184"),]
my.data <- my.data[my.data$clone_subject_id %in% x.clones$Var1,]

# Process example data into proper format, store isotype (optional)
my.clones = formatClones(my.data, clone = "clone_subject_id",traits = c("Cluster.res0.4","Timepoint.num","c_call")) # "Timepoint" 

trees = getTrees(my.clones, build="pml", trait="Timepoint.num", igphyml="/usr/local/share/igphyml/src/igphyml")

saveRDS(trees,file = here::here("analysis","data_objects","10_dowser","B180_timepoint_cluster_ccall_trees_pml.rds"))

#trees <- scaleBranches(trees,edge_type = "genetic_distance")
#trees <- scaleBranches(trees)

# simple tree plotting with ggtree R package with isotypes at tips
plots <- plotTrees(trees, tips="Timepoint.num",tipsize=2,tip_palette = "RdYlBu") #
plots_2 <- plotTrees(trees, tips="Cluster.res0.4",tipsize=2,tip_palette = "Paired")

#check plots
#plots[[1]]

# plots <- lapply(plots, function(x){
#     x + geom_tippoint(aes(colour=Cluster.res0.4,size = Timepoint.num)) 
# })

# make pdf of all trees
treesToPDF(plots, here::here("analysis","plots","10_dowser","B180_timepoint_cluster_ccall_phyigml.pdf"), nrow = 1, ncol = 2)
treesToPDF(plots_2, here::here("analysis","plots","10_dowser","B180_timepoint_cluster_ccall_phyigml.pdf"), nrow = 1, ncol = 2)
