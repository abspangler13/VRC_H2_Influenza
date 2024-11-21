# cd /data/vrc_bsc/Abby/Experiment_316
# srun --pty bash
# apptainer shell --bind /data/vrc_bsc:/mnt ../../containers/immcantation/immcantation_suite-4.5.0.sif
# R
library(dowser)
library(vroom)
library(dplyr)
library(tidyr)
library(ggtree)

## Build trees based on heavy and light chain: https://dowser.readthedocs.io/en/stable/vignettes/Resolve-Light-Chains-Vignette/

##figure out how to add in meta data clone_id prefixes

#Read in output of light_cluster.py and light chain data.
clones <- vroom(list.files(path = here::here("analysis","VDJ_run3_run4_by_subject"),pattern = "*heavy_light_clone-pass_germ-pass.tsv", full.names = TRUE)) 
dim(clones)
# [1] 61492    41

#check for multiple sequence IDs per cell_id?? Each cell should only have one heavy chain. 
length(unique(clones$sequence_id)) #[1] 61492
length(unique(clones$cell_id)) # [1]  61492
## only one H chain per data

lc_clones <- vroom(list.files(path = here::here("analysis","VDJ_run3_run4_by_subject"),pattern = "*light_parse-select.tsv", full.names = TRUE)) 
dim(lc_clones)
#[1] 78150    37
length(unique(lc_clones$sequence_id)) #[1] 78150
length(unique(lc_clones$cell_id)) #[1] 74162
#have cells with multiple light chains, but this is ok. 


# Now combine and H and L chain data according to advice from Ken
#  I think light_cluster.py only returns H chain data, so you might have to map the clone_id of each cell from the H chain onto their L chains before inputting to resolveLightChains. 
map.clones <- clones %>% select(cell_id,clone_id)
lc_clones <- lc_clones %>% left_join(map.clones,by = "cell_id")

table(is.na(lc_clones$clone_id))

# FALSE  TRUE 
# 64343 13807 

### remove germline columnes from H data
clones <- clones[,1:37]

my.data <- rbind(clones,lc_clones)  ## when I do this, I get duplicate sequence ID's 
dim(my.data)
# [1] 139642     37

#  You should only need to run createGermlines once - after LightChains. The concatenation is handled by formatClones.

# find the clone subgroups 
my.data <- resolveLightChains(my.data)

print(my.data$clone_subgroup)

# run createGermlines -- this will create new germline for each locus in each subgroup 
# the directory for the references matches the location on docker
references <- readIMGT("/usr/local/share/germlines/imgt/human/vdj")
my.data <- createGermlines(my.data, references = references, clone = "clone_subgroup_id", nproc = 1)
saveRDS(my.data, file = here::here("analysis","data_objects","10_dowser","heavy_light_germline.rds"))
#my.data <- readRDS(file = here::here("analysis","data_objects","10_dowser","heavy_light_germline.rds"))

## add in meta data so I can include them as traits?
my.meta <- readRDS(file = here::here("analysis","data_objects","10_dowser","seurat_germline.rds"))
my.meta <- my.meta@meta.data

my.meta <- my.meta %>% mutate(Timepoint.num = case_when(Timepoint == "A000" ~ 0,
                                                        Timepoint == "A007" ~ 7,
                                                        Timepoint == "A028" ~ 28,
                                                        Timepoint == "B000" ~ 140,
                                                        Timepoint == "B007" ~ 147,
                                                        Timepoint == "B014" ~ 154,
                                                        Timepoint == "B028" ~ 168,
                                                        Timepoint == "B090" ~ 202,
                                                        Timepoint == "B180" ~ 292))
my.meta <- my.meta %>% mutate(Timepoint.wk = case_when(Timepoint == "A000" ~ "0",
                                                        Timepoint == "A007" ~ "1",
                                                        Timepoint == "A028" ~ "4",
                                                        Timepoint == "B000" ~ "16",
                                                        Timepoint == "B007" ~ "17",
                                                        Timepoint == "B014" ~ "18",
                                                        Timepoint == "B028" ~ "20",
                                                        Timepoint == "B090" ~ "28",
                                                        Timepoint == "B180" ~ "40"))

my.meta <- my.meta %>% mutate(exposure = case_when(vac_grp == "3A" ~ "Naive",
                                                        vac_grp == "3B" ~ "Exposed",
                                                        vac_grp == "4A" ~ "Naive",
                                                        vac_grp == "4B" ~ "Exposed"))
## select just meta data columns we need. Join on cell_id
my.meta <- my.meta %>% select(orig.ident, vac_grp, Subject, Timepoint, CELL, specificity.SFA, Spec.Broad, Cluster.res0.4,Timepoint.num, Timepoint.wk, exposure) 
colnames(my.meta)[5] <- "cell_id"

##figure out how to left join by cell id.  find overlap of cell ids from data and meta data. Remove 81584 sequences without meta data?
my.data <- left_join(my.data,my.meta,by = "cell_id")

my.data <- my.data %>% drop_na(vac_grp)

#concatenate subject to clone_id
my.data$clone_id <- paste0(my.data$Subject, "_", my.data$clone_id)

# filter down to only clones Sarah is interested in
B180 <- read.csv(file = here::here("analysis","data_objects","10_dowser","clones.for.dowser.csv")) #25 cells, cell at last timepoint

# length(unique(my.data$clone_id)) my.data has 9396 clones
my.data <- my.data %>% filter(clone_id %in% B180$clone_id)

## 108 after filtering

## do I need to include other traits here for down stream analysis?
my.clones <- formatClones(my.data, traits = c("Timepoint.num"), chain="HL", text_fields = c("vac_grp", "Subject","Spec.Broad","Timepoint","Timepoint.wk","exposure"),
                       nproc=1, collapse = FALSE, 
                       split_light = TRUE, minseq = 3,filterstop=FALSE ) #can set split light = FALSE

saveRDS(my.clones, file = here::here("analysis","data_objects","10_dowser","HL_formatClones_Timepoint.rds"))
print(my.clones)

#my.clones <- readRDS(file = here::here("analysis","data_objects","10_dowser","HL_formatClones_Timepoint.rds"))

# Building maximum likelihood trees with multiple partitions using IgPhyML 
# Only the newest version of IgPhyML supports this option
# exec here is set to IgPhyML position in the Docker image.
#my.clones <- readRDS(file = here::here("analysis","data_objects","10_dowser","HL_formatClones.rds"))

# my.trees <- getTrees(my.clones[1947,], build="igphyml", nproc=1, partition="hl",
#                    exec="/usr/local/share/igphyml/src/igphyml")

my.trees = getTrees(my.clones, build="raxml", 
    exec="/usr/local/bin/raxml-ng", nproc=1, partition="scaled")

saveRDS(my.trees, file = here::here("analysis","data_objects","10_dowser","HL_trees_Timepoint_raxml.rds"))

my.trees = getTrees(my.clones, build="igphyml", nproc=1, partition="hl",
                   exec="/usr/local/share/igphyml/src/igphyml")
saveRDS(my.trees, file = here::here("analysis","data_objects","10_dowser","HL_trees_Timepoint_igphyml.rds"))

### RAxML make plots ?
#trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","HL_trees_Timepoint_raxml.rds"))

plots <- plotTrees(trees, tips="Timepoint.num",tip_palette = "RdYlBu", tipsize = 3)
plots <- lapply(plots, function(x){
    x + geom_tiplab(aes(label = Timepoint.wk),offset = 0.01)
})

# titles <- trees %>% select(clone_id,vac_grp,Specificity) %>% distinct(clone_id,.keep_all = TRUE)
# for(i in 1:nrow(titles)){
#         plots[[i]] <- plots[[i]] + ggtitle(paste0(titles$clone_id[i]," ",titles$vac_grp[i]," ",titles$Specificity[i]))
# }
treesToPDF(plots, here::here("analysis","plots","10_dowser","HL_tree_plots_timepoint_raxml.pdf"), nrow = 1, ncol = 2, height = 20)

#  [[3]]
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_1401_42_1.pdf"), height = 8)
plots[[23]] + guides(color = guide_legend(title = "Days since Vaccination"))  
dev.off()
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_1401_42_2.pdf"), height = 4)
plots[[115]] + guides(color = guide_legend(title = "Days since Vaccination"))  
dev.off()

# 
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_1605_1_1.pdf"), height = 12)
plots[[12]] + guides(color = guide_legend(title = "Days since Vaccination"))  
dev.off()

# 
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_1064_35_1.pdf"), height = 8)
plots[[21]] + guides(color = guide_legend(title = "Days since Vaccination"))  #318_1064_35
dev.off()

#  [[5]]
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_1091_235_1.pdf"), height = 8)
plots[[82]] + guides(color = guide_legend(title = "Days since Vaccination"))  # 415_1091_235
dev.off()

#  [[5]]
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_1091_235_2.pdf"), height = 8)
plots[[89]] + guides(color = guide_legend(title = "Days since Vaccination"))  # 415_1091_235
dev.off()

#  HL
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_1039_34_1.pdf"), height = 12)
plots[[38]] + guides(color = guide_legend(title = "Days since Vaccination"))  # 415_1091_235
dev.off()


#  H
pdf(file = here::here("analysis","plots","10_dowser","tree_plots_phyigml_timepoint_manuscript_#phi_1039_34.pdf"), height = 12)
plots[[69]] + guides(color = guide_legend(title = "Days since Vaccination"))  # 415_1091_235
dev.off()

#  
pdf(file = here::here("analysis","plots","10_dowser","tree_plots_phyigml_timepoint_manuscript_#phi_340_88 .pdf"), height = 10)
plots[[122]] + guides(color = guide_legend(title = "Days since Vaccination"))  # 415_1091_235
dev.off()

#  
pdf(file = here::here("analysis","plots","10_dowser","HL_tree_plots_phyigml_timepoint_manuscript_#phi_340_88 .pdf"), height = 10)
plots[[60]] + guides(color = guide_legend(title = "Days since Vaccination"))  # 415_1091_235
dev.off()

# RAxML evolution analysis

### measuring evolution ###
# do i need to collapse trees before doing this?
evo.trees <- correlationTest(trees,time = "Timepoint.num",permutations = 10000,perm_type = "uniform") # 
evo.trees <- evo.trees[order(evo.trees$p),]

saveRDS(evo.trees,file = here::here("analysis","data_objects","10_dowser","HL_Evo_Trees_igphyml_Time.rds"))
dat <- evo.trees %>% select(clone_id,seqs,slope,p,correlation)
write.csv(dat,file = here::here("analysis","data_objects","10_dowser","HL_Evo_dat_Time.csv"))