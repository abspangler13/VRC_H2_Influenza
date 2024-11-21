# cd /data/vrc_bsc/Abby/Experiment_316
# srun --pty bash
# apptainer shell --bind /data/vrc_bsc:/mnt ../../containers/immcantation/immcantation_suite-4.5.0.sif
# R
library(dowser)
library(vroom)
library(dplyr)
library(tidyr)
library(ggtree) 


# H igphyml 
H_igphyml <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_Time_meta_filtered.rds"))

# HL igphyml
HL_igphyml <- readRDS(file = here::here("analysis","data_objects","10_dowser","HL_trees_Timepoint_igphyml.rds"))
HL_igphyml2 <- readRDS(file = here::here("analysis","data_objects","10_dowser","HL_trees_Timepoint_igphyml_2.rds"))
HL_igphyml <-bind_rows(HL_igphyml,HL_igphyml2)

# HL raxml 
HL_raxml <- readRDS(file = here::here("analysis","data_objects","10_dowser","HL_trees_Timepoint_raxml.rds"))
HL_raxml2 <- readRDS(file = here::here("analysis","data_objects","10_dowser","HL_trees_Timepoint_raxml_2.rds"))
HL_raxml <-bind_rows(HL_raxml,HL_raxml2)

dim(H_igphyml)
# [1] 139  11
length(unique(H_igphyml$clone_id))
#[1] 139

dim(HL_igphyml)
# [1] 181   6

dim(HL_raxml)
# [1] 181   6

split_HL_ig <- c() ## H clone id that split in HL igphyml model
split_HL_ig_idx <- c() ## HL clone indices
split_HL_rax <- c() # H clone id that split in HL raxml model
split_HL_rax_idx <- c() ## HL clone indices

for(i in 1:length(H_igphyml$clone_id)){
  x <- grep(H_igphyml$clone_id[i],HL_igphyml$clone_id)
  if(length(x) > 1){
    split_HL_ig[i] <- H_igphyml$clone_id[i]
    split_HL_ig_idx <- c(split_HL_ig_idx,x[1],x[2])
  }
  y <- grep(H_igphyml$clone_id[i],HL_raxml$clone_id)
  if(length(y) > 1){
    split_HL_rax[i] <- H_igphyml$clone_id[i]
    split_HL_rax_idx <- c(split_HL_rax_idx,y[1],y[2])
  }
}
split_HL_ig <- split_HL_ig[!is.na(split_HL_ig)]
split_HL_rax <- split_HL_rax[!is.na(split_HL_rax)]

#subset to only the clones that split
HL_igphyml <- HL_igphyml[split_HL_ig_idx,]

##subset H object to only clones that split
H_igphyml <- H_igphyml %>% filter(clone_id %in% split_HL_ig)

#### un-do airrclone objects
unlist.dat <- data.frame()
for(i in 1:nrow(HL_igphyml)){
    dat <- HL_igphyml[i,2]$data[[1]]@data
    dat$clone_id <- HL_igphyml$clone_id[i]
    unlist.dat <- rbind(unlist.dat,dat)
}

### check number of sequences for each sub clone after split. Should add up to the number of sequences in the original clone. 
write.csv(unlist.dat,file = here::here("analysis","data_objects","10_dowser","unlisted_HL_split_clones_igphyml.csv"))

#### un-do airrclone objects
unlist.dat.H <- data.frame()
for(i in 1:nrow(H_igphyml)){
    dat <- H_igphyml[i,1]$data[[1]]@data
    dat$clone_id <- H_igphyml$clone_id[i]
    unlist.dat.H <- rbind(unlist.dat.H,dat)
}

### check number of sequences for each sub clone after split. Should add up to the number of sequences in the original clone. 
write.csv(unlist.dat.H,file = here::here("analysis","data_objects","10_dowser","unlisted_H_split_clones_igphyml.csv"))


### add the CDRH3/CDRL3   HC and LC the v_call, d_call, j_call, c_call,  junction_aa.  You can add the mu_freq as well
## write.csv(All.vdj.M, file = here::here("analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))
seurat <- read.csv(file = here::here("analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))
seurat <- seurat %>% select(v_call,d_call,j_call,c_call,junction_aa,LC_v_call,LC_j_call,LC_c_call,LC_junction_aa,mu_freq,CELL) ## which LC columns to pick?
colnames(seurat)[grep("CELL",colnames(seurat))] <- "cell_id"

split.clones <- read.csv(file = here::here("analysis","data_objects","10_dowser","unlisted_HL_split_clones_igphyml.csv"))

# > dim(split.clones)
# [1] 1036   16
# > length(unique(split.clones$cell_id))
# [1] 1036

split.clones <- left_join(split.clones,seurat,by = "cell_id")

write.csv(split.clones, file = here::here("analysis","data_objects","10_dowser","unlisted_split_clones_igphyml.csv"))
