# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# cd /data/vrc_bsc/Abby/Experiment_316
# srun --pty bash
# apptainer shell --bind /data/vrc_bsc:/mnt ../../containers/immcantation/immcantation_suite-4.5.0.sif
# R

library(dowser)
library(dplyr)

### 6 all no switching constraints #####
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_Cluster_meta.rds"))

#get rid of trees with less than 10 sequences
trees <- trees %>% filter(seqs > 9) ##removes 20 lineages

trait = "Cluster.res0.4"
igphyml_location = "/usr/local/share/igphyml/src/igphyml"

###run find switches on everything together 
switches <- findSwitches(trees,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser","switches_cluster_all.rds"))

### divide trees up by broad specificity
trees.ls <- split(trees,f = trees$broad_specificity)

switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_cluster.rds"))

# sp = lapply(switches.ls, function(x){
#   testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
# }) # decide here whether to permute among trees or within trees. default is within trees.
# saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster.rds"))

# sp.means <- lapply(sp, function(x){
#   x$means %>% arrange(PGT)
# })

# ## combine list of dataframes into one data frame and add broad specificity variable 
# sp.means <- bind_rows(sp.means, .id = "broad_specificity")
# write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster.csv"))


### divide trees up by broad specificity AND exposure
trees$specificity_exposure <- paste0(trees$broad_specificity,"_",trees$Exposure)
trees.ls <- split(trees,f = c(trees$specificity_exposure))

switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_cluster_specificity_exposure.rds"))

# sp = lapply(switches.ls, function(x){
#   testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
# }) # decide here whether to permute among trees or within trees. default is within trees.
# saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure.rds"))

# sp.means <- lapply(sp, function(x){
#   x$means %>% arrange(PGT)
# })

# ## combine list of dataframes into one data frame and add broad specificity variable 
# sp.means <- bind_rows(sp.means, .id = "specificity_exposure")
# write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure.csv"))

## Now do it by lineage in next script


#### make heatmaps in next script outside of the appatainer becuase I can't install pheatmap in the container.

### divide trees up by broad specificity AND vac group
trees$specificity_vac_grp <- paste0(trees$broad_specificity,"_",trees$vac_grp)
trees.ls <- split(trees,f = c(trees$specificity_vac_grp))

switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_cluster_specificity_vac_grp.rds"))


# ### do SP test across dataset
# sp = lapply(switches.ls, function(x){
#   testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
# }) # decide here whether to permute among trees or within trees. default is within trees.
# saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure.rds"))

# sp.means <- lapply(sp, function(x){
#   x$means %>% arrange(PGT)
# })

# ## combine list of dataframes into one data frame and add broad specificity variable 
# sp.means <- bind_rows(sp.means, .id = "specificity_exposure")
# write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure.csv"))

