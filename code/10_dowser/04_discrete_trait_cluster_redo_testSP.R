# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# cd /data/vrc_bsc/Abby/Experiment_316
# srun --pty bash
# apptainer shell --bind /data/vrc_bsc:/mnt ../../containers/immcantation/immcantation_suite-4.5.0.sif
# R

library(dowser)
library(dplyr)


### run sp test for all data together
switches <- readRDS(file = here::here("analysis","data_objects","10_dowser","switches_cluster_all.rds"))
sp <- testSP(switches$switches, alternative="greater",bylineage = TRUE,permuteAll = FALSE)
sp.means <- sp$means %>% arrange(PGT)
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_by_lineage_all.csv"))


## by lineage for just specificity
switches.ls <- readRDS(file = here::here("analysis","data_objects","10_dowser","switches_cluster.rds"))

sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = TRUE,permuteAll = FALSE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_by_lineage.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add broad specificity variable 
sp.means <- bind_rows(sp.means, .id = "broad_specificity")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_by_lineage.csv"))


### divide trees up by broad specificity AND exposure
switches.ls <- readRDS(file = here::here("analysis","data_objects","10_dowser","switches_cluster_specificity_exposure.rds"))

sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = TRUE,permuteAll = FALSE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure_by_lineage.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add broad specificity variable 
sp.means <- bind_rows(sp.means, .id = "specificity_exposure")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure_by_lineage.csv"))


### divide trees up by broad specificity AND vaccine group CHANGE permuteAll = FALSE
switches.ls <- readRDS(file = here::here("analysis","data_objects","10_dowser","switches_cluster_specificity_vac_grp.rds"))

sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = TRUE,permuteAll = FALSE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_vac_grp_by_lineage.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add broad specificity variable 
sp.means <- bind_rows(sp.means, .id = "specificity_vac_grp")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_vac_grp_by_lineage.csv"))

