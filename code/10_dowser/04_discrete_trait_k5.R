# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# module load immcantation/4.4.0
# singularity shell $IMMCAN_SING
# R

library(dowser)
library(pheatmap)
library(dplyr)

t1 = c("Timepoint.num") 
t2 = c("Cluster.res0.4") 
t3 = c("c_call") 
t4 = c("Timepoint.num","Cluster.res0.4")
t5 = c("Timepoint.num",'c_call') 
t6 = c("Timepoint.num","Cluster.res0.4","c_call")
traits <- list(t1,t2,t3,t4,t5,t6)
names <- c("Time","Cluster","c_call","Time_Cluster","Time_ccall","All")

igphyml_location = "/usr/local/share/igphyml/src/igphyml"

#### 5 c_call and timepoint no switching constraints #####
k = 5
igphyml_location = "/usr/local/share/igphyml/src/igphyml"
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_Time_ccall_meta.rds"))

# remove mAb clones
B180 <- read.csv(file = here::here("analysis","data_objects","10_dowser","clones.for.dowser.csv"))

#subset trees
x <-B180 %>% filter(desription2 %in% c("did before","new")) %>% select(clone_id)
trees <- trees %>% filter(clone_id %in% x$clone_id)

trait = "c_call"
switches = findSwitches(trees, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],".rds")))
#switches <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],".rds")))

# Perform PS test on switches
ps = testPS(switches$switches)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],".rds")))

sp = testSP(switches$switches, alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".csv")))


##### 5 incorporate swithcing constraints #### 
trait = "c_call"
igphyml_location = "/usr/local/share/igphyml/src/igphyml"
isotypes = c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2",
  "IGHG4","IGHE","IGHA2")

my.clones <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("my_clones_",names[k],".rds")))

isotype_counts = unlist(lapply(my.clones$data, function(x)
  length(unique(x@data[[trait]]))))

makeModelFile(file=here::here("analysis","data_objects","10_dowser",paste0("isotype_model_",names[k],".txt")), states=isotypes, 
  constraints="irrev", exceptions=c("IGHD,IGHM"))

trees = getTrees(my.clones[isotype_counts > 1,], trait=trait, igphyml=igphyml_location, palette = "Paired",
  modelfile=here::here("analysis","data_objects","10_dowser",paste0("isotype_model_",names[k],".txt")))

saveRDS(trees,file = here::here("analysis","data_objects","10_dowser",paste0("Trees_igphyml_",names[k],"_switch_constraints.rds")))

trait = "c_call"
switches = findSwitches(trees, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],"_switch_constraints.rds")))

# Perform PS test on switches
ps = testPS(switches$switches)
print(ps$means)

# # A tibble: 1 × 6
#   RECON PERMUTE   PLT DELTA STAT   REPS
#   <dbl>   <dbl> <dbl> <dbl> <chr> <int>
# 1  709.    800.     0 -90.1 PS      100
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],"_switch_constraints.rds")))

sp = testSP(switches$switches, alternative="greater")
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_switch_constraints.rds")))
#sp <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_switch_constraints.csv")))

##### Naive vs Exposed H2 only for K = 5 c_call
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("Trees_igphyml_",names[k],"_meta.rds")))
igphyml_location = "/usr/local/share/igphyml/src/igphyml"

##filter trees to get naive clones
trees.naive <- trees %>% filter(vac_grp %in% c("3A","4A") & Specificity == "H2head")

trait = "c_call"
switches = findSwitches(trees.naive, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],"_naive.rds")))

# Perform PS test on switches
ps = testPS(switches$switches)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],"_naive.rds")))

sp = testSP(switches$switches, alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_naive.rds")))
#sp <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_naive.csv")))


##filter trees to get exposed clones
trees.exposed <- trees %>% filter(vac_grp %in% c("3B","4B") & Specificity == "H2head")

trait = "c_call"
switches = findSwitches(trees.exposed, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],"_exposed.rds")))

# Perform PS test on switches
ps = testPS(switches$switches)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],"_exposed.rds")))

sp = testSP(switches$switches, alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_exposed.rds")))
#sp <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_exposed.csv")))

trees.exposed <- trees %>% filter(vac_grp %in% c("3B","4B") & Specificity == "H2head")

trait = "c_call"
switches = findSwitches(trees.exposed, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],"_exposed_constraints.rds")))

# Perform PS test on switches
ps = testPS(switches$switches)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],"_exposed_constraints.rds")))

sp = testSP(switches$switches, alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_exposed_constraints.rds")))
#sp <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_exposed_constraints.csv")))


##### Naive vs Exposed H2 only for K = 5 c_call constrained
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("Trees_igphyml_",names[k],"_meta.rds")))
igphyml_location = "/usr/local/share/igphyml/src/igphyml"

##filter trees to get naive clones
trees.naive <- trees %>% filter(vac_grp %in% c("3A","4A") & Specificity == "H2head")

trait = "c_call"
switches = findSwitches(trees.naive, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],"_naive_constraints.rds")))

# Perform PS test on switches
ps = testPS(switches$switches)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],"_naive_constraints.rds")))

sp = testSP(switches$switches, alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_naive_constraints.rds")))
#sp <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_naive_constraints.csv")))


##filter trees to get exposed clones
trees.exposed <- trees %>% filter(vac_grp %in% c("3B","4B") & Specificity == "H2head")

trait = "c_call"
switches = findSwitches(trees.exposed, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE)

saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],"_exposed_constraints.rds")))

# Perform PS test on switches
ps = testPS(switches$switches)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],"_exposed_constraints.rds")))

sp = testSP(switches$switches, alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_exposed_constraints.rds")))
#sp <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_exposed_constraints.csv")))