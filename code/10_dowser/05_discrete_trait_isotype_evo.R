# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# module load immcantation/4.4.0
# singularity shell $IMMCAN_SING
# R

##dowser screen

library(dowser)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

### 6 all no switching constraints #####
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_All_meta_filtered.rds"))
evo.dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","Evo_dat_all_uniform.csv"),row.names=1)
evo.clones <- evo.dat %>% filter(p < 0.05) %>% select(clone_id)
trees$evolving <- ifelse(trees$clone_id %in% evo.clones$clone_id, TRUE, FALSE)

trait = "c_call"
igphyml_location = "/usr/local/share/igphyml/src/igphyml"

### divide trees up by whether or not they're evolvoing
trees.ls <- split(trees,f = trees$evolving)
switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_c_call_evolving.rds"))

sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_evolving.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add subject variable 
sp.means <- bind_rows(sp.means, .id = "evolving")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_evolving.csv"))


### compare H2_cross vs H2_only
Breaks <- seq(0, 1, length = 100)

dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_evolving.csv"),row.names=1)
dat.evo <- dat %>% filter(PGT < 0.05 & evolving == TRUE)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_c_call_evolving.pdf"))
pheatmap(table(dat.evo$FROM,dat.evo$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

dat.not <- dat %>% filter(PGT < 0.05 & evolving == FALSE)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_c_call_not_evolving.pdf"))
pheatmap(table(dat.not$FROM,dat.not$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20,breaks=Breaks, display_numbers = T,number_format = "%.0f")
dev.off()


##now do it for cluster screen 7827.dowser_plots
trait = "Cluster.res0.4"

### divide trees up by whether or not they're evolvoing
switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_cluster_evolving.rds"))

sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_evolving.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add subject variable 
sp.means <- bind_rows(sp.means, .id = "evolving")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_evolving.csv"))


### compare H2_cross vs H2_only
Breaks <- seq(0, 1, length = 100)

dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_evolving.csv"),row.names=1)
dat.evo <- dat %>% filter(PGT < 0.05 & evolving == TRUE)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_evolving.pdf"))
pheatmap(table(dat.evo$FROM,dat.evo$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

dat.not <- dat %>% filter(PGT < 0.05 & evolving == FALSE)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_not_evolving.pdf"))
pheatmap(table(dat.not$FROM,dat.not$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20,breaks=Breaks, display_numbers = T,number_format = "%.0f")
dev.off()