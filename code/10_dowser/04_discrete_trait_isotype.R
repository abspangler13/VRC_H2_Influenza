# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# module load immcantation/4.4.0
# singularity shell $IMMCAN_SING
# R

library(dowser)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

### 6 all no switching constraints #####
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_All_meta_filtered.rds"))

trait = "c_call"
igphyml_location = "/usr/local/share/igphyml/src/igphyml"

### divide trees up by broad specificity 
trees.ls <- split(trees,f = trees$broad_specificity)

switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_c_call_broad_spec.rds"))
#switches.ls <- readRDS(file = here::here("analysis","data_objects","10_dowser","switches_c_call_broad_spec.rds"))

# sp = lapply(switches.ls, function(x){
#   testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = FALSE)
# }) # decide here whether to permute among trees or within trees. default is within trees.
# saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_within_c_call_broad_spec.rds"))

# sp.means <- lapply(sp, function(x){
#   x$means %>% arrange(PGT)
# })

# ## combine list of dataframes into one data frame and add subject variable 
# sp.means <- bind_rows(sp.means, .id = "Subject")
# write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_within_c_call_broad_spec.csv"))


sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_broad_spec.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add subject variable 
sp.means <- bind_rows(sp.means, .id = "broad_specificity")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_broad_spec.csv"))


### compare H2_cross vs H2_only
Breaks <- seq(0, 1, length = 100)

dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_broad_spec.csv"),row.names=1)
dat.h2 <- dat %>% filter(PGT < 0.05 & broad_specificity == "H2_Only")
pdf(file = here::here("analysis","plots","10_dowser","switches_among_c_call_H2.pdf"))
pheatmap(table(dat.h2$FROM,dat.h2$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

dat.cr <- dat %>% filter(PGT < 0.05 & broad_specificity == "H2_Cross")
pdf(file = here::here("analysis","plots","10_dowser","switches_among_c_call_cross.pdf"))
pheatmap(table(dat.cr$FROM,dat.cr$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20,breaks=Breaks, display_numbers = T,number_format = "%.0f")
dev.off()

### compare naive vs exposed
trees <- trees %>% mutate(exposure = case_when(vac_grp == "3A" ~ "Naive",
                                                        vac_grp == "3B" ~ "Exposed",
                                                        vac_grp == "4A" ~ "Naive",
                                                        vac_grp == "4B" ~ "Exposed"))

### divide trees up by broad specificity 
trees.ls <- split(trees,f = trees$exposure)

switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_c_call_exposure.rds"))

sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_exposure.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add subject variable 
sp.means <- bind_rows(sp.means, .id = "exposture")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_exposure.csv"))


### exposure
Breaks <- seq(0, 1, length = 100)

dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_c_call_exposure.csv"),row.names=1)
dat.h2 <- dat %>% filter(PGT < 0.05 & exposture == "Naive")
pdf(file = here::here("analysis","plots","10_dowser","switches_among_c_call_Naive.pdf"))
pheatmap(table(dat.h2$FROM,dat.h2$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

dat.cr <- dat %>% filter(PGT < 0.05 & exposture == "Exposed")
pdf(file = here::here("analysis","plots","10_dowser","switches_among_c_call_exposed.pdf"))
pheatmap(table(dat.cr$FROM,dat.cr$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20,breaks=Breaks, display_numbers = T,number_format = "%.0f")
dev.off()