# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# cd /data/vrc_bsc/Abby/Experiment_316
# srun --pty bash
# apptainer shell --bind /data/vrc_bsc:/mnt ../../containers/immcantation/immcantation_suite-4.5.0.sif
# R

library(dowser)
library(pheatmap)
library(dplyr)

### 6 all no switching constraints #####
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_All_meta_filtered.rds"))

trait = "Cluster.res0.4"
igphyml_location = "/usr/local/share/igphyml/src/igphyml"

# trees <- trees %>% mutate(exposure = case_when(vac_grp == "3A" ~ "Naive",
#                                                         vac_grp == "3B" ~ "Exposed",
#                                                         vac_grp == "4A" ~ "Naive",
#                                                         vac_grp == "4B" ~ "Exposed"))

##filter for H2 only
trees <- trees %>% filter(broad_specificity == "H2_Only")

### divide trees up by exposure 
trees.ls <- split(trees,f = trees$exposure)

switches.ls <- lapply(trees.ls,function(x){
  findSwitches(x,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
})

saveRDS(switches.ls, file = here::here("analysis","data_objects","10_dowser","switches_cluster_naive_exposed_H2.rds"))

# sp = lapply(switches.ls, function(x){
#   testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = FALSE)
# }) # decide here whether to permute among trees or within trees. default is within trees.
# saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_within_cluster.rds"))

# sp.means <- lapply(sp, function(x){
#   x$means %>% arrange(PGT)
# })

# ## combine list of dataframes into one data frame and add subject variable 
# sp.means <- bind_rows(sp.means, .id = "Subject")
# write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_within_cluster.csv"))


sp = lapply(switches.ls, function(x){
  testSP(x$switches, alternative="greater",bylineage = FALSE,permuteAll = TRUE)
}) # decide here whether to permute among trees or within trees. default is within trees.
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_H2.rds"))

sp.means <- lapply(sp, function(x){
  x$means %>% arrange(PGT)
})

## combine list of dataframes into one data frame and add subject variable 
sp.means <- bind_rows(sp.means, .id = "exposure")
write.csv(sp.means, file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_H2.csv"))


### compare Naive vs Exposed
Breaks <- seq(0, 1, length = 100)

dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_H2.csv"),row.names=1)
dat.n <- dat %>% filter(PGT < 0.05 & exposure == "Naive")
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_naive.pdf"))
pheatmap(table(dat.n$FROM,dat.n$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

dat.e <- dat %>% filter(PGT < 0.05 & exposure == "Exposed")
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_exposed.pdf"))
pheatmap(table(dat.e$FROM,dat.e$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20,breaks=Breaks, display_numbers = T,number_format = "%.0f")
dev.off()


Breaks <- seq(0, 4, length = 100)

dat.a.n <- dat.a %>% filter(Subject %in% #phi)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_naive.pdf"))
pheatmap(table(dat.a.n$FROM,dat.a.n$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20,breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

dat.a.e <- dat.a %>% filter(Subject %in% #phi)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_exposed.pdf"))
pheatmap(table(dat.a.e$FROM,dat.a.e$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20,breaks=Breaks,display_numbers=T,number_format = "%.0f")#
dev.off()


#####
tab <- trees %>% select(clone_id,vac_grp)
colnames(tab) <- c("CLONE","vac_grp")

sp$means <- left_join(sp$means,tab, by = "CLONE")
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_lineage",names[k],".csv")))

sig <- sp$means %>% filter(PGT < 0.05)

table(sig$vac_grp)
table(sig$FROM,sig$TO)
pdf(file = here::here("analysis","plots","10_dowser",paste0("sig_switches_lineage_pheatmap_",names[k],".pdf")))
pheatmap(table(sig$FROM,sig$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()

x <- sig %>% filter(vac_grp == "3A")
pdf(file = here::here("analysis","plots","10_dowser",paste0("sig_switches_lineage_pheatmap_",names[k],"_3A.pdf")))
pheatmap(table(x$FROM,x$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()

x <- sig %>% filter(vac_grp == "3B")
pdf(file = here::here("analysis","plots","10_dowser",paste0("sig_switches_lineage_pheatmap_",names[k],"_3B.pdf")))
pheatmap(table(x$FROM,x$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()

x <- sig %>% filter(vac_grp == "4A")
pdf(file = here::here("analysis","plots","10_dowser",paste0("sig_switches_lineage_pheatmap_",names[k],"_4A.pdf")))
pheatmap(table(x$FROM,x$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()

x <- sig %>% filter(vac_grp == "4B")
pdf(file = here::here("analysis","plots","10_dowser",paste0("sig_switches_lineage_pheatmap_",names[k],"_4B.pdf")))
pheatmap(table(x$FROM,x$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()


# Perform PS and SP tests on switches by repertoire (divide by vaccine group)
grp3A <- trees %>% filter(vac_grp == "3A") %>% select(clone_id)
grp3B <- trees %>% filter(vac_grp == "3B") %>% select(clone_id)
grp4A <- trees %>% filter(vac_grp == "4A") %>% select(clone_id)
grp4B <- trees %>% filter(vac_grp == "4B") %>% select(clone_id)

# ps = testPS(switches$switches)
# print(ps$means)
# saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_by_repertoire",names[k],".rds")))

sp = testSP(switches$switches %>% filter(CLONE %in% grp3A$clone_id), alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_repertoire_3A_",names[k],".rds")))

write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_repertoire_3A_",names[k],".csv")))

sp = testSP(switches$switches %>% filter(CLONE %in% grp3B$clone_id), alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_repertoire_3B_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_repertoire_3B_",names[k],".csv")))

sp = testSP(switches$switches %>% filter(CLONE %in% grp4A$clone_id), alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_repertoire_4A_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_repertoire_4A_",names[k],".csv")))

sp = testSP(switches$switches %>% filter(CLONE %in% grp4B$clone_id), alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_repertoire_4B_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_repertoire_4B_",names[k],".csv")))


##### Naive vs Exposed H2 only for K = 3 c_call
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

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()