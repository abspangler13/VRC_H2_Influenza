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

k = 3

igphyml_location = "/usr/local/share/igphyml/src/igphyml"

### 3 c_call no constraints ###
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_c_call_meta.rds"))

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

# Perform PS and SP tests on switches by lineage
ps = testPS(switches$switches,bylineage=TRUE)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_by_lineage",names[k],".rds")))

sp = testSP(switches$switches, alternative="greater",bylineage=TRUE)
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_lineage",names[k],".rds")))

tab <- trees %>% select(clone_id,vac_grp)
colnames(tab) <- c("CLONE","vac_grp")
sp$means <- left_join(sp$means,tab, by = "CLONE")
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_lineage",names[k],".csv")))

sig <- sp$means %>% filter(PGT < 0.05)

table(sig$vac_grp)
table(sig$FROM,sig$TO)
pdf(file = here::here("analysis","plots","10_dowser","sig_switches_lineage.pdf"))
pheatmap(table(sig$FROM,sig$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 20,fontsize_col = 20)
dev.off()

x <- sig %>% filter(vac_grp == "3A")
pdf(file = here::here("analysis","plots","10_dowser","sig_switches_lineage_3A.pdf"))
pheatmap(table(x$FROM,x$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()

x <- sig %>% filter(vac_grp == "3B")
pdf(file = here::here("analysis","plots","10_dowser","sig_switches_lineage_3B.pdf"))
pheatmap(table(x$FROM,x$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()

x <- sig %>% filter(vac_grp == "4A")
pdf(file = here::here("analysis","plots","10_dowser","sig_switches_lineage_4A.pdf"))
pheatmap(table(x$FROM,x$TO),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20)
dev.off()

x <- sig %>% filter(vac_grp == "4B")
pdf(file = here::here("analysis","plots","10_dowser","sig_switches_lineage_4B.pdf"))
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


#### do it for naive vs exposed
# Perform PS and SP tests on switches by exposure
naive <- trees %>% filter(vac_grp %in% c("3A","4A")) %>% select(clone_id)
exposed <- trees %>% filter(vac_grp %in% c("3B","4B")) %>% select(clone_id)

sp = testSP(switches$switches %>% filter(CLONE %in% naive$clone_id), alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_naive_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_naive_",names[k],".csv")))

sp = testSP(switches$switches %>% filter(CLONE %in% exposed$clone_id), alternative="greater")
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_exposed_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_exposed_",names[k],".csv")))


##### 3 incorporate swithcing constraints #### 
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
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_",names[k],"_switch_constraints.rds")))

sp = testSP(switches$switches, alternative="greater")
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_switch_constraints.rds")))
#sp <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],".rds")))
write.csv(sp$means, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_",names[k],"_switch_constraints.csv")))


##### Naive vs Exposed H2 only for K = 3 c_call constrained 
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