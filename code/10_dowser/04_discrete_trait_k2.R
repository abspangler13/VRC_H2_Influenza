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

k = 2

igphyml_location = "/usr/local/share/igphyml/src/igphyml"

# calculate switches along trees compared to 100 random permutations 
# this may take a while, and can be parallelized using nproc

#### 2 ####
trees <- readRDS(file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_Cluster.rds"))

# remove mAb clones
B180 <- read.csv(file = here::here("analysis","data_objects","10_dowser","clones.for.dowser.csv"))

#subset trees
x <-B180 %>% filter(desription2 %in% c("did before","new")) %>% select(clone_id)
trees <- trees %>% filter(clone_id %in% x$clone_id)

trait = "Cluster.res0.4"
switches = findSwitches(trees, permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
saveRDS(switches, file = here::here("analysis","data_objects","10_dowser",paste0("switches_igphyml_",names[k],".rds")))

# Perform PS and SP tests on switches by lineage
ps = testPS(switches$switches,bylineage=TRUE)
print(ps$means)
saveRDS(ps, file = here::here("analysis","data_objects","10_dowser",paste0("switches_ps_by_lineage",names[k],".rds")))

tab <- trees %>% select(clone_id,vac_grp)
colnames(tab) <- c("CLONE","vac_grp")

sp = testSP(switches$switches, alternative="greater",bylineage=TRUE)
print(sp$means)
saveRDS(sp, file = here::here("analysis","data_objects","10_dowser",paste0("switches_sp_by_lineage",names[k],".rds")))

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
