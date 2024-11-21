library(Seurat)
library(tidyseurat)
library(tidyverse)


#load 3 and 4 vdj seurats 
A316.vdj  <- readRDS(file=here::here("analysis","data_objects","03_vdj","A316_final_vdj_all.rds"))

subject.counts <- as.data.frame(A316.vdj %>% tidyseurat::count(Subject))
write.csv(subject.counts, file = here::here("analysis","data_objects","subject_counts_run_all.csv"))

# join HA's and select for vaccine group 3
DefaultAssay(object = A316.vdj) <- "Probes"
grp3.HA <- as.data.frame(A316.vdj %>% 
    join_features(features=rownames(A316.vdj)) %>% 
    pivot_wider(names_from = .feature, values_from = .abundance_Probes) %>% 
    filter(vac_grp %in% c("3A","3B") & run %in% c("run3","run4")))


# join HA's and selece for vaccine group 4 
grp4.HA <- as.data.frame(A316.vdj %>% 
    join_features(features=rownames(A316.vdj)) %>% 
    pivot_wider(names_from = .feature, values_from = .abundance_Probes) %>% 
    filter(vac_grp %in% c("4A","4B") & run %in% c("run3","run4")))

DefaultAssay(object = A316.vdj) <- "Prot"
grp3.Prot <- as.data.frame(A316.vdj %>% 
    join_features(features=rownames(A316.vdj)) %>% 
    pivot_wider(names_from = .feature, values_from = .abundance_Prot) %>% 
    filter(vac_grp %in% c("3A","3B") & run %in% c("run3","run4")))


# join HA's and selece for vaccin group 4 
grp4.Prot <- as.data.frame(A316.vdj %>% 
    join_features(features=rownames(A316.vdj)) %>% 
    pivot_wider(names_from = .feature, values_from = .abundance_Prot) %>% 
    filter(vac_grp %in% c("4A","4B")& run %in% c("run3","run4")))


grp3 <- cbind(grp3.HA,grp3.Prot[,c(73:132)])
grp4 <- cbind(grp4.HA,grp4.Prot[,c(73:132)])

write.csv(grp3, file = here::here("analysis","data_objects","vac_group_3_HA_vdj.csv"))
write.csv(grp4, file = here::here("analysis","data_objects","vac_group_4_HA_vdj.csv"))
