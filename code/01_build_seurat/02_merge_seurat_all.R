#qrsh -l mem_free=20G,h_vmem=20G

library(Seurat)
library(here)
library(tidyverse)
library(gridExtra)
library(sessioninfo)    
library(tidyseurat)

run_info <- read.csv(file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"))


######## Run 2 and Run5 #############
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phiB.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phiB","#phiB_pos.rds"))
#phiB.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phiB","#phiB_pos.rds"))

#load all seurat objects containing negative cells for each sample
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phiB.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phiB","#phiB_neg.rds"))
#phiB.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phiB","#phiB_neg.rds"))

#renames probes before merge
##phi.p
HA.df <- GetAssayData(#phi.p, assay = "Probes", slot = "counts")
rownames(HA.df)
rownames(HA.df)<-c("H2-oligo","SA-PE-0951","SA-0972","SA-PE-0953","SA-PE-0954","SA-PE-0955")
rownames(HA.df)
HA.assay<- CreateAssayObject(counts = HA.df)
rownames(HA.assay)
#phi.p[["HAs"]] <- HA.assay
#phi.p[['Probes']] <- NULL
#phi.p <- RenameAssays(#phi.p, HAs = "Probes")

##phi.n
HA.df <- GetAssayData(#phi.n, assay = "Probes", slot = "counts")
rownames(HA.df)
rownames(HA.df)<-c("H2-oligo","SA-PE-0951","SA-0972","SA-PE-0953","SA-PE-0954","SA-PE-0955")
rownames(HA.df)
HA.assay<- CreateAssayObject(counts = HA.df)
rownames(HA.assay)
#phi.n[["HAs"]] <- HA.assay
#phi.n[['Probes']] <- NULL
#phi.n <- RenameAssays(#phi.n, HAs = "Probes")

#renames probes before merge
##phi.p
HA.df <- GetAssayData(#phi.p, assay = "Probes", slot = "counts")
rownames(HA.df)
rownames(HA.df)<-c("H2-oligo","SA-PE-0951","SA-0972","SA-PE-0953","SA-PE-0954","SA-PE-0955")
rownames(HA.df)
HA.assay<- CreateAssayObject(counts = HA.df)
rownames(HA.assay)
#phi.p[["HAs"]] <- HA.assay
#phi.p[['Probes']] <- NULL
#phi.p <- RenameAssays(#phi.p, HAs = "Probes")

#renames probes before merge
##phi.p
HA.df <- GetAssayData(#phi.n, assay = "Probes", slot = "counts")
rownames(HA.df)
rownames(HA.df)<-c("H2-oligo","SA-PE-0951","SA-0972","SA-PE-0953","SA-PE-0954","SA-PE-0955")
rownames(HA.df)
HA.assay<- CreateAssayObject(counts = HA.df)
rownames(HA.assay)
#phi.n[["HAs"]] <- HA.assay
#phi.n[['Probes']] <- NULL
#phi.n <- RenameAssays(#phi.n, HAs = "Probes")

#merge objects into one and save
A316.p2 <- merge(#phi.p, y = c(#phi.p,#phiB.p,#phiB.p), 
                add.cell.ids = c("#phi", "#phi","#phiB","#phiB"))
plyr::count(A316.p2@meta.data$orig.ident)
#       x freq
# 1  #phi 1#phi
# 2 #phiB 2952
# 3  #phi 1457
# 4 #phiB  741
#merge into one object and save
A316.n2 <- merge(#phi.n, y = c(#phi.n,#phiB.n,#phiB.n), 
                add.cell.ids = c("#phi", "#phi","#phiB","#phiB"))
             
A316.p2$Subject <- "#phi"

A316.p2 <- A316.p2%>%mutate(Timepoint = case_when(MULTI_ID == "HTO-0251" & orig.ident %in% c("#phiB","#phiB") ~ "HA-",
                                                       MULTI_ID == "HTO-0254" & orig.ident %in% c("#phiB","#phiB") ~ "B180", 
                                                       MULTI_ID == "HTO-0260" & orig.ident %in% c("#phiB","#phiB") ~ "B090",
                                                       MULTI_ID == "HTO-0259" & orig.ident %in% c("#phiB","#phiB") ~ "B028",
                                                       MULTI_ID == "HTO-0258" & orig.ident %in% c("#phiB","#phiB") ~ "B014",
                                                       MULTI_ID == "HTO-0257" & orig.ident %in% c("#phiB","#phiB") ~ "B007",
                                                       MULTI_ID == "HTO-0256" & orig.ident %in% c("#phiB","#phiB") ~ "B000",
                                                       MULTI_ID == "HTO-0252" & orig.ident %in% c("#phiB","#phiB") ~ "A000",

                                                       MULTI_ID == "HTO-0257" & orig.ident %in% c("#phi","#phi") ~ "B007",
                                                       MULTI_ID == "HTO-0258" & orig.ident %in% c("#phi","#phi") ~ "B014",
                                                       MULTI_ID == "HTO-0259" & orig.ident %in% c("#phi","#phi") ~ "B028",
                                                       MULTI_ID == "HTO-0260" & orig.ident %in% c("#phi","#phi") ~ "B090",
                                                       MULTI_ID == "HTO-0252" & orig.ident %in% c("#phi","#phi") ~ "HA-"))


run_info$neg.drops[2] <- dim(A316.n2)[2]
run_info$wt.cells[2] <- dim(A316.p2)[2]

saveRDS(A316.p2, file = here::here("analysis","data_objects","01_build_seurat","A316_p_run2.rds"))
saveRDS(A316.n2, file = here::here("analysis","data_objects","01_build_seurat","A316_n_run2.rds"))
#A316.p2 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_p_run2.rds"))
#A316.n2 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_n_run2.rds"))

######## Run 3 ############
#load all seurat objects containing positive cells from each sample
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))


#merge objects into one and save
A316.p3 <- merge(#phi.p , y = c(#phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p),
                add.cell.ids = c("#phi","#phi", "#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi"))

plyr::count(A316.p3@meta.data$orig.ident)
#        x freq
# 1   #phi 4863
# 2  #phi 6376
# 3  #phi 6061
# 4   #phi 6436
# 5  #phi 3744
# 6  #phi 3356
# 7   #phi 7771
# 8   #phi #phi0
# 9   #phi  504
# 10  #phi 2876
# 11  #phi 1485
# 12  #phi 2286
# 13  #phi 2279
# 14  #phi  631


run_info$wt.cells[3] <- dim(A316.p3)[2]

#load all seurat objects containing negative cells for each sample
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))


#merge into one object and save
A316.n3 <- merge(#phi.n, y = c(#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n), 
                add.cell.ids = c("#phi", "#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi"))


run_info$neg.drops[3] <- dim(A316.n3)[2]

A316.p3 <-A316.p3 %>% tidyseurat::mutate(Subject = substr(strtrim(A316.p3$orig.ident,4),2,4)) 

A316.p3 <- A316.p3%>%mutate(Timepoint = case_when(MULTI_ID == "HTO-0251" ~ "HA-", 
                                                       MULTI_ID == "HTO-0257" ~ "B007",
                                                       MULTI_ID == "HTO-0258" ~ "B014",
                                                       MULTI_ID == "HTO-0259" ~ "B028",
                                                       MULTI_ID == "HTO-0260" ~ "B090",
                                                       MULTI_ID == "HTO-0252" ~ "A000",
                                                       MULTI_ID == "HTO-0253" ~ "A007",
                                                       MULTI_ID == "HTO-0255" ~ "A028",
                                                       MULTI_ID == "HTO-0256" ~ "B000"))

A316.p3 <- A316.p3%>%mutate(vac_grp = case_when(Subject == "#phi" ~ "4A", 
                                                Subject == "#phi" ~ "4B",
                                                Subject == "#phi" ~ "4A",
                                                Subject == "#phi" ~ "4A",
                                                Subject == "#phi" ~ "4B",
                                                Subject == "#phi" ~ "4B"
                                                ))
#A316.p3
HA.df <- GetAssayData(A316.p3, assay = "Probes", slot = "counts")
rownames(HA.df)
rownames(HA.df)<-c("H2-oligo","SA-PE-0951","SA-0972","SA-PE-0953","SA-PE-0954","SA-PE-0955")
rownames(HA.df)
HA.assay<- CreateAssayObject(counts = HA.df)
rownames(HA.assay)
A316.p3[["HAs"]] <- HA.assay
A316.p3[['Probes']] <- NULL
A316.p3 <- RenameAssays(A316.p3, HAs = "Probes")

#A316.n3
HA.df <- GetAssayData(A316.n3, assay = "Probes", slot = "counts")
rownames(HA.df)
rownames(HA.df)<-c("H2-oligo","SA-PE-0951","SA-0972","SA-PE-0953","SA-PE-0954","SA-PE-0955")
rownames(HA.df)
HA.assay<- CreateAssayObject(counts = HA.df)
rownames(HA.assay)
A316.n3[["HAs"]] <- HA.assay
A316.n3[['Probes']] <- NULL
A316.n3 <- RenameAssays(A316.n3, HAs = "Probes")

saveRDS(A316.p3, file = here::here("analysis","data_objects","01_build_seurat","A316_p_run3.rds"))
saveRDS(A316.n3, file = here::here("analysis","data_objects","01_build_seurat","A316_n_run3.rds"))
#A316.p3 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_p_run3.rds"))
#A316.n3 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_n_run3.rds"))

######### Run 4 and Run5 ###############
#load all seurat objects containing positive cells from each sample
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi#phi","#phi#phi_pos.rds"))
#phiB.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phiB","#phiB_pos.rds"))

#merge objects into one and save
A316.p4 <- merge(#phi.p, y = c(#phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi.p, #phi#phi.p,#phiB.p), 
                add.cell.ids = c("#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi","#phi","#phi","#phi#phi","#phiB"))
plyr::count(A316.p4@meta.data$orig.ident)
#          x freq
# 1     #phi 1963
# 2  #phi#phi 7150
# 3     #phi 2539
# 4    #phiB 1223
# 5     #phi 4003
# 6     #phi 5015
# 7     #phi 4925
# 8     #phi 4363
# 9     #phi 5709
# 10    #phi 5437
# 11    #phi 5180
# 12    #phi  746
# 13    #phi 3616
# 14    #phi 1271
# 15    #phi #phi0
# 16    #phi 4602
# 17    #phi 5150
# 18    #phi 1904

run_info$wt.cells[4] <- dim(A316.p4)[2]

#load all seurat objects containing negative cells for each sample
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi#phi","#phi#phi_neg.rds"))
#phiB.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phiB","#phiB_neg.rds"))

#merge into one object and save
A316.n4 <- merge(#phi.n, y = c(#phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi.n, #phi#phi.n,#phiB.n), 
                add.cell.ids = c("#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi", "#phi","#phi","#phi","#phi#phi","#phiB"))


run_info$neg.drops[4] <- dim(A316.n4)[2]

# Label cells based on Pool and Hash Tag #
#A316.p$Subject = A316.p$orig.ident
A316.p4 <- A316.p4%>%mutate(Subject = case_when(orig.ident == "#phi" & MULTI_ID == "HTO-0256" ~ "#phi", 
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0255" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0254" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0251" ~ "M#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0258" ~ "M#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0257" ~ "M#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0259" ~ "M#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0260" ~ "M#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0256" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0255" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0254" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0251" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0258" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0257" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0259" ~ "#phi",
                                                orig.ident == "#phi" & MULTI_ID == "HTO-0260" ~ "#phi",
                                                orig.ident == "#phi#phi" & MULTI_ID == "HTO-0252" ~ "#phiB",
                                                orig.ident == "#phi#phi" & MULTI_ID == "HTO-0256" ~ "#phiB",
                                                orig.ident == "#phi#phi" & MULTI_ID == "HTO-0253" ~ "#phiB",
                                                orig.ident == "#phi#phi" & MULTI_ID == "HTO-0255" ~ "#phiB",
                                                orig.ident == "#phi#phi" & MULTI_ID == "HTO-0254" ~ "#phiB",
                                                orig.ident == "#phi#phi" & MULTI_ID == "HTO-0258" ~ "#phiB",
                                                orig.ident == "#phi#phi" & MULTI_ID == "HTO-0251" ~ "#phiB",
                                                TRUE ~ orig.ident
                                                ))


A316.p4 <- A316.p4%>%mutate(Timepoint = case_when(MULTI_ID == "HTO-0251" ~ "HA-",
                                                
                                                MULTI_ID == "HTO-0252" ~ "A000",
                                                
                                                MULTI_ID == "HTO-0253" & Subject %in% c("#phiB") ~ "A000",
                                                MULTI_ID == "HTO-0253" & Subject %in% c("#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi") ~ "A007",
                                                
                                                MULTI_ID == "HTO-0254" ~ "B180",
                                                
                                                MULTI_ID == "HTO-0255" & Subject %in% c("#phiB") ~ "B000",
                                                MULTI_ID == "HTO-0255" & Subject %in% c("#phi","#phi")  ~ "B014",
                                                MULTI_ID == "HTO-0255" & Subject %in% c("#phi","#phi","#phi","#phi")  ~ "A028",
                                                
                                                MULTI_ID == "HTO-0256" & Subject %in% c("#phi","#phi") ~ "A000",
                                                MULTI_ID == "HTO-0256" & Subject %in% c("#phi","#phi","#phi","#phi","#phiB","#phiB")  ~ "B000",
                                                
                                                MULTI_ID == "HTO-0257" & Subject %in% c("#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi") ~ "B007",
                                                MULTI_ID == "HTO-0257" & Subject %in% c("M#phi","#phi") ~ "B180",

                                                MULTI_ID == "HTO-0258" ~ "B014",

                                                MULTI_ID == "HTO-0259" & Subject %in% c("M#phi","#phi") ~ "B014",
                                                MULTI_ID == "HTO-0259" & Subject %in% c("#phi","#phi","#phi","#phi","#phi","#phi","#phi")  ~ "B028",

                                                MULTI_ID == "HTO-0260" & Subject %in% c("#phi","#phi","#phi","#phi","#phi","#phi","#phi")  ~ "B090",
                                                MULTI_ID == "HTO-0260" & Subject %in% c("M#phi","#phi")  ~ "B180"
                                                ))



A316.p4 <- A316.p4%>%mutate(Subject = substr(A316.p4$Subject,2,4))

#fill in vac_grp
A316.p4 <- A316.p4%>%mutate(vac_grp = case_when(Subject == "#phi" ~ "4A", 
                                                Subject == "#phi" ~ "3B",
                                                Subject == "#phi" ~ "3B",
                                                Subject == "#phi" ~ "3A",
                                                Subject == "#phi" ~ "4B",
                                                Subject == "#phi" ~ "4A",
                                                Subject == "#phi" ~ "4B",
                                                Subject == "#phi" ~ "4B",
                                                Subject == "#phi" ~ "4A",
                                                Subject == "#phi" ~ "4A",
                                                Subject == "#phi" ~ "4A",
                                                Subject == "#phi" ~ "4B",
                                                Subject == "#phi" ~ "4B"
                                                ))

table(A316.p4$vac_grp,useNA = "always")

saveRDS(A316.p4, file = here::here("analysis","data_objects","01_build_seurat","A316_p.rds"))
saveRDS(A316.n4, file = here::here("analysis","data_objects","01_build_seurat","A316_n.rds"))
#A316.p4 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_p.rds"))
#A316.n4 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_n.rds"))

message("finished run 4")

######## Run 5 ############
#load all seurat objects containing positive cells from each sample
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))
#phi.p <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_pos.rds"))

#merge objects into one and save
A316.p5 <- merge(#phi.p , y = c(#phi.p,#phi.p,#phi.p,#phi.p,#phi.p,#phi.p,#phi.p,#phi.p,#phi.p,#phi.p),
                add.cell.ids = c("#phi","#phi", "#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi"))

table(A316.p5@meta.data$orig.ident,A316.p5@meta.data$MULTI_ID)

run_info$wt.cells[5] <- dim(A316.p5)[2]

#load all seurat objects containing negative cells for each sample
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))
#phi.n <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","#phi","#phi_neg.rds"))

#merge into one object and save
A316.n5 <- merge(#phi.n , y = c(#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n,#phi.n),
                add.cell.ids = c("#phi","#phi", "#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi","#phi"))


run_info$neg.drops[5] <- dim(A316.n5)[2]

A316.p5 <-A316.p5 %>% tidyseurat::mutate(Subject = substr(strtrim(A316.p5$orig.ident,4),2,4)) 

A316.p5 <- A316.p5%>%mutate(Timepoint = case_when(MULTI_ID == "HTO-0251" ~ "HA-", 
                                                       MULTI_ID == "HTO-0257" ~ "B007",
                                                       MULTI_ID == "HTO-0258" ~ "B014",
                                                       MULTI_ID == "HTO-0259" ~ "B028",
                                                       MULTI_ID == "HTO-0260" ~ "B090",
                                                       MULTI_ID == "HTO-0252" ~ "A000",
                                                       MULTI_ID == "HTO-0256" ~ "B000",
                                                       MULTI_ID == "HTO-0254" ~ "B180"))


saveRDS(A316.p5, file = here::here("analysis","data_objects","01_build_seurat","A316_p_run5.rds"))
saveRDS(A316.n5, file = here::here("analysis","data_objects","01_build_seurat","A316_n_run5.rds"))
#A316.p5 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_p_run5.rds"))
#A316.n5 <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_n_run5.rds"))

# ### run5 hashtag check
# flow.dat <- read.csv(file = here::here("analysis","data_objects","01_build_seurat","hashtag_sorted_cells_run5.csv"))
# flow.dat$id <- paste0(flow.dat$orig.ident, "_", flow.dat$MULTI_ID)
# sc.dat <- A316.p5@meta.data %>% group_by(orig.ident,MULTI_ID) %>% count()
# sc.dat$id <- paste0(sc.dat$orig.ident,"_",sc.dat$MULTI_ID)
# #use left_join instead
# dat <- full_join(flow.dat,sc.dat, by = "id")
# write.csv(dat,file = here::here("analysis","data_objects","01_build_seurat","hashtag_table_for_sarah.csv"))
# samples <- unique(dat$orig.ident)
# pdf(file = here::here("analysis","plots","01_build_seurat","flow_vs_10x_hashtags.pdf"))
# for(i in 1:length(samples)){
#     my.plot <-dat %>% filter(orig.ident == samples[i]) %>% ggplot(aes(x=MULTI_ID, y=n, fill=assay)) +
#     geom_bar(stat="identity", position=position_dodge()) + ggtitle(samples[i]) +
#     theme_minimal()
#     plot(my.plot)
# }
# dev.off()

# fix HA names for all objects before merge
# seurats <- list(A316.p2,A316.p3,A316.p4,A316.p5,A316.n2,A316.n3,A316.n4,A316.n5)
# seurats <- lapply(seurats, function(i){
#     HA.df <- GetAssayData(i, assay = "Probes", slot = "counts") 
#     rownames(HA.df)<-c("H2FL", "H5FL", "NC99FL", "H3st", "MichFL", "H2st")
#     HA.assay<- CreateAssayObject(counts = HA.df)
#     i[["HAs"]] <- HA.assay
#     i[['Probes']] <- NULL
#     i <- RenameAssays(i, HAs = "Probes")
# })

seurats <- list(A316.p2,A316.p3,A316.p4,A316.p5,A316.n2,A316.n3,A316.n4,A316.n5)

### merge all together
A316_p_all <- merge(seurats[[1]], y = c(seurats[[2]],seurats[[3]],seurats[[4]]))
A316_n_all <- merge(seurats[[5]], y = c(seurats[[6]],seurats[[7]],seurats[[8]]))


### add line to run_info
run_info$wt.cells[5] <- dim(A316_p_all)[2]
run_info$neg.drops[5] <- dim(A316_n_all)[2]
write.csv(run_info,file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"), row.names = FALSE)

saveRDS(A316_p_all, file = here::here("analysis","data_objects","01_build_seurat","A316_p_all.rds"))
#A316_p_all <- readRDS(file = here::here("analysis","data_objects","01_build_seurat","A316_p_all.rds"))
saveRDS(A316_n_all, file = here::here("analysis","data_objects","01_build_seurat","A316_n_all.rds"))

#A316_p_all %>% filter(Subject == "#phi", run == "run4") %>% group_by(Timepoint) %>% dplyr::summarise(Freq = n())

###
# add real HA names instead of SA-PE 
# remove nCount_HA columns from meta.data
# what X column in meta.data?

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()