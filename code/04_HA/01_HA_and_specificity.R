library(here)
library(ggplot2)
library(plyr)
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(sessioninfo)

#load data after vdj data has been added and we've filtered out cells that don't have a vdj. 
A316 <- readRDS(file = here::here("analysis","data_objects","035_integration","integrated_A316.rds"))

pdf(file = here::here("analysis","plots","04_HA","CD19_v_CD56_all.pdf"))
A316 %>%
  tidyseurat::join_features(features = c("P-CD19","P-CD56")) %>%
  pivot_wider(names_from = .feature, values_from = .abundance_Prot) %>%
  ggplot(aes(x= `P-CD19`, y=`P-CD56`)) + 
  geom_point(aes(color = Timepoint), size=0.5) +
  scale_x_continuous(limits = c(-2, 30)) +
  scale_y_continuous(limits = c(-2, 20)) +
  facet_wrap(.~Timepoint) +
  theme(aspect.ratio=1)
dev.off()

# Make table of the number of cells present at each time point for each subject 
Cell.Subject <- A316 %>% group_by (Timepoint, Subject) %>% 
  dplyr::summarise(Freq = n()) %>% 
  pivot_wider(names_from = Timepoint, values_from = Freq)

write.csv(Cell.Subject,file = here::here("analysis","data_objects","04_HA","Cells_per_timepoint_all.csv"))

############################################# Look at HA data and add in specificity ################################

#change Probe Names
HA.df <- GetAssayData(A316, assay = "Probes", slot = "counts")
rownames(HA.df)<-c("H2FL", "H5FL", "NC99FL", "H3st", "MichFL", "H2st")
HA.assay<- CreateAssayObject(counts = HA.df)
rownames(HA.assay)
A316[["HAs"]] <- HA.assay
A316[['Probes']] <- NULL
A316 <- RenameAssays(A316, HAs = "Probes")

#re-save object 
saveRDS(A316, file = here::here("analysis","data_objects","04_HA","A316_HA.rds"))
# A316 <- readRDS(file = here::here("analysis","data_objects","04_HA","A316_HA.rds"))

HAs <-rownames(A316@assays$Probes)
HA.combos <- combn(unique(HAs),2)
# HA.combos <-expand.grid(x = HAs,y = HAs)
# HA.combos <- data.frame(lapply(HA.combos, as.character),stringsAsFactors=FALSE)

# deal with counts so probably don't have to adjust scale
pdf(file = here::here("analysis","plots","04_HA","HA_scatterplots_all.pdf"))
for(i in 1:ncol(HA.combos)){
    my.plot <- A316 %>% 
      tidyseurat::join_features(features = c("H2FL", "H5FL", "NC99FL", "H3st", "MichFL", "H2st")) %>%
      pivot_wider(names_from=.feature,values_from=.abundance_Probes) %>%
      ggplot((aes_string(x = HA.combos[1,i], y = HA.combos[2,i]))) + 
        geom_point(aes(color = Timepoint), size = 0.2, ratio = 1) + 
        theme_classic() +
        scale_x_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 16, 32, 64, 128, 256, 5000)) + 
        scale_y_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 16, 32, 64, 128, 256, 5000)) +
        theme(legend.position = "right", 
            axis.title.x = element_text(size = 10), 
            axis.title.y = element_text(size = 10), 
            axis.text.x = element_text(size = 8), 
            axis.text.y = element_text(size = 8), 
            strip.text.x = element_text(size = 8)) + 
        facet_wrap(.~Timepoint) + 
        theme(aspect.ratio=1)
    plot(my.plot)
}
dev.off()

pdf(file = here::here("analysis","plots","04_HA","H5FL_NC99L_all.pdf"))
A316 %>% tidyseurat::join_features(features = c("H2FL", "H5FL", "NC99FL", "H3st", "MichFL", "H2st")) %>%
  pivot_wider(names_from=.feature,values_from=.abundance_Probes) %>%
  filter(H2FL > 32) %>%
  ggplot(aes(x = H5FL + 1, y = NC99FL + 1)) + 
  geom_point(aes(color = MULTI_ID), size = 0.5, ratio = 1) + 
  theme_classic() +
  scale_x_continuous(trans = "log2", limits = c(1, 10000)) + 
  scale_y_continuous(trans = "log2", limits = c(1, 10000)) +
  theme(legend.position = "right", 
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        strip.text.x = element_text(size = 8)) + 
  facet_wrap(.~MULTI_ID) + 
  xlab("H5FL") + 
  ylab("H1Mich") +
  theme(aspect.ratio=1)
dev.off()

#use cell selector to select populations of cells I want. Separate into membcells, plasmablasts and junk
#then save vectors of cell names for each category and transfer them to locus

#load vectors of cell names for each population
cell_selector <- readRDS(file = here::here("analysis","data_objects","04_HA","cell_selector.rds"))

#use mutate function to add broad celltype column to meta.data
A316 <- A316 %>% tidyseurat::mutate(broad.celltype = case_when(CELL %in% cell_selector[[1]] ~ "MemBCell",
                                                                       CELL %in% cell_selector[[2]] ~ "PlasmaBlast",
                                                                       CELL %in% cell_selector[[3]] ~ "Junk"))

#remove junk 817 cells
A316 <- subset(A316, subset = broad.celltype %in% c("MemBCell","PlasmaBlast"))

### Need graph for every HA with Mem and PB populations ####
pdf(file = here::here("analysis","plots","04_HA","Timepoints_All_HAs_M_all.pdf"))
for(i in 1:length(HAs)){
  my.plot <- A316 %>% tidyseurat::join_features(features=HAs) %>%
    pivot_wider(names_from=.feature,values_from=.abundance_Probes) %>%
    filter(broad.celltype == "MemBCell") %>% ggplot(aes_string(x = "Timepoint", y = HAs[i])) + 
    geom_jitter(aes(color = Timepoint), size = 0.3, ratio = 1) + 
    theme_classic() +
    scale_y_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 12, 16, 24, 32, 64, 128, 256, 5000)) +
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8), 
          strip.text.x = element_text(size = 8)) + 
    xlab("Timepoint") + 
    theme(aspect.ratio=1)
  plot(my.plot)
}
dev.off()

pdf(file = here::here("analysis","plots","04_HA","Timepoints_All_HAs_P_all.pdf"))
for(i in 1:length(HAs)){
  my.plot <- A316 %>% tidyseurat::join_features(features=HAs) %>%
    pivot_wider(names_from=.feature,values_from=.abundance_Probes) %>%
    filter(broad.celltype == "PlasmaBlast") %>% ggplot(aes_string(x = "Timepoint", y = HAs[i])) + 
    geom_jitter(aes(color = Timepoint), size = 0.3, ratio = 1) + 
    theme_classic() +
    scale_y_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 12, 16, 24, 32, 64, 128, 256, 5000)) +
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8), 
          strip.text.x = element_text(size = 8)) + 
    xlab("Timepoint") + 
    theme(aspect.ratio=1)
  plot(my.plot)
}
dev.off()

##### Add in specificity based on cutoff determined by graphs ########

### call specificity for Memory bcells only first. For subjects 
odd.subjects <- c("307","313","314")

Specificity.M <- A316 %>% tidyseurat::join_features(features=HAs) %>% filter(broad.celltype == "MemBCell" & Subject %in% setdiff(unique(A316$Subject),odd.subjects)) %>%
pivot_wider(names_from = .feature, values_from = .abundance_Probes) %>% 
tidyseurat::mutate(Specificity = case_when(MichFL <= 16 & NC99FL <= 24 & H2st <= 16 & H5FL <= 32 & H2FL > 24 ~ "H2head", 
                                            MichFL > 16 & H2st < 16 & H2FL > 24 ~ "H2H1 Head",
                                            NC99FL > 24 & H2st < 16 & H2FL > 24  ~ "H2H1 Head",
                                            H5FL > 32 & H2st < 16 & H2FL > 8 ~ "H2H5 Head",
                                            H2st > 16 & H3st <= 8 ~ "H2st",
                                            H2st > 16  & H3st > 8 ~"H2H3st",
                                            H2FL <= 24 & H2st <= 16 ~"H2neg",
                                            TRUE~"unclear")) %>%
                                            select("Specificity",".cell")
colnames(Specificity.M) <- c("Specificity","CELL")

Specificity.M.2 <- A316 %>% tidyseurat::join_features(features=HAs) %>% filter(broad.celltype == "MemBCell" & Subject %in% odd.subjects) %>%
pivot_wider(names_from = .feature, values_from = .abundance_Probes) %>% 
tidyseurat::mutate(Specificity = case_when(MichFL <= 16 & NC99FL <= 24 & H2st <= 16 & H5FL <= 32 & H2FL > 64 ~ "H2head", 
                                 MichFL > 16 & H2st < 16 & H2FL > 64 ~ "H2H1 Head",
                                 NC99FL > 24 & H2st < 16 & H2FL > 64  ~ "H2H1 Head",
                                 H5FL > 32 & H2st < 16 & H2FL > 8 ~ "H2H5 Head",
                                 H2st > 16 & H3st <= 8 ~ "H2st",
                                 H2st > 16  & H3st > 8 ~"H2H3st",
                                 H2FL <= 64 & H2st <= 16 ~"H2neg",
                                 TRUE~"unclear")) %>%
                                 select("Specificity",".cell")
colnames(Specificity.M.2) <- c("Specificity","CELL")

Specificity.M <- rbind(Specificity.M,Specificity.M.2)

### call specificity for plasmablasts next 
Specificity.P <- A316 %>% tidyseurat::join_features(features=HAs) %>% filter(broad.celltype == "PlasmaBlast") %>%
pivot_wider(names_from = .feature, values_from = .abundance_Probes) %>% 
tidyseurat::mutate(Specificity = case_when(MichFL <= 12 & NC99FL <= 12 & H2st <= 4 & H5FL <= 16 & H2FL > 8 ~ "H2head", 
                                    MichFL > 12 & H2st < 4 & H2FL > 8 ~ "H2H1 Head",
                                    NC99FL > 12 & H2st < 4 & H2FL > 8  ~ "H2H1 Head",
                                    H5FL > 16 & H2st < 8 & H2FL > 4 ~ "H2H5 Head",
                                    H2st > 4 & H3st < 4 ~ "H2st",
                                    H2st > 4  & H3st > 4 ~"H2H3st",
                                    H2FL <= 8 & H2st <= 4 ~"H2neg",
                                    TRUE~"unclear")) %>%
                                    select("Specificity",".cell")

colnames(Specificity.P) <- c("Specificity","CELL")

Specificity <- rbind(Specificity.M,Specificity.P)

A316 <- A316 %>% left_join(Specificity, by = "CELL")

### adjust specificity so that all cells with the same clone id have the same specficity (which is the most common specificity for all cells of that clone id)
# table of clone_id vs specificity, choose the max in the row, 
adj.spec.P <- as.data.frame(A316 %>% filter(broad.celltype == "PlasmaBlast") %>% 
group_by(clone_subject_id,Specificity) %>% 
dplyr::summarise(Freq = n()) %>% 
pivot_wider(names_from = Specificity, values_from = Freq) %>% 
replace(is.na(.),0))

rownames(adj.spec.P) <- adj.spec.P$clone_subject_id
adj.spec.P <- adj.spec.P[,-1]
adj.spec.P <- adj.spec.P[,c("H2st","H2head","H2H1 Head","H2H3st","H2H5 Head","H2neg","unclear")]
adj.spec.P$adj.Specificity.P<-colnames(adj.spec.P)[apply(adj.spec.P,1,which.max)]
adj.spec.P$clone_subject_id <- rownames(adj.spec.P)
adj.spec.P <- adj.spec.P[,c("adj.Specificity.P","clone_subject_id")]
A316 <- A316 %>% left_join(adj.spec.P,by="clone_subject_id")


adj.spec.M <- as.data.frame(A316 %>% filter(broad.celltype == "MemBCell")%>% group_by(clone_subject_id,Specificity) %>%
    dplyr::summarise(Freq = n()) %>% pivot_wider(names_from = Specificity, values_from = Freq) %>% replace(is.na(.), 0))
rownames(adj.spec.M) <- adj.spec.M$clone_subject_id
adj.spec.M <- adj.spec.M[,-1]
adj.spec.M <- adj.spec.M[,c("H2st","H2head","H2H1 Head","H2H3st","H2H5 Head","H2neg","unclear")]
adj.spec.M$adj.Specificity.M<-colnames(adj.spec.M)[apply(adj.spec.M,1,which.max)]
adj.spec.M$clone_subject_id <- rownames(adj.spec.M)
adj.spec.M <- adj.spec.M[,c("adj.Specificity.M","clone_subject_id")]
A316 <- A316 %>% left_join(adj.spec.M,by="clone_subject_id")

#### For clone ids that have a mix of memory and plasmabasts chose memory adj.specificity. for clone id's that are all plasmablasts chose P.adj.spec
both <- A316 %>% group_by(clone_subject_id,broad.celltype) %>% dplyr::summarise(Freq = n()) %>% pivot_wider(names_from = broad.celltype, values_from = Freq) %>%
filter(!is.na(PlasmaBlast) & !is.na(MemBCell)) %>% pull("clone_subject_id") 

memB <- A316 %>% group_by(clone_subject_id,broad.celltype) %>% dplyr::summarise(Freq = n()) %>% pivot_wider(names_from = broad.celltype, values_from = Freq) %>%
filter(is.na(PlasmaBlast)) %>% pull("clone_subject_id")

pb <- A316 %>% group_by(clone_subject_id,broad.celltype) %>% dplyr::summarise(Freq = n()) %>% pivot_wider(names_from = broad.celltype, values_from = Freq) %>%
filter(is.na(MemBCell)) %>% pull("clone_subject_id")

A316 <- A316 %>% tidyseurat::mutate(specificity.final = case_when(clone_subject_id %in% both ~ adj.Specificity.M,
                                                                       clone_subject_id %in% memB ~ adj.Specificity.M,
                                                                       clone_subject_id %in% pb ~ adj.Specificity.P))
A316 <- A316 %>% tidyseurat::mutate(clone.cell.type = case_when(clone_subject_id %in% both ~ "both",
                                                                       clone_subject_id %in% memB ~ "memB",
                                                                       clone_subject_id %in% pb ~ "PB"))


# shared clones between two plasmablast timepoints. mark clones that exist in both plasmablasts timepoints A007PB B007PB
pb.time <- A316 %>% 
    filter(broad.celltype == "PlasmaBlast") %>%
    group_by(clone_subject_id,Timepoint) %>% 
    dplyr::summarise(Freq = n()) %>% 
    pivot_wider(names_from = Timepoint, values_from = Freq) %>%
    select(A007,B007) %>%
    mutate(clone.multi.tp.pb = case_when(is.na(A007) ~ FALSE,
                                      is.na(B007) ~ FALSE,
                                      TRUE ~ TRUE)) %>%
    select(clone_subject_id,clone.multi.tp.pb)
A316 <- A316 %>% left_join(pb.time,by="clone_subject_id")

# remove designation for membory b cells
A316$clone.multi.tp.pb[which(A316$broad.celltype == "MemBCell")] <- NA

# mark clones that exist in more than one memory bcell timepoint 
mb.time <- A316 %>% 
    filter(broad.celltype == "MemBCell") %>%
    group_by(clone_subject_id,Timepoint) %>% 
    dplyr::summarise(Freq = n()) %>% 
    pivot_wider(names_from = Timepoint, values_from = Freq) %>%
    select(clone_subject_id,B090,B014,B028,B000,B007,A000,A028,A007,B180)

mb.time$is.na <- rowSums(is.na(mb.time[,-1]))
mb.time <- mb.time %>% mutate(clone.multi.tp.mb = case_when(is.na < 8 ~ TRUE,
                                                  TRUE ~ FALSE)) %>%
select(clone_subject_id,clone.multi.tp.mb)

A316 <- A316 %>% left_join(mb.time,by="clone_subject_id")
A316$clone.multi.tp.mb[which(A316$broad.celltype == "PlasmaBlast")] <- NA


# Create column with both Sample and timepoint#
A316$Sample <- paste0(A316$Subject, "_", A316$Timepoint)
head(A316@meta.data)

saveRDS(A316, file = here::here("analysis","data_objects","04_HA","A316_specificity.rds"))
#A316 <- readRDS(file = here::here("analysis","data_objects","04_HA","A316_specificity.rds"))

## for clones that have both membcells and plasmaplasts, do they have the same specificity?
x <- A316 %>% filter(clone_subject_id %in% both) %>% distinct(clone_subject_id,.keep_all = TRUE) %>% select(c(clone_subject_id,adj.Specificity.M,adj.Specificity.P,specificity.final)) %>% filter(adj.Specificity.M != adj.Specificity.P)
y <- A316 %>% group_by(clone_subject_id,broad.celltype) %>% dplyr::summarise(Freq = n()) %>% pivot_wider(names_from = broad.celltype, values_from = Freq)
z <- left_join(x,y, by = "clone_subject_id")
write.csv(z, file = here::here("analysis","data_objects","04_HA","M_P_specificities_by_clone.csv"))

#make heatmap of m specificities vs p specifities 

#make clonal subject overlap table 
tab <- A316 %>% group_by(clone_subject_id, Subject) %>% dplyr::summarise(Freq = n()) %>% pivot_wider(names_from = Subject, values_from = Freq) 
tab$na_count <- rowSums(is.na(tab))
tab <- tab %>% filter(na_count <19)
write.csv(tab,file = here::here("analysis","data_objects","04_HA","clonal_overlap.csv"))


#### Make plots of final specificity
pdf(file = here::here("analysis","plots","04_HA","HA_scatterplots_specificity_final.pdf"))
for(i in 1:ncol(HA.combos)){
    my.plot <- A316 %>% 
      tidyseurat::join_features(features = c("H2FL", "H5FL", "NC99FL", "H3st", "MichFL", "H2st")) %>%
      pivot_wider(names_from=.feature,values_from=.abundance_Probes) %>%
      ggplot((aes_string(x = HA.combos[1,i], y = HA.combos[2,i]))) + 
        geom_point(aes(color = specificity.final), size = 0.2, ratio = 1) + 
        theme_classic() +
        scale_x_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 16, 32, 64, 128, 256, 5000)) + 
        scale_y_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 16, 32, 64, 128, 256, 5000)) +
        theme(legend.position = "right", 
            axis.title.x = element_text(size = 10), 
            axis.title.y = element_text(size = 10), 
            axis.text.x = element_text(size = 8), 
            axis.text.y = element_text(size = 8), 
            strip.text.x = element_text(size = 8)) + 
        facet_wrap(.~specificity.final) + 
        theme(aspect.ratio=1)
    plot(my.plot)
}
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

