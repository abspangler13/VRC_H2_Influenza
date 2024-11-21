library(tidyverse)
library(Seurat)
library(tidyseurat)
library(sessioninfo)

####script to add vdj info to seurat object 
run_info <- read.csv(file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"))

#load seurat object 
A316 <- readRDS(file = here::here("analysis","data_objects","02_dsb_normalization","A316_dsb_all.rds"))

#load final product of vdj
vdj <- read.csv(file = here::here("analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))

# remove duplicated columns #
vdj <- vdj %>% distinct(CELL,.keep_all=TRUE)
rownames(vdj) <- vdj$CELL

# see number of cells with no Pr or no HT #
apply(A316@meta.data, 2, FUN=function(x) length(which(is.na(x))))
apply(vdj, 2, FUN=function(x) length(which(is.na(x))))

#remove unwanted columns from vdj meta data
vdj <- vdj[,-which(colnames(vdj)%in%c("rev_comp","productive",
                                                        "v_cigar",
                                                        "d_cigar",
                                                        "j_cigar",
                                                        "junction_length",
                                                        "np1_length",
                                                        "np2_length",
                                                        "v_sequence_start",
                                                        "v_sequence_end",
                                                        "v_germline_start",
                                                        "v_germline_end",
                                                        "d_sequence_start",
                                                        "d_sequence_end",
                                                        "d_germline_start",
                                                        "d_germline_end",
                                                        "j_sequence_start",
                                                        "j_sequence_end",
                                                        "j_germline_start",
                                                        "j_germline_end",
                                                        "chain",
                                                        "v_gene",
                                                        "d_gene",
                                                        "j_gene",
                                                        "c_gene",
                                                        "fwr1",
                                                        "fwr1_nt",
                                                        "cdr1_nt",
                                                        "fwr2",
                                                        "fwr2_nt",
                                                        "cdr2_nt",
                                                        "fwr3",
                                                        "fwr3_nt",
                                                        "cdr3",
                                                        "cdr3_nt",
                                                        "fwr4",
                                                        "fwr4_nt",
                                                        "exact_subclonotype_id",
                                                        "length.l",
                                                        "chain.l",
                                                        "v_gene.l",
                                                        "d_gene.l",
                                                        "j_gene.l",
                                                        "c_gene.l",
                                                        "fwr1.l",
                                                        "fwr1_nt.l",
                                                        "cdr1.l",
                                                        "cdr1_nt.l",
                                                        "fwr2.l",
                                                        "fwr2_nt.l",
                                                        "cdr2.l",
                                                        "cdr2_nt.l",
                                                        "fwr3.l",
                                                        "fwr3_nt.l",
                                                        "cdr3.l",
                                                        "cdr3_nt.l",
                                                        "fwr4.l",
                                                        "fwr4_nt.l",
                                                        "exact_subclonotype_id.l",
                                                        "LC_v_cigar",
                                                        "LC_j_cigar",
                                                        "LC_junction_length",
                                                        "LC_np1_length",
                                                        "LC_v_sequence_start",
                                                        "LC_v_sequence_end",
                                                        "LC_v_germline_end",
                                                        "LC_j_sequence_start",
                                                        "LC_j_sequence_end",
                                                        "LC_j_germline_start",
                                                        "LC_j_germline_end"))]



run_info$A316[5] <- dim(A316)[2]

run_info$all.vdj[5] <- dim(vdj)[1]


#add vdj meta data into seurat object 
A316@meta.data$barcode <- rownames(A316@meta.data)
colnames(A316@meta.data)[grep("barcode",colnames(A316@meta.data))] <- "CELL"
run_info$wt.vdj[5] <- length(intersect(A316@meta.data$CELL,vdj$CELL))
A316 <-left_join(A316,vdj, by="CELL")

apply(A316@meta.data, 2, FUN=function(x) length(which(is.na(x))))

#save whole seurat object
saveRDS(A316, file=here::here("analysis","data_objects","03_vdj","A316_final_all.rds"))

#remove cells that don't have vdj and save it as a separate object.
A316.vdj <- A316 %>% filter(!is.na(sequence_id))

#concatenate clone_id and subject because some clone_ids are repeated for different subjects
A316.vdj <- A316.vdj %>% tidyseurat::mutate(clone_subject_id = paste(Subject, clone_id,sep = '_'))

saveRDS(A316.vdj, file=here::here("analysis","data_objects","03_vdj","A316_final_vdj_all.rds"))
#A316.vdj <- readRDS(file=here::here("analysis","data_objects","03_vdj","A316_final_vdj_all.rds"))

write.csv(run_info,file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"), row.names=FALSE)

# Make csv with HA and VDJ data for Sarah
dat <- A316.vdj %>% tidyseurat::join_features(features=rownames(A316.vdj@assays$Probes)) %>% 
    pivot_wider(names_from=.feature,values_from=.abundance_Probes)
dat <- as.data.frame(dat)

write.csv(dat,file = here::here("analysis","data_objects","03_vdj","VDJ_HA_metadata.csv"), row.names=FALSE)

y<-A316.vdj %>% filter(Subject == #phi, orig.ident %in% c(#phi),run == "run4")
table(y$Subject,y$Timepoint)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


