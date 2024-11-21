library(tidyverse)
library(Seurat)
library(readr)
library(plyr)
library(alakazam)
library(shazam)
library(devtools)
library(sessioninfo)
library(vroom)

# helpful function combined <- AddMetaData(combined, metadata = all_cln_collapsed)
run_info <- read.csv(file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"))

####### Add VDJ information ###############
#Read in all annotated vdj from 10X and combine them#
mem.vdj.filter_34 <- vroom(paste0("analysis/VDJ_run3_run4_by_subject/",list.files(path = here::here("analysis","VDJ_run3_run4_by_subject"),pattern = "*annotations.csv")))
mem.vdj.filter_34 <- mem.vdj.filter_34[,-1]
mem.vdj.filter_1 <-vroom(paste0("analysis/VDJ_run1_by_subject/",list.files(path = here::here("analysis","VDJ_run1_by_subject"),pattern = "*annotations.csv"))[c(4,7:11)])
mem.vdj.filter_1 <- mem.vdj.filter_1[,-1]
mem.vdj.filter_2 <- vroom(paste0("analysis/VDJ_run2_by_subject/",list.files(path = here::here("analysis","VDJ_run2_by_subject"),pattern = "*annotations.csv")))
mem.vdj.filter <-rbind(mem.vdj.filter_1,mem.vdj.filter_2)
mem.vdj.filter <-rbind(mem.vdj.filter,mem.vdj.filter_34)

#[1] 201214     31

#using VDJpair function to select only barcodes with 1 HC paired with 1 LC#
files <- mem.vdj.filter
files <- list(files)
temp.file <- list()
paired.files <- list()
i = 1

  temp.file <- files[[i]] %>% dplyr::select(-is_cell, -contig_id, -high_confidence, -full_length, -productive, -raw_clonotype_id, -raw_consensus_id)
  
  temp.file.h <- temp.file %>% dplyr::filter(chain == "IGH")
  #temp.file.k <- temp.file %>% dplyr::filter(chain == "IGK")
  temp.file.l <- temp.file %>% dplyr::filter(chain %in% c("IGL","IGK"))

  #remove duplicate barcodes
  temp.file.h <- temp.file.h[!duplicated(temp.file.h$barcode),]
  #temp.file.k <- temp.file.k[!duplicated(temp.file.k$barcode),]
  temp.file.l <- temp.file.l[!duplicated(temp.file.l$barcode),]

  #temp.HK <- dplyr::left_join(temp.file.h, temp.file.k, by = "barcode")
  temp.HL <- dplyr::left_join(temp.file.h, temp.file.l, by = "barcode")
  
  my.colnames <- c("barcode","length", "chain","v_gene","d_gene",
  "j_gene","c_gene","fwr1","fwr1_nt","cdr1","cdr1_nt","fwr2",
  "fwr2_nt","cdr2","cdr2_nt","fwr3","fwr3_nt","cdr3",
   "cdr3_nt","fwr4","fwr4_nt","reads" ,"umis", "exact_subclonotype_id")
  #colnames(temp.HK) <- c(my.colnames, paste(my.colnames[-1], ".l", sep = ""))
  colnames(temp.HL) <- c(my.colnames, paste(my.colnames[-1], ".l", sep = ""))
  
  #temp.merge <- rbind(temp.HK, temp.HL)
  temp.HL <- temp.HL %>% distinct(barcode, .keep_all = TRUE)

  temp.merge <- temp.HL
  #temp.merge <- temp.merge[!(temp.merge$barcode %in% temp.merge[duplicated(temp.merge$barcode),]$barcode),]
  idents = NULL
  if(length(idents) == length(files)) {
    temp.merge$orig <- paste(idents[i])
  } else {
    temp.merge$orig <- paste("file", i, sep = "_")
  }
  
  paired.files[[i]] <- temp.merge

rtrn <- (do.call(rbind, paired.files))
All.M.vdj <- rtrn

#run_info[nrow(run_info)+1,] <- NA
run_info$vdj.ann.pair[5] <- dim(All.M.vdj)[1]

All.M.vdj$CELL <- All.M.vdj$barcode

write.csv(All.M.vdj, file = here::here("analysis","data_objects","03_vdj", "All.M.vdj_all.csv")) 

# Read in clones determined by Immcantation. Used heavy and light chain to determine clones. 
clones_34 <- vroom(paste0("analysis/VDJ_run3_run4_by_subject/",list.files(path = here::here("analysis","VDJ_run3_run4_by_subject"),pattern = "*heavy_light_clone-pass.tsv"))) #  29091 clones, 61496 cells/sequenes
clones_1 <- vroom(paste0("analysis/VDJ_run1_by_subject/",list.files(path = here::here("analysis","VDJ_run1_by_subject"),pattern = "*heavy_light_clone-pass.tsv"))) #  9030 clones, 21626 cells/sequenes
clones_2 <- vroom(paste0("analysis/VDJ_run2_by_subject/",list.files(path = here::here("analysis","VDJ_run2_by_subject"),pattern = "*heavy_light_clone-pass.tsv"))) #  1562 clones, 3293 cells/sequenes
clones <- rbind(clones_1,clones_2)
clones <- rbind(clones,clones_34)

colnames(clones)[grep("cell_id",colnames(clones))] <- "CELL"
run_info$clones[5] <- dim(clones)[1]

#Merging 10X annotation sheet with clonal assignment# left join so only keep cells that have a clone id from immcantation. get rid of remaining contigs from 10x data
All.vdj.M.clones <- left_join(clones, All.M.vdj, type = "left", by = "CELL") #each naive b-cell_id has unique vdj. daughter cells are clones

# Add in LC info from Imcantation #
#have to put in the output of light chain parse Select function

All.LC_34 <- vroom(paste0("analysis/VDJ_run3_run4_by_subject/",list.files(path = here::here("analysis","VDJ_run3_run4_by_subject"),pattern = "light_parse-select.tsv")))
All.LC_1 <-vroom(paste0("analysis/VDJ_run1_by_subject/",list.files(path = here::here("analysis","VDJ_run1_by_subject"),pattern = "light_parse-select.tsv")))
All.LC_2 <-vroom(paste0("analysis/VDJ_run2_by_subject/",list.files(path = here::here("analysis","VDJ_run2_by_subject"),pattern = "light_parse-select.tsv")))
All.LC <- rbind(All.LC_1,All.LC_2)
All.LC <- rbind(All.LC,All.LC_34)

#All L_ to all column names to distinquish from HC#
colnames(All.LC) <- paste0("LC_", colnames(All.LC))

#Add barcode column
colnames(All.LC)[grep("LC_cell_id",colnames(All.LC))] <- "CELL"

# Remove LCs where more than one per barcode (cell) #
All.LC.singlet <- All.LC[!duplicated(All.LC$CELL),]
run_info$LC.singlet[5] <- dim(All.LC.singlet)[1]

#Join LC data with all other data #
All.vdj.M <- left_join(All.vdj.M.clones, All.LC.singlet, by = "CELL")

#remove unwanted *
All.vdj.M$v_call <- gsub("[*].*$", "", All.vdj.M$v_call)
All.vdj.M$LC_v_call <- gsub("[*].*$", "", All.vdj.M$LC_v_call)

#function from shazam. calculates somatic hyper- mutation of vdj. compares sequence of vdj to germline. somatic mutations in vdj. naive bcell has undergone recombination but not somatic hypermutation
All.vdj.M <- observedMutations(All.vdj.M, sequenceColumn = "sequence_alignment", 
                                   germlineColumn = "germline_alignment", 
                                   regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = TRUE)


All.vdj.M  <- All.vdj.M %>% select(-contains("LC_d"))

#final product of vdj
write.csv(All.vdj.M, file = here::here("analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))

# check for where HC, but no LC #
apply(All.vdj.M, 2, FUN=function(x) length(which(is.na(x))))

write.csv(run_info,file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"),row.names = FALSE) 

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
