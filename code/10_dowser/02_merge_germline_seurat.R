library(tidyseurat)
library(tidyverse)
library(vroom)

# read in new tsv files
clones_34 <- vroom(list.files(path = here::here("analysis","VDJ_run3_run4_by_subject"),pattern = "*heavy_light_clone-pass_germ-pass.tsv", full.names = TRUE)) #61492 cells, 29087 clones 
clones_1 <- vroom(list.files(path = here::here("analysis","VDJ_run1_by_subject"),pattern = "*heavy_light_clone-pass_germ-pass.tsv", full.names = TRUE)) # 21624 cells, 9028 clones
clones_2 <- vroom(list.files(path = here::here("analysis","VDJ_run2_by_subject"),pattern = "*heavy_light_clone-pass_germ-pass.tsv", full.names = TRUE)) # 3293 cells, 1562 clones
clones <- rbind(clones_1,clones_2)
clones <- rbind(clones,clones_34)

clones <- clones %>% select(junction_length,sequence_id,germline_alignment_d_mask,germline_v_call,germline_d_call,germline_j_call)

seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering", "H2_Mem_G_clusters.rds"))

#check that there aren't multiple sequences per cell since already been paired up heavy and light per cell
length(unique(clones$sequence_id))
# [1] 86409
length(unique(clones$cell_id)) 
# [1] 86409

# add these columns based on sequence_id germline_alignment_d_mask	germline_v_call	germline_d_call	germline_j_call
seurat <- left_join(seurat,clones, by = "sequence_id")

saveRDS(seurat, file = here::here("analysis","data_objects","10_dowser","seurat_germline.rds"))