# module load immcantation/4.4.0
# singularity shell $IMMCAN_SING
# R

# Load required packages
library(alakazam)
library(dowser)
library(dplyr)
library(ggtree)
library(readr)

# load example AIRR tsv data
my.data <- readRDS(file = here::here("analysis","data_objects","10_dowser","seurat_germline.rds"))
my.data <- my.data@meta.data

#"B007" "B014" "B028" "B090" "B000" "B180" "A000" "A028" "A007"
my.data <- my.data %>% mutate(Timepoint.num = case_when(Timepoint == "A000" ~ 0,
                                                        Timepoint == "A007" ~ 7,
                                                        Timepoint == "A028" ~ 28,
                                                        Timepoint == "B000" ~ 140,
                                                        Timepoint == "B007" ~ 147,
                                                        Timepoint == "B014" ~ 154,
                                                        Timepoint == "B028" ~ 168,
                                                        Timepoint == "B090" ~ 202,
                                                        Timepoint == "B180" ~ 292))
my.data <- my.data %>% mutate(Timepoint.wk = case_when(Timepoint == "A000" ~ "0",
                                                        Timepoint == "A007" ~ "1",
                                                        Timepoint == "A028" ~ "4",
                                                        Timepoint == "B000" ~ "16",
                                                        Timepoint == "B007" ~ "17",
                                                        Timepoint == "B014" ~ "18",
                                                        Timepoint == "B028" ~ "20",
                                                        Timepoint == "B090" ~ "28",
                                                        Timepoint == "B180" ~ "40"))

my.data <- my.data %>% mutate(exposure = case_when(vac_grp == "3A" ~ "Naive",
                                                        vac_grp == "3B" ~ "Exposed",
                                                        vac_grp == "4A" ~ "Naive",
                                                        vac_grp == "4B" ~ "Exposed"))


B180 <- read.csv(file = here::here("analysis","data_objects","10_dowser","clones.for.dowser.csv")) #25 cells, cell at last timepoint

# my.data has 12,648 clones
my.data <- my.data %>% filter(clone_subject_id %in% B180$clone_id)
# my.data has 168 clones

table(my.data$c_call, useNA = "always")
# IGHA1 IGHA2  IGHD IGHG1 IGHG2 IGHG3 IGHG4  IGHM  <NA> 
#  1659     3    54  4263   295   603    30  1127   147 


##### remove cells that don't have c_call
my.data <- my.data %>% filter(!is.na(c_call))

length(unique(my.data$clone_id))
# still have 168 clones

# ### remove cells that are H2_cross
# my.data <- my.data %>% filter(!Spec.Broad == "H2_Cross")

# dim(my.data)
# #[1] 5689   92
# ## down to 121 clones after removing H2_cross clones

x <- as.data.frame(table(my.data$clone_subject_id))
x <- x[order(x$Freq),]
write.csv(x,file = here::here("analysis","data_objects","10_dowser","cells_per_dowser_clone.csv"))

### not necessary
# x.clones <- x %>% select(Var1)
# my.data <- my.data[my.data$clone_subject_id %in% x.clones$Var1,]
#### 

# Process example data into proper format, store isotype (optional)
t1 = c("Timepoint.num") 
t2 = c("Cluster.res0.4") 
t3 = c("c_call") 
t4 = c("Timepoint.num","Cluster.res0.4")
t5 = c("Timepoint.num",'c_call') 
t6 = c("Timepoint.num","Cluster.res0.4","c_call","specificity.SFA")
traits <- list(t1,t2,t3,t4,t5,t6)
names <- c("Time","Cluster","c_call","Time_Cluster","Time_ccall","All")

#k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k = 3
my.clones = formatClones(my.data, clone = "clone_subject_id",traits =traits[[k]],text_fields = c("vac_grp", "Subject","Spec.Broad","Timepoint","Timepoint.wk","exposure")) # include exposure and weeks. 
saveRDS(my.clones,file = here::here("analysis","data_objects","10_dowser",paste0("my_clones_",names[k],".rds")))
#my.clones <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("my_clones_",names[k],".rds")))

## play with using different traits. can use more than one at a time. collapses identical sequnces together unless they differ by a trait value, then they will  will be kept separate.
## numerical fields will be added together if their sequences are to be collapsed. 

# Build maxmimum parsimony trees for first two clones using 
# phangorn package in R
### account for uncertainty in tree topology by setting fixtrees=FALSE
trees <- getTrees(my.clones, build = "igphyml", exec="/usr/local/share/igphyml/src/igphyml", nproc=1)
saveRDS(trees,file = here::here("analysis","data_objects","10_dowser",paste0("Trees_igphyml_",names[k],".rds")))
#trees <- readRDS(file = here::here("analysis","data_objects","10_dowser",paste0("Trees_igphyml_",names[k],".rds")))


###add meta data to trees
my.data <- my.data %>% select(clone_subject_id,vac_grp,specificity.SFA,Subject,Spec.Broad,exposure) %>% distinct(clone_subject_id, .keep_all = TRUE)
colnames(my.data) <- c("clone_id","vac_grp","Specificity","Subject","broad_specificity", "Exposure")
trees <- left_join(trees,my.data, by = "clone_id")

saveRDS(trees,file = here::here("analysis","data_objects","10_dowser",paste0("Trees_igphyml_",names[k],"_meta.rds")))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()