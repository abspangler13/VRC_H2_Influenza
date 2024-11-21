# module load immcantation/4.4.0
# singularity shell $IMMCAN_SING
# R

library(dowser)
library(dplyr)

#load trees with meta data
trees <- readRDS(file = "analysis/data_objects/10_dowser/Trees_igphyml_Time_meta.rds")

## remove mAb clones
B180 <- read.csv(file = here::here("analysis","data_objects","10_dowser","clones.for.dowser.csv"))
B180 <- B180 %>% filter(!desription2 == "doesn't fit criteria")
trees <- trees %>% filter(clone_id %in% B180$clone_id)
saveRDS(trees, file = here::here("analysis","data_objects","10_dowser","Trees_igphyml_Time_meta_filtered.rds"))

### measuring evolution ###
# do i need to collapse trees before doing this?
evo.trees <- correlationTest(trees,time = "Timepoint.num",permutations = 10000,perm_type = "uniform") # 
evo.trees <- evo.trees[order(evo.trees$p),]

saveRDS(evo.trees,file = here::here("analysis","data_objects","10_dowser","Evo_Trees_igphyml_Time_meta_filtered.rds"))
#evo.trees <-readRDS(file = here::here("analysis","data_objects","10_dowser","Evo_Trees_igphyml_Time_meta_filtered.rds"))

dat <- select(evo.trees, clone_id,slope,correlation,p,vac_grp,Specificity,Subject,broad_specificity,Exposure)
dat <- dat %>% mutate(sig = case_when(p < 0.05 ~ TRUE,
                              p > 0.05 ~ FALSE))

write.csv(dat,file = here::here("analysis","data_objects","10_dowser","Evo_dat_Time_uniform.csv"))
#dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","Evo_dat_Time_uniform.csv"))

dat <-dat %>% filter(broad_specificity == "H2_Cross")

table(dat$sig,dat$Exposure)

pdf(file = here::here("analysis","plots","10_dowser","Evo_scatterplot_uniform.pdf"))
ggplot(dat,aes(x=vac_grp,y=slope,colour=sig,shape=Specificity)) + 
geom_point(position=position_jitter(),size = 5) +
theme_bw() +
theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=18,face="bold"))

ggplot(dat,aes(x=vac_grp,y=correlation,colour=sig,shape=Specificity)) + 
geom_point(position=position_jitter(),size = 5) +
theme_bw() +
theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=18,face="bold"))
dev.off()

# pdf(file = here::here("analysis","plots","10_dowser","test.pdf"))
# ggplot(dat,aes(x=vac_grp,y=slope,colour=sig,shape=Specificity)) + 
# geom_point(position=position_jitter(),size = 5) +
# theme_bw() +
# theme(axis.text=element_text(size=16,face="bold"),
#         axis.title=element_text(size=18,face="bold"))
# dev.off()

### Evoltuion plots
# pdf(file = here::here("analysis","plots","10_dowser",paste0("Evo_violin_",names[k],".pdf")))
# ggplot(dat, aes(x = vac_grp,y = slope, fill = vac_grp)) +
# geom_violin() +
# geom_jitter() +
# theme_bw() +
# theme(axis.text=element_text(size=16,face="bold"),
#         axis.title=element_text(size=18,face="bold"))


# ggplot(dat, aes(x = Specificity,y = slope, fill = Specificity)) +
# geom_violin() +
# geom_jitter() +
# theme_bw() +
# theme(axis.text=element_text(size=16,face="bold"),
#         axis.title=element_text(size=18,face="bold"))


# ggplot(dat, aes(x = vac_grp,y = correlation, fill = vac_grp)) +
# geom_violin() +
# geom_jitter() +
# theme_bw() +
# theme(axis.text=element_text(size=16,face="bold"),
#         axis.title=element_text(size=18,face="bold"))


# ggplot(dat, aes(x = Specificity,y = correlation, fill = Specificity)) +
# geom_violin() +
# geom_jitter() +
# theme_bw() +
# theme(axis.text=element_text(size=16,face="bold"),
#         axis.title=element_text(size=18,face="bold"))
# dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()