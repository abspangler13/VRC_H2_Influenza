library(Seurat)
library(tidyverse)
library(sessioninfo)
library(tidyseurat)

A316.dsb <-readRDS(file = here::here("analysis","data_objects","02_dsb_normalization","A316_dsb_all.rds"))
####################################### Look at DSB Normalized Protein Expression################################

run_info <- read.csv(file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"))

A316.dsb %>% count(Timepoint)

pdf(file = here::here("analysis","plots","02_dsb_normalization","CD14_v_CD19_all.pdf"))
A316.dsb %>%
  join_features(features = c("P-CD14","P-CD19")) %>% 
  select(one_of(c(".cell", ".feature",".abundance_Prot","Timepoint"))) %>% 
  pivot_wider(names_from = .feature, values_from = .abundance_Prot) %>% 
  tidyseurat::ggplot(aes(x= `P-CD14`, y=`P-CD19`)) + 
  geom_point(aes(color = Timepoint), size=0.2) +
  scale_x_continuous(limits = c(-2, 25)) +
  scale_y_continuous(limits = c(-2, 25)) +
  facet_wrap(.~Timepoint) +
  theme(aspect.ratio=1) +
  theme_bw()
dev.off()


#Label Bcells in Meta data
A316.dsb$Bcell <- A316.dsb %>% 
  tidyseurat::join_features(features = c("P-CD19","P-CD3","P-CD14","P-CD56")) %>% 
  pivot_wider(names_from = .feature, values_from = .abundance_Prot) %>% 
  mutate(Bcell = if_else(`P-CD19` > 2.5 & `P-CD3` < 7.5 & `P-CD14` < 5 & `P-CD56` < 10, TRUE, FALSE)) %>%
  pull("Bcell")

run_info$bcells[5] <- table(A316.dsb$Bcell)[2]

write.csv(run_info,file = here::here("analysis","data_objects","01_build_seurat","run_info_stats.csv"),row.names=FALSE)
saveRDS(A316.dsb, file = here::here("analysis","data_objects","02_dsb_normalization","A316_dsb_all.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()