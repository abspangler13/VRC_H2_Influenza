# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# cd /data/vrc_bsc/Abby/Experiment_316
# srun --pty bash
# module load r/4.3.0-a22xr47
# R

library(pheatmap)
library(dplyr)

### look at cluster swtiching for all H2 specific and H2 cross data 
Breaks <- seq(0, 1, length = 100)

dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster.csv"),row.names=1)
dat.s <- dat %>% filter(PGT < 0.05 & broad_specificity == "H2_Only")
x <- as.matrix(table(dat.s$FROM,dat.s$TO))
x <- rbind(x,rep(0,6))
rownames(x)[6] <- "AM2-0"
x <- rbind(x,rep(0,6))
rownames(x)[7] <- "AM1a-6"
x <- cbind(x,rep(0,7))
colnames(x)[7] <- "AM1a-4"
x <- x[c(4,5,1,7,6,2,3),c(5,6,7,1,2,3,4)]

# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"           

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at cluster swtiching for all H2 specific and H2 cross data 
Breaks <- seq(0, 1, length = 100)

dat.c <- dat %>% filter(PGT < 0.05 & broad_specificity == "H2_Cross")
x <- as.matrix(table(dat.c$FROM,dat.c$TO))
x <- rbind(x,rep(0,6))
rownames(x)[6] <- "AM2-0"
x <- rbind(x,rep(0,6))
rownames(x)[7] <- "AM1a-6"
x <- cbind(x,rep(0,7))
colnames(x)[7] <- "AM1a-4"
x <- x[c(4,5,1,7,6,2,3),c(5,6,7,1,2,3,4)]

# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"           

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_cross.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()


### now divided by specificity and exposure
Breaks <- seq(0, 1, length = 100)
# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"    

### look at H2 specific naive
dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure.csv"),row.names=1)
dat.sn <- dat %>% filter(PGT < 0.05 & specificity_exposure == "H2_Only_Naive")
x <- as.matrix(table(dat.sn$FROM,dat.sn$TO))
x <- rbind(x,rep(0,5))
rownames(x)[5] <- "AM2-0"
x <- rbind(x,rep(0,5))
rownames(x)[6] <- "AM1a-6"
x <- rbind(x,rep(0,5))
rownames(x)[7] <- "AM1a-4"
x <- cbind(x,rep(0,7))
x <- cbind(x,rep(0,7))
colnames(x)[6:7] <- c("AM1a-4","RM-2")

x <- x[c(3,4,7,6,5,1,2),c(5,7,1,2,6,3,4)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_naive.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at H2 specific exposed
dat.se <- dat %>% filter(PGT < 0.05 & specificity_exposure == "H2_Only_Exposed")
x <- as.matrix(table(dat.se$FROM,dat.se$TO))
x <- rbind(x,rep(0,4))
x <- rbind(x,rep(0,4))
rownames(x)[6:7] <- c("AM3-acute-1","AM1a-6")


x <- cbind(x,rep(0,7))
x <- cbind(x,rep(0,7))
x <- cbind(x,rep(0,7))
colnames(x)[5:7] <- c("AM1a-4","RM-2","AM3-chronic-3")

x <- x[c(4,5,1,7,2,6,3),c(4,6,5,1,2,3,7)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_exposed.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at H2 cross naive
dat.cn <- dat %>% filter(PGT < 0.05 & specificity_exposure == "H2_Cross_Naive")
x <- as.matrix(table(dat.cn$FROM,dat.cn$TO))
x <- rbind(x,rep(0,4))
x <- rbind(x,rep(0,4))
rownames(x)[6:7] <- c("AM3-acute-1","AM2-0")

# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"   

x <- cbind(x,rep(0,7))
x <- cbind(x,rep(0,7))
x <- cbind(x,rep(0,7))
colnames(x)[5:7] <- c("IgM-5","AM1a-6","AM3-chronic-3")

x <- x[c(4,5,1,2,7,6,3),c(5,4,1,6,2,3,7)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_cross_naive.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at H2 cross exposed
dat.ce <- dat %>% filter(PGT < 0.05 & specificity_exposure == "H2_Cross_Exposed")
x <- as.matrix(table(dat.ce$FROM,dat.ce$TO))
x <- rbind(x,rep(0,6))
x <- rbind(x,rep(0,6))
x <- rbind(x,rep(0,6))
rownames(x)[5:7] <- c("IgM-5","RM-2","AM1a-6")

# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"   

x <- cbind(x,rep(0,7))
colnames(x)[7] <- c("AM3-chronic-3")

x <- x[c(5,6,1,7,2,3,4),c(5,6,1,2,3,4,7)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_cross_exposed.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

#### now do everything by lineage