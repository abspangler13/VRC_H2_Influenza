# ##### https://dowser.readthedocs.io/en/stable/vignettes/Discrete-Trait-Vignette/

# cd /data/vrc_bsc/Abby/Experiment_316
# srun --pty bash
# module load r/4.3.0-a22xr47
# R

library(pheatmap)
library(dplyr)

### look at cluster swtiching for all data together
Breaks <- seq(0, 12, length = 100)
dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_by_lineage_all.csv"),row.names=1)
n <- length(unique(dat$CLONE))
# [1] 147

dat.a <- dat %>% filter(PGT < 0.05)
x <- as.matrix(table(dat.a$FROM,dat.a$TO))

#re order matrix
x <- x[c(6,7,1,2,3,4,5),c(6,7,1,2,3,4,5)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_by_lineage_all.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### divide matrix by the total number of lineages
y <- (x/n) * 100 

Breaks <- seq(0, 8, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_by_lineage_all_proportion.pdf"))
pheatmap(y,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at cluster swtiching for all H2 specific 
dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_by_lineage.csv"),row.names=1)
dat.s <- dat %>% filter(broad_specificity == "H2_Only")
n <- length(unique(dat.s$CLONE))
# [1] 101
dat.s <- dat.s %>% filter(PGT < 0.05)
x <- as.matrix(table(dat.s$FROM,dat.s$TO))

x <- x[c(6,7,1,2,3,4,5),c(6,7,1,2,3,4,5)]

# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"           
Breaks <- seq(0, 10, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### divide matrix by the total number of lineages


y <- (x/n) * 100 

Breaks <- seq(0, 10, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_by_lineage_proportion.pdf"))
pheatmap(y,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at cluster swtiching for H2 cross data 
dat.c <- dat %>% filter(broad_specificity == "H2_Cross")
n <- length(unique(dat.c$CLONE))
# [1] 46
dat.c <- dat.c %>% filter(PGT < 0.05)
x <- as.matrix(table(dat.c$FROM,dat.c$TO))
x <- rbind(x,rep(0,7))
rownames(x)[7] <- "IgM-5"

x <- x[c(7,6,1,2,3,4,5),c(6,7,1,2,3,4,5)]

# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"           
Breaks <- seq(0, 10, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_cross_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### divide matrix by the total number of lineages

y <- (x/n) * 100

Breaks <- seq(0, 10, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_cross_by_lineage_proportion.pdf"))
pheatmap(y,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()


### now divided by specificity and exposure
# [1]"IgM-5" "RM-2" "AM1a-4"        "AM1a-6"        "AM2-0"         "AM3-acute-1"  
# [5] "AM3-chronic-3"    

### look at H2 specific naive
dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_exposure_by_lineage.csv"),row.names=1)
dat.sn <- dat %>% filter(specificity_exposure == "H2_Only_Naive")
n <- length(unique(dat.sn$CLONE))
# [1] 30

dat.sn <- dat.sn %>% filter(PGT < 0.05 )
x <- as.matrix(table(dat.sn$FROM,dat.sn$TO))

x <- x[c(6,7,1,2,3,4,5),c(6,7,1,2,3,4,5)]

Breaks <- seq(0, 6, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_naive_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### divide matrix by the total number of lineages

y <- (x/n) * 100 

Breaks <- seq(0, 17, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_naive_by_lineage_proportion.pdf"))
pheatmap(y,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at H2 specific exposed
dat.se <- dat %>% filter(specificity_exposure == "H2_Only_Exposed")
n <- length(unique(dat.se$CLONE))
# [1] 71
dat.se <- dat.se %>% filter(PGT < 0.05)
x <- as.matrix(table(dat.se$FROM,dat.se$TO))

x <- x[c(6,7,1,2,3,4,5),c(6,7,1,2,3,4,5)]

Breaks <- seq(0, 6, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_exposed_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### divide matrix by the total number of lineages


y <- (x/n) * 100

Breaks <- seq(0, 17, length = 100)
pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_exposed_by_lineage_proportion.pdf"))
pheatmap(y,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at H2 cross naive
dat.cn <- dat %>% filter(PGT < 0.05 & specificity_exposure == "H2_Cross_Naive")
x <- as.matrix(table(dat.cn$FROM,dat.cn$TO))
x <- rbind(x,rep(0,6))
rownames(x)[7] <- "IgM-5"

x <- cbind(x,rep(0,7))
colnames(x)[7] <- "IgM-5"

x <- x[c(7,6,1,2,3,4,5),c(7,6,1,2,3,4,5)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_cross_naive_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### look at H2 cross exposed
dat.ce <- dat %>% filter(PGT < 0.05 & specificity_exposure == "H2_Cross_Exposed")
x <- as.matrix(table(dat.ce$FROM,dat.ce$TO))
x <- rbind(x,rep(0,6)) # x3
rownames(x)[5:7] <- c("IgM-5","AM1a-4","AM1a-6")

x <- cbind(x,rep(0,7))
colnames(x)[7] <- "IgM-5"

x <- x[c(5,4,6,7,1,2,3),c(7,6,1,2,3,4,5)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_cross_exposed_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()


#########
# make heat maps for H2 only and vaccine group
#########
Breaks <- seq(0, 5, length = 100)

dat <- read.csv(file = here::here("analysis","data_objects","10_dowser","switches_sp_among_cluster_specificity_vac_grp_by_lineage.csv"),row.names=1)

dat.s3a <- dat %>% filter(PGT < 0.05 & specificity_vac_grp == "H2_Only_3A")
x <- as.matrix(table(dat.s3a$FROM,dat.s3a$TO))
x <- rbind(x,rep(0,6))
rownames(x)[7] <- "IgM-5"

x <- cbind(x,rep(0,7))
colnames(x)[7] <- "IgM-5"

x <- x[c(7,6,1,2,3,4,5),c(7,6,1,2,3,4,5)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_3A_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### 3B
dat.s3b <- dat %>% filter(PGT < 0.05 & specificity_vac_grp == "H2_Only_3B")
x <- as.matrix(table(dat.s3b$FROM,dat.s3b$TO))
x <- rbind(x,rep(0,5))
x <- rbind(x,rep(0,5))
rownames(x)[6:7] <- c("IgM-5","RM-2")

x <- cbind(x,rep(0,7))
x <- cbind(x,rep(0,7))
colnames(x)[6:7] <- c("IgM-5","AM1a-6")

x <- x[c(6,7,1,2,3,4,5),c(6,5,1,7,2,3,4)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_3B_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### 4A
dat.s4a <- dat %>% filter(PGT < 0.05 & specificity_vac_grp == "H2_Only_4A")
x <- as.matrix(table(dat.s4a$FROM,dat.s4a$TO))
x <- rbind(x,rep(0,1)) #X3

rownames(x)[5:7] <- c("AM1a-6","AM3-chronic-3","AM2-0")

x <- cbind(x,rep(0,7)) #3

colnames(x)[5:7] <- c("AM1a-4", "AM2-0","AM3-acute-1")

x <- x[c(3,4,1,5,7,2,6),c(3,4,5,1,6,7,2)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_4A_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()

### 4B
dat.s4b <- dat %>% filter(PGT < 0.05 & specificity_vac_grp == "H2_Only_4B")
x <- as.matrix(table(dat.s4b$FROM,dat.s4b$TO))

x <- x[c(6,7,1,2,3,4,5),c(6,7,1,2,3,4,5)]

pdf(file = here::here("analysis","plots","10_dowser","switches_among_cluster_H2_specific_4B_by_lineage.pdf"))
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 20, breaks=Breaks,display_numbers = T,number_format = "%.0f")
dev.off()