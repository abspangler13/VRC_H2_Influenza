## Multinomial test for distribution of clones across clusters
## Paul L. Maurizio
## 2023-10-04

library(dplyr)
library(ggplot2)
library(qvalue)
## library(EMT) ## tooslow
library(ExactMultinom)
set.seed(1234)
library(sessioninfo)

##--------------------------------------------------------------
## START: GENERATE TOY DATA
##--------------------------------------------------------------
seurat <- readRDS(file = here::here("analysis","data_objects","05_clustering","H2_Mem_G_clusters.rds"))

## number of clones
N <- length(unique(seurat$clone_subject_id))

## number of clusters
K <- length(unique(seurat$Cluster.res0.4))

## total number of cells
ntot <- dim(seurat)[2]

cells.per.cluster <- as.matrix(table(seurat@meta.data$Cluster.res0.4))

## Randomly generate clone sums
#cells.per.clone <- rmultinom(1,ntot,rep(ntot/N,N))
#rownames(cells.per.clone) <- paste0("clone_",c(1:N))

## Alternative for generating clone sums:
## Ideally, randomly assign cells to clone such that
## the cell count per clone ranges 1-446, mean 2.5, median 1
## For now, use a high probability for clone 1, and flat probability from 2 to N
cells.per.clone <- as.matrix(table(seurat@meta.data$clone_subject_id))
mean(cells.per.clone)
median(cells.per.clone)
range(cells.per.clone)

## What do the larger cluster counts look like?
head(sort(cells.per.clone, decreasing = TRUE), n=100)

## Generate toy data frame
data1 <- seurat@meta.data %>% select(barcode,clone_subject_id,Cluster.res0.4)
colnames(data1) <- c("cell_name", "clone_ID", "cluster_ID")

saveRDS(data1, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","data1.rds"))

data1_counts <- data1 %>% count(clone_ID, cluster_ID)
data1_counts$clone_ID <- factor(data1_counts$clone_ID, levels=unique(data1_counts$clone_ID))

saveRDS(data1_counts, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","data1_counts.rds"))
data1_counts <- readRDS(file = here::here("analysis","data_objects","08_cluster_clonal_overlap","data1_counts.rds"))

##-----------------------------------------------------------------
## END: TOY DATA GENERATION; Insert real data in code block above
##-----------------------------------------------------------------

## Plot relative proportion of clones per cluster
g1 <- ggplot(data1_counts, aes(x=cluster_ID, y=n, fill=clone_ID)) + 
  geom_bar(position="fill", stat="identity") + theme(legend.position="none") + coord_flip()
pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","prop_clones_per_clusters.pdf"))
plot(g1)
dev.off()


## Plot relative cluster representation per clone, for larger cell count clones
## Select first 10 clusters
clones <- unique(data1_counts$clone_ID)
clusters <- unique(data1_counts$cluster_ID)

data1_counts_sub <- data1_counts[data1_counts$clone_ID %in% clones[1:10],]
saveRDS(data1_counts_sub, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","data1_counts_sub.rds"))

g2 <- ggplot(data1_counts_sub, aes(x=clone_ID, y=n, fill=cluster_ID)) + 
  geom_bar(position="fill", stat="identity") + coord_flip()

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","prop_clones_per_clusters_sub.pdf"))
plot(g2)
dev.off()

## Decide on threshold for clone size (# cells per clone) for proceeding with testing
min_cells <- 10
keep_bool <- table(data1$clone_ID)>10
keeps <- names(table(data1$clone_ID)[keep_bool])
length(keeps)

## For each clone, test whether the distribution differs from what is expected by chance using the multinomial test.
res_list <- list()
prob <- unname((cells.per.cluster/ntot)[,1])
pval_table <- data.frame(matrix(NA, nrow=N, ncol=3))
colnames(pval_table) <- c("pval.probmass", "pval.Chisq", "pval.LLR")
for(n in c(1:N)){
  thisclone <-clones[n]
  if(keep_bool[thisclone]==TRUE){
    samplevec <- rep(0,K)
    names(samplevec) <- clusters
    data1_sub <- data1_counts[data1_counts$clone_ID %in% clones[n],]
    for(k in 1:K){
      if(clusters[k] %in% data1_sub$cluster_ID){
        samplevec[k] <- data1_sub[data1_sub$cluster_ID %in% clusters[k],"n"]
      }else{samplevec[k] <- 0}
    }
    ## mtest <- multinomial.test(samplevec, prob) ## can be slow
    Sys.time()
    mtest <- multinom.test(x=samplevec, p=prob, timelimit=240, theta=1e-04) ## default is 1e-04, but using larger theta to speed it up
    ## The first p-value (e.g., pvals_ex[1]) is obtained from the probability mass, the second from Pearson's chi-square and the third from the log-likelihood ratio.
    res_list[[n]] <- mtest$pvals_ex
    pval_table[n,] <- res_list[[n]]
    cat("Finished clone", n, "!\n")
  }
  
  else{
    ##cat("Skipping clone", n, " (too few cells)!\n")
  }
}

res_list <- Filter(Negate(is.null), res_list) 
saveRDS(res_list, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","res_list.rds"))

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","pval_hist.pdf"))
hist(pval_table$pval.probmass, n=100)
dev.off()

pval_table$clone <- clones
pval_table_sub <- pval_table[pval_table$clone %in% keeps,]

saveRDS(pval_table_sub, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","pval_table_sub.rds"))


## Generate permuted data sets
nperms <- 100
data1_perm_list <- list()
for(p in 1:nperms){
  data1_perm_list[[p]] <- data1
  data1_perm_list[[p]]$cluster_ID <- sample(data1_perm_list[[p]]$cluster_ID, length(data1_perm_list[[p]]$cluster_ID), replace=FALSE)
}

saveRDS(data1_perm_list, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","data1_perm_list.rds"))

## Generate pvalues from permuted data sets
## For each clone, test whether the distribution differs from what is expected by chance using the multinomial test.
prob <- unname((cells.per.cluster/ntot)[,1])
res_list_perm <- list()
pval_table_perm <- list()

for(p in 1:nperms){
  cat("\n\n=========================\n")
  cat("Starting Perm",p,"!\n")
  data1_perm_counts <- data1_perm_list[[p]] %>% count(clone_ID, cluster_ID)
  data1_perm_counts$clone_ID <- factor(data1_perm_counts$clone_ID, levels=unique(data1_perm_counts$clone_ID))

  res_list_perm[[p]] <- list()
  pval_table_perm[[p]] <- data.frame(matrix(NA, nrow=N, ncol=3))
  colnames(pval_table_perm[[p]]) <- c("pval.probmass", "pval.Chisq", "pval.LLR")
  for(n in c(1:N)){
    thisclone <- clones[n]
    if(keep_bool[thisclone]==TRUE){
      samplevec <- rep(0,K)
      names(samplevec) <- clusters
      data1_sub <- data1_perm_counts[data1_perm_counts$clone_ID %in% clones[n],]
      for(k in 1:K){
        if(clusters[k] %in% data1_sub$cluster_ID){
          samplevec[k] <- data1_sub[data1_sub$cluster_ID %in% clusters[k],"n"]
        }else{samplevec[k] <- 0}
      }
      ## mtest <- multinomial.test(samplevec, prob) ## can be slow
      Sys.time()
      mtest <- multinom.test(x=samplevec, p=prob, timelimit=240, theta=1e-04) ## default is 1e-04, but using larger theta to speed it up
      Sys.time()
      ## The first p-value (e.g., pvals_ex[1]) is obtained from the probability mass, the second from Pearson's chi-square and the third from the log-likelihood ratio.
      res_list_perm[[p]][[n]] <- mtest$pvals_ex
      pval_table_perm[[p]][n,] <- res_list_perm[[p]][[n]]
      ##cat("Perm",p,": Finished clone", n, "!\n")
    }else{
      ##cat("Perm",p,": Skipping clone", n, " (too few cells)!\n")
      }
    }
}

# ### remove NULL
# for(i in 1:100){
#   res_list_perm[[i]] <- Filter(Negate(is.null), res_list_perm[[i]]) 
# }

# saveRDS(res_list_perm, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","res_list_perm.rds"))

# for(i in 1:nperms){
#   pval_table_perm[[i]]$clone <- clones
#   pval_table_perm[[i]] <- pval_table_perm[[i]][pval_table_perm[[i]]$clone %in% keeps,]
# }

# saveRDS(pval_table_perm, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","pval_table_perm.rds"))
 
## View the distribution of pvalues for the permuted data sets (individually and/or pooled)
##hist(pval_table_perm[[1]]$pval.probmass, n=100)
pval_perm_pooled <- unlist(lapply(as.list(c(1:nperms)), FUN=function(x){pval_table_perm[[x]]$pval.probmass}))
pval_perm_pooled <- pval_perm_pooled[!is.na(pval_perm_pooled)]

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","perm_hist.pdf"))
hist(pval_perm_pooled, n=100)
dev.off()

## Do FDR correction based on p-values from permuted data sets
pval_table_sub$emp_p <- empPvals(stat=pval_table_sub$pval.probmass, stat0=pval_perm_pooled, pool=TRUE)
pval_table_sub$q <- qvalue(p=pval_table_sub$emp_p)$qvalue
#pval_table_sub$q <- qvalue_truncp(p=pval_table_sub$emp_p)$qvalues
nrow(pval_table_sub[pval_table_sub$q<0.01,])
nrow(pval_table_sub[pval_table_sub$q<0.05,])
nrow(pval_table_sub[pval_table_sub$q<0.1,])
nrow(pval_table_sub[pval_table_sub$q<0.2,])

# > nrow(pval_table_sub[pval_table_sub$q<0.01,])
# [1] 0
# > nrow(pval_table_sub[pval_table_sub$q<0.05,])
# [1] 0
# > nrow(pval_table_sub[pval_table_sub$q<0.1,])
# [1] 0
# > nrow(pval_table_sub[pval_table_sub$q<0.2,])
# [1] 334

## Pick some of the 'significant' clones (or those that are most skewed), and show their cluster distribution
## Rank the data frame based on lower empirical p values, then pick the top 10 results
pval_table_sub1 <- pval_table_sub[with(pval_table_sub, order(emp_p, decreasing=FALSE)),]
saveRDS(pval_table_sub1, file = here::here("analysis","data_objects","08_cluster_clonal_overlap","pval_table_sub1.rds"))
pval_table_sub1 <- readRDS(file = here::here("analysis","data_objects","08_cluster_clonal_overlap","pval_table_sub1.rds"))
 
topclones <- pval_table_sub1[pval_table_sub1$pval.probmass<0.05,"clone"][1:60]  ## why top 10 and not bottom 10?
##topclones <- pval_table_sub1[pval_table_sub1$q<0.01,"clone"][1:10]
## counts for topclones:
topclones_counts <- table(data1$clone_ID)[topclones]
topclones_counts

data1_counts_top <- data1_counts[data1_counts$clone_ID %in% topclones,]
g2 <- ggplot(data1_counts_top, aes(x=clone_ID, y=n, fill=cluster_ID)) + 
  geom_bar(position="fill", stat="identity") + coord_flip()

pdf(file = here::here("analysis","plots","08_cluster_clonal_overlap","top_clones.pdf"))
plot(g2)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()