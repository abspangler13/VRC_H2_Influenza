## Multinomial test for distribution of clones across clusters
## Paul L. Maurizio
## 2023-10-04

library(dplyr)
library(qvalue)
## library(EMT) ## tooslow
library(ExactMultinom)
set.seed(1234)

##--------------------------------------------------------------
## START: GENERATE TOY DATA
##--------------------------------------------------------------
## number of clones
N <- 12600

## number of clusters
K <- 7

## total number of cells
ntot <- 31800

## Randomly generate cluster sums
cells.per.cluster <- rmultinom(1,ntot,rep(ntot/K,K))
rownames(cells.per.cluster) <- paste0("cluster_",c(1:K))

## Randomly generate clone sums
#cells.per.clone <- rmultinom(1,ntot,rep(ntot/N,N))
#rownames(cells.per.clone) <- paste0("clone_",c(1:N))

## Alternative for generating clone sums:
## Ideally, randomly assign cells to clone such that
## the cell count per clone ranges 1-446, mean 2.5, median 1
## For now, use a high probability for clone 1, and flat probability from 2 to N
cells.per.clone <- rmultinom(1,ntot-N,c(rep(1,300),rep(0.005,N-300))) + 1
rownames(cells.per.clone) <- paste0("clone_",c(1:N))
mean(cells.per.clone)
median(cells.per.clone)
range(cells.per.clone)

## What do the larger cluster counts look like?
head(sort(cells.per.clone, decreasing = TRUE), n=100)

## Generate toy data frame
data1 <- data.frame(matrix(NA, nrow=ntot,ncol=3))
colnames(data1) <- c("cell_name", "clone_ID", "cluster_ID")
for(i in c(1:ntot)){
  data1[i,"cell_name"] <- paste0("cell_",i)
}

## Randomly assign cluster names
cluster_ids <- NULL
for(k in c(1:K)){
  ncells <- as.vector(cells.per.cluster[rownames(cells.per.cluster)==paste0("cluster_",k),1])
  cluster_ids <- c(cluster_ids,rep(rownames(cells.per.cluster)[k],ncells))
}

## Randomly assign clone names
clone_ids <- NULL
for(n in c(1:N)){
  ncells <- as.vector(cells.per.clone[rownames(cells.per.clone)==paste0("clone_",n),1])
  clone_ids <- c(clone_ids,rep(rownames(cells.per.clone)[n],ncells))
}
clone_ids1 <- sample(clone_ids, length(clone_ids), replace=F)

## ADD ID's to data frame
data1$cluster_ID <- cluster_ids
data1$clone_ID <- clone_ids1

data1_counts <- data1 %>% count(clone_ID, cluster_ID)
data1_counts$clone_ID <- factor(data1_counts$clone_ID, levels=paste0("clone_",c(1:N)))

##-----------------------------------------------------------------
## END: TOY DATA GENERATION; Insert real data in code block above
##-----------------------------------------------------------------

## Plot relative proportion of clones per cluster
g1 <- ggplot(data1_counts, aes(x=cluster_ID, y=n, fill=clone_ID)) + 
  geom_bar(position="fill", stat="identity") + theme(legend.position="none") + coord_flip()
plot(g1)

## Plot relative cluster representation per clone, for larger cell count clones
## Select first 10 clusters

data1_counts_sub <- data1_counts[data1_counts$clone_ID %in% paste0("clone_",c(1:10)),]
g2 <- ggplot(data1_counts_sub, aes(x=clone_ID, y=n, fill=cluster_ID)) + 
  geom_bar(position="fill", stat="identity") + coord_flip()
plot(g2)

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
  thisclone <- paste0("clone_",n)
  if(keep_bool[thisclone]==TRUE){
    samplevec <- rep(0,K)
    names(samplevec) <- paste0("cluster_",c(1:K))
    data1_sub <- data1_counts[data1_counts$clone_ID %in% paste0("clone_",n),]
    for(k in 1:K){
      if(paste0("cluster_",k) %in% data1_sub$cluster_ID){
        samplevec[k] <- data1_sub[data1_sub$cluster_ID %in% paste0("cluster_",k),"n"]
      }else{samplevec[k] <- 0}
    }
    ## mtest <- multinomial.test(samplevec, prob) ## can be slow
    mtest <- multinom.test(x=samplevec, p=prob, timelimit=20, theta=1e-03) ## default is 1e-04, but using larger theta to speed it up
    ## The first p-value (e.g., pvals_ex[1]) is obtained from the probability mass, the second from Pearson's chi-square and the third from the log-likelihood ratio.
    res_list[[n]] <- mtest$pvals_ex
    pval_table[n,] <- res_list[[n]]
    cat("Finished clone", n, "!\n")
  }else{
    ##cat("Skipping clone", n, " (too few cells)!\n")
  }
}

hist(pval_table$pval.probmass, n=100)
pval_table$clone <- paste0("clone_",c(1:N))
pval_table_sub <- pval_table[pval_table$clone %in% keeps,]
#samplevec <- c(300,10,10,10,10,10,10)

## Generate permuted data sets
nperms <- 100
data1_perm_list <- list()
for(p in 1:nperms){
  data1_perm_list[[p]] <- data1
  data1_perm_list[[p]]$cluster_ID <- sample(data1_perm_list[[p]]$cluster_ID, length(data1_perm_list[[p]]$cluster_ID), replace=FALSE)
}

## Generate pvalues from permuted data sets
## For each clone, test whether the distribution differs from what is expected by chance using the multinomial test.
prob <- unname((cells.per.cluster/ntot)[,1])
res_list_perm <- list()
pval_table_perm <- list()

for(p in 1:nperms){
  cat("\n\n=========================\n")
  cat("Starting Perm",p,"!\n")
  data1_perm_counts <- data1_perm_list[[p]] %>% count(clone_ID, cluster_ID)
  data1_perm_counts$clone_ID <- factor(data1_perm_counts$clone_ID, levels=paste0("clone_",c(1:N)))

  res_list_perm[[p]] <- list()
  pval_table_perm[[p]] <- data.frame(matrix(NA, nrow=N, ncol=3))
  colnames(pval_table_perm[[p]]) <- c("pval.probmass", "pval.Chisq", "pval.LLR")
  for(n in c(1:N)){
    thisclone <- paste0("clone_",n)
    if(keep_bool[thisclone]==TRUE){
      samplevec <- rep(0,K)
      names(samplevec) <- paste0("cluster_",c(1:K))
      data1_sub <- data1_perm_counts[data1_perm_counts$clone_ID %in% paste0("clone_",n),]
      for(k in 1:K){
        if(paste0("cluster_",k) %in% data1_sub$cluster_ID){
          samplevec[k] <- data1_sub[data1_sub$cluster_ID %in% paste0("cluster_",k),"n"]
        }else{samplevec[k] <- 0}
      }
      ## mtest <- multinomial.test(samplevec, prob) ## can be slow
      mtest <- multinom.test(x=samplevec, p=prob, timelimit=20, theta=1e-03) ## default is 1e-04, but using larger theta to speed it up
      ## The first p-value (e.g., pvals_ex[1]) is obtained from the probability mass, the second from Pearson's chi-square and the third from the log-likelihood ratio.
      res_list_perm[[p]][[n]] <- mtest$pvals_ex
      pval_table_perm[[p]][n,] <- res_list_perm[[p]][[n]]
      ##cat("Perm",p,": Finished clone", n, "!\n")
    }else{
      ##cat("Perm",p,": Skipping clone", n, " (too few cells)!\n")
      }
    }
}

## View the distribution of pvalues for the permuted data sets (individually and/or pooled)
##hist(pval_table_perm[[1]]$pval.probmass, n=100)

pval_perm_pooled <- unlist(lapply(as.list(c(1:nperms)), FUN=function(x){pval_table_perm[[x]]$pval.probmass}))
pval_perm_pooled <- pval_perm_pooled[!is.na(pval_perm_pooled)]
hist(pval_perm_pooled, n=100)

## Do FDR correction based on p-values from permuted data sets
pval_table_sub$emp_p <- empPvals(stat=pval_table_sub$pval.probmass, stat0=pval_perm_pooled, pool=TRUE)
pval_table_sub$q <- qvalue(p=pval_table_sub$emp_p)$qvalue
nrow(pval_table_sub[pval_table_sub$q<0.01,])
nrow(pval_table_sub[pval_table_sub$q<0.05,])
nrow(pval_table_sub[pval_table_sub$q<0.1,])
nrow(pval_table_sub[pval_table_sub$q<0.2,])

## Pick some of the 'significant' clones (or those that are most skewed), and show their cluster distribution
## Rank the data frame based on lower empirical p values, then pick the top 10 results
pval_table_sub1 <- pval_table_sub[with(pval_table_sub, order(emp_p, decreasing=FALSE)),]
topclones <- pval_table_sub1[pval_table_sub1$pval.probmass<0.05,"clone"][1:10]
##topclones <- pval_table_sub1[pval_table_sub1$q<0.01,"clone"][1:10]
## counts for topclones:
topclones_counts <- table(data1$clone_ID)[topclones]
topclones_counts

data1_counts_top <- data1_counts[data1_counts$clone_ID %in% topclones,]
g2 <- ggplot(data1_counts_top, aes(x=clone_ID, y=n, fill=cluster_ID)) + 
  geom_bar(position="fill", stat="identity") + coord_flip()

plot(g2)

