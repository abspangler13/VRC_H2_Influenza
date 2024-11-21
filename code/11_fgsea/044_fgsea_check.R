library(tidyverse)
library(sessioninfo)

file_paths <- list.files(path = here::here("analysis","data_objects","06_DE"), pattern = "_CollapsePathways.rds", full.names = TRUE)

# Make a vector of file names
file_names <-  gsub(pattern = "_CollapsePathways.rds$", replacement = "", x = basename(file_paths))

# Read all your data into a list
data_list <- lapply(file_paths, readRDS)

# Assign file names to list elements
names(data_list) <- file_names  

total_paths <- lapply(data_list, function(x){
    length(x$parentPathways)
})

nas <- lapply(data_list, function(x){
    sum(is.na(x$parentPathways))
})

after <- lapply(data_list,function(x){
    length(x$mainPathways) + length(unique(x$parentPathways))
})

y <- unlist(total_paths)
z <- unlist(nas)
a <- unlist(after)

df <- data.frame(y,z,a)

df <- df[grep("subset",rownames(df)),]
colnames(df) <- c("total","not_collapsed","after_collapse")

write.csv(df, file = here::here("analysis","data_objects","06_DE","num_paths_collapsed.csv"))
