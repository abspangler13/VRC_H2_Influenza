library(here)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(sessioninfo)
library(tidyseurat)
library(patchwork)
library(scuttle)

## Define some info for the samples
sample_info <- data.frame(
  sample_id = #code removed. contains PHI.,
  cellranger_id = #code removed. contains PHI.,
vac_grp = c("3B",
            "3B",
            "3B",
            "3A",
            "3A",
            "3B",
            "3B",
            "3A",
            "3A",
            "NA",
            "3A",
            "3B",
            "3B",
            "3A",
            "3A")
)

sample_info$cellranger_raw_dir <- file.path(here::here("CellRanger_run5"),
                                        sample_info$cellranger_id,"outs",
                                        "raw_feature_bc_matrix")

sample_info$cellranger_filtered_dir <- file.path(here::here("CellRanger_run5"),
                                        sample_info$cellranger_id,"outs",
                                        "filtered_feature_bc_matrix")

# sample_info$protein_csv_file = file.path(here::here("analysis","data_objects","01_build_seurat"),
#                              sample_info$sample_id,
#                              paste0(sample_info$sample_id,"_sample_specific_oligos.csv"))

sample_info$log.txt = file.path(here::here("analysis","data_objects","01_build_seurat"),
                                sample_info$sample_id,
                                paste0(sample_info$sample_id,"_log.txt"))

sample_info$run = "run5"

for (i in 1:nrow(sample_info)){
  dir.create(file.path(here::here("analysis","data_objects","01_build_seurat"),
                       sample_info$sample_id[i]))
  dir.create(file.path(here::here("analysis","plots","01_build_seurat"),
                       sample_info$sample_id[i]))
  file.create(file.path(here::here("analysis","data_objects","01_build_seurat"),
               sample_info$sample_id,
               paste0(sample_info$sample_id,"_log.txt")))
}

#stopifnot(all(file.exists(sample_info$protein_csv_file)))

QC.plots.b4 <- list()
QC.plots.after <- list()
background.plots <- list()
HTO.vln.plots <- list()
HTO.ridge.plots <- list()
HTO.scatter.plots <- list()
RNA.vln.plots <- list()
RNA.count.vlnplot <- list()

for (i in 1:nrow(sample_info)){
  message("Starting sample ",sample_info$sample_id[i])
  # Read in raw matrix file only so we can get the barcodes that are empty and use them to make a negative seurat object for dsb normalization
  raw.data <- Read10X(data.dir = file.path(sample_info$cellranger_raw_dir[i]))
  message("genes = ", dim(raw.data$`Gene Expression`)[1],", cell barcodes = ",dim(raw.data$`Gene Expression`)[2], " in raw matrix.")
  sample_info$raw.barcodes[i] <- dim(raw.data$`Gene Expression`)[2]

  # Read in cell ranger output. creates a list of two matrices. one is feature barcode 
  # for gene expression and the other is for the CSO custom assay
  filtered.data <- Read10X(data.dir = file.path(sample_info$cellranger_filtered_dir[i]))
  message("genes = ", dim(filtered.data$`Gene Expression`)[1],", cell barcodes = ",dim(filtered.data$`Gene Expression`)[2], " in filtered matrix.")
  sample_info$filtered.barcodes[i] <- dim(filtered.data$`Gene Expression`)[2]

  #get cell barcodes for cells and empty drops
  stained_cells = colnames(filtered.data$`Gene Expression`)
  background = setdiff(colnames(raw.data$`Gene Expression`), stained_cells)
  message(length(background), " empty droplets according to CellRanger")
  sample_info$CR.empty.drops[i] <- length(background)
  
  # Create Seurat Object using whole transcriptome gene expression matrix #
  seurat <- CreateSeuratObject(counts = raw.data$`Gene Expression`, project = sample_info$sample_id[i])

  # add indicator for barcodes Cell Ranger called as cells
  seurat@meta.data$CR.drop.class = ifelse(rownames(seurat@meta.data) %in% stained_cells, 'cell', 'background')
  
  # CITE SEQ/Probes and HTOs must be added in separately (as assays) #
  # look st list of oligos to figure out which to include #
  oligoslist <- rownames(as.matrix(raw.data$`Antibody Capture`)) 

  #make csv of oligos
  write.csv(oligoslist, file = here::here("analysis","data_objects","01_build_seurat",sample_info$sample_id[i],"Oligos_list.csv"))

  #load oligo csv that lists which hashtags are present in each sample
  my.csv <- read.csv(here::here("analysis","data_objects","01_build_seurat","Run5_sample_hashtags.csv"))
  hto.oligos <- my.csv$Oligo[which(my.csv$sample == sample_info$sample_id[i] & my.csv$in.sample == TRUE)]
  
  probes.oligos <- oligoslist[11:16]
  
  prots.oligos <- oligoslist[17:77]
  
  # Divide up custom assay matrix so we can add hashtags, probes, and proteins as separate assays to the seurat object
  # - excluded some oligos that didn't actually use.
  HTO <- as.matrix(raw.data$`Antibody Capture`[hto.oligos,])
  message("HTO assay dim ",dim(HTO)[1]," ",dim(HTO)[2] )
  
  Probes <- as.matrix(raw.data$`Antibody Capture`[probes.oligos,])
  message("Probes assay dim ",dim(Probes)[1]," ",dim(Probes)[2] )
  
  Prot <- as.matrix(raw.data$`Antibody Capture`[prots.oligos,])
  
  #### remove CD-20
  Prot <- Prot[!(row.names(Prot)%in% "P_CD20"),]
  message("Prot assay dim ",dim(Prot)[1]," ",dim(Prot)[2])
  
  #Add matrices as separate assays to seurat object 
  seurat[["HTO"]] <- CreateAssayObject(counts = HTO)
  seurat[["Prot"]] <- CreateAssayObject(counts = Prot)
  seurat[["Probes"]] <- CreateAssayObject(counts = Probes)

  #log normalize for visualization purposes only
  seurat@meta.data$log_nCount_RNA <- log10(Matrix::colSums(seurat@assays$RNA@counts))
  seurat@meta.data$log_nCount_Prot <- log10(Matrix::colSums(seurat@assays$Prot@counts))

  # remove barcodes with no evidence of capture in the experiment. Need to do this so HTO demux runs
  barcodes_keep <- rownames(seurat@meta.data)[which(seurat@meta.data$log_nCount_RNA >0 & seurat@meta.data$log_nCount_Prot > 0)]
  sample_info$zero.capture[i] <- sample_info$raw.barcodes[i] - length(barcodes_keep)

  #add some meta data
  seurat$vac_grp <- sample_info$vac_grp[i]
  seurat$run <- sample_info$run[i]

  #save raw object before QC
  seurat_raw<-seurat

  # down select so that QC plots look better 
  # 
  # do another of sliding first then doing drop cells to see if I can recreate the numbers from before
  seurat.old <- seurat %>% arrange(desc(nCount_RNA))  %>% dplyr::slice(1:20000)
  source(file = here::here("analysis","code","QCGenes_FUN.R"))
  seurat.old <- QCGenes(seurat.old)[[1]]
  source(file = here::here("analysis","code","DropCells_FUN.R"))
  seurat_pos.old <- DropCells(seurat.old, min.features = 100, min.hk = 50, max.mito = 15, min.count = 100, max.count = 50000, path.out = here::here("analysis","data_objects","01_build_seurat",sample_info$sample_id[i],"QC_stats.csv"))
  sample_info$after.old.DropCells[i] <- dim(seurat_pos.old)[2]


  seurat <- subset(seurat, subset = nCount_RNA > 100 & nCount_RNA < 50000)
  message("NumCells after dropping cells with less than 50 and more the 5000 UMIs = ", dim(seurat)[2])
  sample_info$UMI_drop_100_5000[i] <- dim(seurat)[2]

  #Do QC on transcriptome and get rid of bad cells
  #used QCGenes function to calculate percent mito and hk genes. 
  source(file = here::here("analysis","code","QCGenes_FUN.R"))
  seurat<- QCGenes(seurat)[[1]] #adds columns "percent.mito" and "n.exphkgenes" to meta.data

  sample_info$avg.mito[i] <- mean(seurat@meta.data$percent.mito)
  sample_info$avg.hk[i] <- mean(seurat@meta.data$n.exp.hkgenes)
  sample_info$avg.rna[i] <- mean(seurat@meta.data$nCount_RNA)
  sample_info$avg.feature.rna[i] <- mean(seurat@meta.data$nFeature_RNA)
  sample_info$avg.hto[i] <- mean(seurat@meta.data$nCount_HTO)
  sample_info$avg.prot[i] <- mean(seurat@meta.data$nCount_Prot)
  sample_info$avg.probes[i] <- mean(seurat@meta.data$nCount_Probes)
  sample_info$feature.drop[i] <- dim(seurat %>% filter(nFeature_RNA < 100))[2]
  sample_info$hk.drop[i] <- dim(seurat %>% filter(n.exp.hkgenes < 50))[2]
  sample_info$mito.drop[i] <- dim(seurat %>% filter(percent.mito > 15))[2]

  #Change this to just use the filter function from tidyseurat. Don't need dropcells function. 
  source(file = here::here("analysis","code","DropCells_FUN.R"))
  seurat_pos <- DropCells(seurat, min.features = 100, min.hk = 50, max.mito = 15, min.count = 100, max.count = 50000, 
    path.out = here::here("analysis","data_objects","01_build_seurat",sample_info$sample_id[i],"QC_stats.csv"))

  sample_info$after.DropCells[i] <- dim(seurat_pos)[2]
  sample_info$percent.DropCells[i] <- ((dim(seurat)[2]- dim(seurat_pos)[2]) / dim(seurat)[2])*100

  seurat@meta.data$cell_barcode <- colnames(seurat)

  ##use scran too in order to compare
  sce <- seurat@assays$RNA@counts
  is.mito <- grep("MT-",rownames(sce))
  names(is.mito) <- rownames(sce)[is.mito]
  
  df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  df$cell_barcode <- rownames(df)

  sample_info$med_umi[i] <- summary(df$sum)[3]
  sample_info$mean_umi[i] <- summary(df$sum)[4]
  sample_info$max_umi[i] <- summary(df$sum)[6]
  sample_info$med_feature[i] <- summary(df$detected)[3]
  sample_info$mean_feature[i] <- summary(df$detected)[4]
  sample_info$max_feature[i] <- summary(df$detected)[6]
  sample_info$med_mito[i] <- summary(df$subsets_Mito_percent)[3]
  sample_info$mean_mito[i] <- summary(df$subsets_Mito_percent)[4]
  sample_info$max_mito[i] <- summary(df$subsets_Mito_percent)[6]

  reasons <- perCellQCFilters(df, 
    sub.fields=c("subsets_Mito_percent"))

colSums(as.matrix(reasons))

  sample_info$low_lib_size[i] <- colSums(as.matrix(reasons))[1]
  sample_info$low_n_features[i] <- colSums(as.matrix(reasons))[2]
  sample_info$high_subsets_Mito_percent[i] <- colSums(as.matrix(reasons))[3]
  sample_info$discard[i] <- colSums(as.matrix(reasons))[4]

  scran.dat <- cbind(df,reasons)
  seurat <- seurat %>% left_join(scran.dat,by = "cell_barcode", copy = TRUE)

  seurat_scran <- seurat %>% filter(discard == "FALSE")
  sample_info$after.ScranDrop[i] <- dim(seurat_scran)[2]
  sample_info$percent.ScranDrop[i] <- ((dim(seurat)[2]- dim(seurat_scran)[2]) / dim(seurat)[2])*100


  #OC on cell barcodes and drop bad cells 
  source(file = here::here("analysis","code","QCplots_FUN.R"))
  pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("QCplots_before_DropCells",sample_info$sample_id[i],".pdf")))
  QCplots(seurat)
  dev.off()

  QC.plots.b4[[i]] <- QCplots(seurat)
  
  # make note of which cells are dropped by drop cells and save in raw seurat object
  seurat_raw@meta.data$DropCell = ifelse(rownames(seurat_raw@meta.data) %in% rownames(seurat_pos@meta.data), FALSE, TRUE)
  save(seurat_raw, file=here::here("analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0("seurat_raw_",sample_info$sample_id[i],".rds")))

  #look at background drops vs cells
  pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("Background_v_cells",sample_info$sample_id[i],".pdf")))
  background.plots[[i]] <- FeatureScatter(seurat_raw,feature1 = "log_nCount_Prot",feature2 = "log_nCount_RNA",group.by = "CR.drop.class") + labs(title = sample_info$sample_id[i])
  print(background.plots[[i]])
  dev.off()

  #QC on background drops. change this to MAD approach from dsb normalization script. 
  seurat_neg <- seurat_raw %>% tidyseurat::filter(DropCell==TRUE & CR.drop.class=="background" & log_nCount_Prot>1.5 & log_nCount_Prot <3 & log_nCount_RNA <2.5)
  sample_info$seurat.neg[i] <- dim(seurat_neg)[2]
  
  seurat_neg@meta.data$barcode <- colnames(seurat_neg)
  saveRDS(seurat_neg, file = here::here("analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0(sample_info$sample_id[i],"_neg.rds")))

      # Normalizing hashtag data by centered log ratio transformation # 
    seurat_pos <- NormalizeData(seurat_pos, assay = "HTO", normalization.method = "CLR")
    # Demultiplex HT with MULTIseqDemux For each cell it's designating it with one of the hashtags or as a doublet or negative
    seurat_pos <- MULTIseqDemux(seurat_pos, assay = "HTO", quantile = 0.4) # , autoThresh = TRUE
  #}

  #if negative and above threshold for 257 call it 257 
  if(sample_info$sample_id[i] == #code removed. contains PHI. ){
    hto.dat <- as.data.frame(as.matrix(t(seurat_pos@assays$HTO@data)))
    hto.dat$barcode <- rownames(hto.dat)
    meta.dat <- seurat_pos@meta.data
    meta.dat$barcode <- rownames(meta.dat)
    dat <- left_join(meta.dat,hto.dat, by = "barcode")
    dat <- dat %>% filter(`HTO-0257` > 0.5) %>% filter(MULTI_ID == "Negative")
    seurat_pos$MULTI_ID <- as.character(seurat_pos$MULTI_ID)
    seurat_pos$MULTI_ID[which(rownames(seurat_pos@meta.data) %in% dat$barcode )] <- "HTO-0257"
    seurat_pos$MULTI_ID <- as.factor(seurat_pos$MULTI_ID)
  }

  #not sure if this line was necessary
  Idents(seurat_pos) <- "MULTI_ID"

  # Adding barcodes to meta data as separate column apart from rownames (sommetimes rownames in metadata are converted to numbers) #
  seurat_pos@meta.data$barcode <- colnames(seurat_pos)
  
  #make some plots. Go back and make plos with and with negatives and doublets so we see what dropping them looks like. 
  message("Making some plots")
  pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("VlnPlot_nCountHTO_",sample_info$sample_id[i],".pdf")))
  HTO.vln.plots[[i]] <- VlnPlot(seurat_pos, features = "nCount_HTO", group.by = "MULTI_ID", pt.size = 0.01, log = TRUE) + labs(title = sample_info$sample_id[i])
  print(HTO.vln.plots[[i]])
  dev.off()
  
  pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("RidgePlot_nCountHTO_",sample_info$sample_id[i],".pdf")))
  HTO.ridge.plots[[i]] <- RidgePlot(seurat_pos[,!(seurat_pos@meta.data$MULTI_ID %in% c("Doublet"))], features = rownames(seurat_pos@assays$HTO@counts), group.by = "MULTI_ID", ncol = 2) + labs(title = sample_info$sample_id[i])
  print(HTO.ridge.plots[[i]])
  dev.off()
  
  htos<-rownames(seurat_pos@assays$HTO)
  htos <- combn(unique(htos),2)
  pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("FeatureScatter_",sample_info$sample_id[i],".pdf")))
  for(j in 1:ncol(htos)){
    my.plot <- FeatureScatter(seurat_pos, feature1 = htos[1,j], feature2 = htos[2,j]) + labs(title = sample_info$sample_id[i]) ##chose two that are in all samples (1 and 7)
    print(my.plot)
  }
  dev.off()

  pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("VlnPlot_nCount_RNA_",sample_info$sample_id[i],".pdf")))
  RNA.vln.plots[[i]] <- VlnPlot(seurat_pos, features = "nCount_RNA", pt.size = 0.1, log = TRUE) + labs(title = sample_info$sample_id[i])
  print(RNA.vln.plots[[i]])
  dev.off()
  
  message("Frequency of each hashtag identity including doublets and negatives\n",HTO.table <- plyr::count(seurat_pos@meta.data$MULTI_classification))
  write.csv(HTO.table, file = here::here("analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0(sample_info$sample_id[i],"_HTO_counts_preQC.csv")))
  
  #add doublet and negative counts to sample_info
  multi.id.tab <- plyr::count(seurat_pos@meta.data$MULTI_ID)

  # sample_info$doublets <- NA
  if(length(grep("Doublet",multi.id.tab$x)) > 0){
    sample_info$doublets[i] <- multi.id.tab[grep("Doublet",multi.id.tab$x),2]
  } else {
    sample_info$doublets[i] <- 0
  }

  if(length(grep("Negative",multi.id.tab$x)) > 0){
    sample_info$negatives[i] <- multi.id.tab[grep("Negative",multi.id.tab$x),2]
  } else {
    sample_info$negatives[i] <- 0
  }
  # Subset only singlets. get rid of doublets and negatives
  seurat_pos <- seurat_pos %>% filter(!MULTI_ID %in% c("Negative","Doublet"))
  # seurat_pos <- subset(seurat_pos, idents %in% c("Negative", "Doublet"), invert = TRUE)
  sample_info$num.final.cells[i] <- dim(seurat_pos)[2]

  RNA.count.vlnplot[[i]] <- (VlnPlot(seurat_pos, features = "nCount_RNA",pt.size = 0) - VlnPlot(seurat_scran, features = "nCount_RNA",pt.size = 0)) / VlnPlot(seurat, features = "nCount_RNA",pt.size = 0,group.by = "orig.ident" )
  
  pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("QCplots_after_DropCells",sample_info$sample_id[i],".pdf")))
  QCplots(seurat_pos)
  QCplots(seurat_scran)
  dev.off()

  QC.plots.after[[i]] <- QCplots(seurat_pos)
  #save both postive and neative seurat objects
  saveRDS(seurat_pos, file = here::here("analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0(sample_info$sample_id[i],"_pos.rds")))
}

#save sample info table 
write.csv(sample_info, file = here::here("analysis","data_objects","01_build_seurat","sample_info_stats_run5.csv"))

#save grob lists 
# saveRDS(QC.plots.b4, file = here::here("analysis","plots","01_build_seurat","qc_b4_grobs_run4.rds")) 
# saveRDS(QC.plots.after, file = here::here("analysis","plots","01_build_seurat","qc_after_grobs_run4.rds"))
# saveRDS(background.plots, file = here::here("analysis","plots","01_build_seurat","background_grobs_run4.rds"))
# saveRDS(HTO.vln.plots,file = here::here("analysis","plots","01_build_seurat","HTO_vln_grobs_run4.rds"))
# saveRDS(HTO.ridge.plots, file = here::here("analysis","plots","01_build_seurat","HTO_ridge_grobs_run4.rds"))
# saveRDS(HTO.scatter.plots, file = here::here("analysis","plots","01_build_seurat","HTO_scatter_grobs_run4.rds"))
# saveRDS(RNA.vln.plots, file = here::here("analysis","plots","01_build_seurat","rna_vln_grobs_run4.rds"))

#plot grob lists
pdf(file = here::here("analysis","plots","01_build_seurat","qc_b4_grobs_run5.pdf"))
lapply(QC.plots.b4,plot)
dev.off()

pdf(file = here::here("analysis","plots","01_build_seurat","qc_after_grobs_run5.pdf"))
lapply(QC.plots.after,plot)
dev.off()

pdf(file = here::here("analysis","plots","01_build_seurat","background_grobs_run5.pdf"))
marrangeGrob(background.plots, ncol = 2, nrow = 2)
dev.off()

pdf(file = here::here("analysis","plots","01_build_seurat","HTO_vln_grobs_run5.pdf"))
lapply(HTO.vln.plots, plot)
dev.off()

pdf(file = here::here("analysis","plots","01_build_seurat","HTO_ridge_grobs_run5.pdf"))
lapply(HTO.ridge.plots,plot)
dev.off()

# pdf(file = here::here("analysis","plots","01_build_seurat","HTO_scatter_grobs_run5.pdf"))
# marrangeGrob(HTO.scatter.plots, ncol = 2, nrow = 2)
# dev.off()

pdf(file = here::here("analysis","plots","01_build_seurat","rna_vln_grobs_run5.pdf"))
lapply(RNA.vln.plots, plot)
dev.off()

pdf(file = here::here("analysis","plots","01_build_seurat","rna_count_vln_grobs_run5.pdf"))
lapply(RNA.count.vlnplot, plot)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()