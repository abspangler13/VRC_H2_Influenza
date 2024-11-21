library("VennDiagram") 

load(file = here::here("analysis","data_objects","06_DE","fgsea_padj_full_glist","Hallmark_sig_pathways_data.Rdata"))
padj.sig.path <- sig.path

load(file = here::here("analysis","data_objects","06_DE","fgsea_fc_only","Hallmark_sig_pathways_data.Rdata"))
fc.sig.path <- sig.path

padj <- nrow(padj.sig.path)
fc <- nrow(fc.sig.path)
int <- length(intersect(padj.sig.path$pathway,fc.sig.path$pathway))

grid.newpage() 

pdf(file = here::here("analysis","plots","06_DE","hallmark_venn.pdf"))
draw.pairwise.venn(area1=padj, area2=fc,cross.area=int, 
                   category=c("padj_fc","fc"),fill=c("Red","Yellow"), cex = 3)
dev.off()

###############

load(file = here::here("analysis","data_objects","06_DE","fgsea_padj_full_glist","Immunologic Signature_sig_pathways_data.Rdata"))
padj.sig.path <- sig.path

load(file = here::here("analysis","data_objects","06_DE","fgsea_fc_only","Immunologic Signature_sig_pathways_data.Rdata"))
fc.sig.path <- sig.path

padj <- nrow(padj.sig.path)
fc <- nrow(fc.sig.path)
int <- length(intersect(padj.sig.path$pathway,fc.sig.path$pathway))

grid.newpage() 

pdf(file = here::here("analysis","plots","06_DE","Immunologic Signature_venn.pdf"))
draw.pairwise.venn(area1=padj, area2=fc,cross.area=int, 
                   category=c("padj_fc","fc"),fill=c("Red","Yellow"),cex = 3)
dev.off()

###############

load(file = here::here("analysis","data_objects","06_DE","fgsea_padj_full_glist","Reactome_sig_pathways_data.Rdata"))
padj.sig.path <- sig.path

load(file = here::here("analysis","data_objects","06_DE","fgsea_fc_only","Reactome_sig_pathways_data.Rdata"))
fc.sig.path <- sig.path

padj <- nrow(padj.sig.path)
fc <- nrow(fc.sig.path)
int <- length(intersect(padj.sig.path$pathway,fc.sig.path$pathway))

grid.newpage() 

pdf(file = here::here("analysis","plots","06_DE","Reactome_venn.pdf"))
draw.pairwise.venn(area1=padj, area2=fc,cross.area=int, 
                   category=c("padj_fc","fc"),fill=c("Red","Yellow"),cex = 3)
dev.off()

###############

load(file = here::here("analysis","data_objects","06_DE","fgsea_padj_full_glist","Ontology_sig_pathways_data.Rdata"))
padj.sig.path <- sig.path

load(file = here::here("analysis","data_objects","06_DE","fgsea_fc_only","Ontology_sig_pathways_data.Rdata"))
fc.sig.path <- sig.path

padj <- nrow(padj.sig.path)
fc <- nrow(fc.sig.path)
int <- length(intersect(padj.sig.path$pathway,fc.sig.path$pathway))

grid.newpage() 

pdf(file = here::here("analysis","plots","06_DE","Ontology_venn.pdf"))
draw.pairwise.venn(area1=padj, area2=fc,cross.area=int, 
                   category=c("padj_fc","fc"),fill=c("Red","Yellow"),cex = 3)
dev.off()