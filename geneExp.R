library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)

mysample = args[1]
myRDS <- paste(mysample, "_annotated.rds", sep="")
mysample
myRDS

myObject <- readRDS(myRDS)

myObject[["cells"]] <- Idents(object = myObject)

DEGs <- FindAllMarkers(object = myObject, only.pos = TRUE, min.pct = 0.1,test.use ='wilcox', logfc.threshold = 0.5)

top_genes <- DEGs[order(DEGs$p_val), ][1:50, ]  # Top 10 genes based on smallest p-value
top_gene_list <- rownames(top_genes)

top_gene_list 

figure_name <- ""
figure_name <- paste(mysample, "Cellsheatmap.pdf", sep="")
pdf(file =figure_name, width =20)
DoHeatmap(object = myObject, features = top_gene_list,size = 2.5, group.by ="cells",  slot="data") + theme(axis.text.y = element_text(size = 2.5)+ scale_fill_gradientn(colors = c("blue", "white", "red")))  
dev.off()


myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(myObject, file = myRDS)








