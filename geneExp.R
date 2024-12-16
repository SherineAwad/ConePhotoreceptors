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
library(RColorBrewer)



args <- commandArgs(trailingOnly = TRUE)
myRDS <- args[1]

split_string <- strsplit(myRDS, "_")[[1]]
mysample <- split_string[1]
print(mysample)


myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"



myObject[["cells"]] <- Idents(object = myObject)

DEGs <- FindAllMarkers(object = myObject, only.pos = TRUE, min.pct = 0.1,test.use ='wilcox', logfc.threshold = 0.5)
top_genes<- DEGs%>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top_gene_list <- top_genes$gene
top_gene_list 


figure_name <- ""
figure_name <- paste(mysample, "allDEGsDP.pdf", sep="")
pdf(file =figure_name, width=18, height=12)
DotPlot(myObject, features = top_gene_list,  col.min = -100, col.max =100,group.by="cells", dot.scale = 12)  + RotatedAxis()
dev.off()





figure_name <- ""
figure_name <- paste(mysample, "allDEGsDPsplitSample.pdf", sep="")
pdf(file =figure_name, width=18, height=12)
DotPlot(myObject, features = top_gene_list,  col.min = -100, col.max =100, split.by ="sample", group.by="cells",cols =  brewer.pal(9, "RdYlBu"),  dot.scale = 12)  + RotatedAxis()
dev.off()





