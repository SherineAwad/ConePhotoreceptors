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

mysample = args[1]
myRDS <- paste(mysample, "_annotated.rds", sep="")
mysample
myRDS

myObject <- readRDS(myRDS)

myObject[["cells"]] <- Idents(object = myObject)

DEGs <- FindAllMarkers(object = myObject, only.pos = TRUE, min.pct = 0.1,test.use ='wilcox', logfc.threshold = 0.5)

top_genes<- DEGs%>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top_gene_list <- top_genes$gene



top_gene_list 

figure_name <- ""
figure_name <- paste(mysample, "Cellsheatmap.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = top_gene_list,size = 2.5, group.by ="cells",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red")) 
dev.off()

myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(myObject, file = myRDS)








