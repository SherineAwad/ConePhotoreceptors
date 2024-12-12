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

cones <- readRDS(myRDS)
filename = ""
filename = paste(mysample,"conesDGEsPerSamples.csv", sep="")
DEGs <- FindAllMarkers(object = cones, only.pos = TRUE, min.pct = 0.1,test.use ='wilcox', logfc.threshold = 0.5)
write.csv(DEGs, filename) 

cones @meta.data[,'sample'] <- recode(cones@meta.data[,'sample'], "15dayS1" = "Ctrl", "30dayS1" = "Ctrl")

cones @meta.data[,'sample'] <- recode(cones@meta.data[,'sample'], "15dayS2" = "KO", "30dayS2" = "KO")

cones @active.ident  <- recode(cones @active.ident, "15dayS1" = "Ctrl", "30dayS1" = "Ctrl")

cones @active.ident  <- recode(cones @active.ident, "15dayS2" = "KO", "30dayS2" = "KO")


DEGs <- FindAllMarkers(object = cones, only.pos = TRUE, min.pct = 0.1,test.use ='wilcox', logfc.threshold = 0.5)
filename =""
filename <- paste(mysample, "conesDEGs.csv", sep="")
write.csv(DEGs, filename)


top_genes<- DEGs%>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top_gene_list <- top_genes$gene
top_gene_list 

figure_name <- ""
figure_name <- paste(mysample, "DEGsDP.pdf", sep="")
pdf(file =figure_name, width=15, height=5)
DotPlot(cones, features = top_gene_list,  col.min = -100, col.max =100,group.by="sample", dot.scale = 12) + RotatedAxis()
dev.off()


cones <- NormalizeData(object = cones, normalization.method = "LogNormalize", scale.factor = 10000)
AverageExpression <- AverageExpression(object = cones,assays = "RNA",  group.by ="sample")
filename =""
filename <- paste(mysample, "_coneAvgExp.csv", sep="")
write.csv(AverageExpression, filename)

myRDS <- paste(mysample, "_cones.rds", sep="")
saveRDS(cones, file = myRDS)








