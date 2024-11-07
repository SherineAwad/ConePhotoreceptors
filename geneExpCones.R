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


cones <- subset(myObject, idents = "Cone", invert = FALSE)
table(cones@active.ident)
table(cones@meta.data[,'sample'])
cones@active.ident <- as.factor(cones@meta.data[,'sample'])
names(cones@active.ident) <- rownames(cones@meta.data)

DEGs <- FindAllMarkers(object = cones, only.pos = TRUE, min.pct = 0.1,test.use ='wilcox', logfc.threshold = 0.5)
write.csv(DEGs, "photoreceptors_conesDGESper4.csv") 

top_genes <- DEGs[order(DEGs$p_val), ][1:50, ]  # Top 10 genes based on smallest p-value
top_gene_list <- rownames(top_genes)

figure_name <- ""
figure_name <- paste(mysample, "heatmapPerSamples.pdf", sep="")
pdf(file =figure_name, width =12)
DoHeatmap(cones,features = top_gene_list,size = 4, group.by = "sample", slot="data") + NoLegend()+ theme(axis.text.y = element_text(size = 5)+ scale_fill_gradientn(colors = c("blue", "white", "red")))
dev.off()

cones <- NormalizeData(object = cones, normalization.method = "LogNormalize", scale.factor = 10000)
AverageExpression <- AverageExpression(object = cones,assays = "RNA",  group.by ="sample")
write.csv(AverageExpression,"photoreceptors_conesAvgExpper4.csv") 
DefaultAssay(myObject) <- "RNA"


cones @meta.data[,'sample'] <- recode(cones@meta.data[,'sample'], "15dayS1" = "Ctrl", "30dayS1" = "Ctrl")

cones @meta.data[,'sample'] <- recode(cones@meta.data[,'sample'], "15dayS2" = "Ctrl", "30dayS2" = "KO")

cones @active.ident  <- recode(cones @active.ident, "15dayS1" = "Ctrl", "30dayS1" = "Ctrl")

cones @active.ident  <- recode(cones @active.ident, "15dayS2" = "Ctrl", "30dayS2" = "KO")


DEGs <- FindAllMarkers(object = cones, only.pos = TRUE, min.pct = 0.1,test.use ='wilcox', logfc.threshold = 0.5)
write.csv(DEGs, "photoreceptors_conesDGES.csv")

top_genes <- DEGs[order(DEGs$p_val), ][1:50, ]  # Top 10 genes based on smallest p-value
top_gene_list <- rownames(top_genes)

figure_name <- ""
figure_name <- paste(mysample, "heatmap.pdf", sep="")
pdf(file =figure_name, width =12)
DoHeatmap(cones,features = top_gene_list, size = 4,group.by = "sample", slot="data") + NoLegend()+ theme(axis.text.y = element_text(size = 5)+ scale_fill_gradientn(colors = c("blue", "white", "red")))
dev.off()

cones <- NormalizeData(object = cones, normalization.method = "LogNormalize", scale.factor = 10000)
AverageExpression <- AverageExpression(object = cones,assays = "RNA",  group.by ="sample")
write.csv(AverageExpression,"photoreceptors_conesAvgExp.csv")

myRDS <- paste(mysample, "_cones.rds", sep="")
saveRDS(cones, file = myRDS)








