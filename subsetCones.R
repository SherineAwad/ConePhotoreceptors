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

cones <- FindClusters(cones, resolution = 2.0)
cones <- RunUMAP(cones, dims = 1:8, reduction.name = "umap")


figure_name <- ""
figure_name <- paste(mysample, "_reclusteredConesUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(cones, reduction = "umap",  repel = TRUE) + ggtitle("UMAP")
DimPlot(cones, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


cell_values <- c("0", "1", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
mySubset <- subset(cones, idents = cell_values, invert = FALSE)

mySubset@active.ident <- as.factor(mySubset@meta.data[,'sample'])
names(mySubset@active.ident) <- rownames(mySubset@meta.data)


myRDS <- paste(mysample, "_reClusteredCones.rds", sep="")
saveRDS(mySubset, file = myRDS)








