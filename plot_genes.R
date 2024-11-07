library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
library(RColorBrewer)
library(viridis)


args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
myRDS <- paste(mysample, "_cones.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

figure_name <- ""
figure_name <- paste(mysample, "Cones_UMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP") 
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()

myGenes <- c("Gls")
figure_name <- paste(mysample, "GLSFP.pdf", sep="")
pdf(file = figure_name, width = 12)
FeaturePlot(myObject , features =myGenes, split.by = "sample", reduction = "umap", pt.size = 0.4)
dev.off() 


figure_name <- paste(mysample, "GLSVP.pdf", sep="")
pdf(file = figure_name, width = 12)
VlnPlot(myObject , features =myGenes, split.by = "sample",cols =  brewer.pal(9, "RdYlBu"),pt.size=0.4)
dev.off()



myGenes <- c("Gls","Arr3",
"Gnat2",
"Opn1mw",
"Opn1sw",
"Nrxn3",
"Wdfy1",
"Iqgap2",
"Gulp1",
"Kcng1",
"Unc5c")
figure_name <- ""
figure_name <- paste(mysample, "GLSDP.pdf", sep="")
pdf(file = figure_name, width = 12)
DotPlot(myObject, features = myGenes, split.by ="sample", group.by="sample", dot.scale = 12) + RotatedAxis()  
dev.off()



