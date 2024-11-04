library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

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
figure_name <- paste(mysample, "GLSFP.png", sep="")
png(filename = figure_name, width = 1500, height = 1500, units = "px")
FeaturePlot(myObject , features =myGenes, split.by = "sample", reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
dev.off() 


figure_name <- paste(mysample, "GLSVP.png", sep="")
png(filename = figure_name, width = 1500, height = 1500, units = "px")
VlnPlot(myObject , features =myGenes, split.by = "sample")
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
figure_name <- paste(mysample, "GLSDP.png", sep="")
png(filename = figure_name, width = 1500, height = 1500, units = "px")
DotPlot(myObject, features = myGenes, split.by ="sample") + RotatedAxis()
dev.off()


