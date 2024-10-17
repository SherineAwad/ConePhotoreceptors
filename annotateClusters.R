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

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

myRDS <- paste(mysample, "_split.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

#head(myObject) 
#levels(myObject)
#head(Idents(myObject) 


myObject <- RenameIdents(
  object = myObject,
  "20" = 'AC', 
  "23" = 'AC',
  "13" = 'AC',  
  "27" = 'AC',
  "0" = 'Rod', 
  "1" = 'Rod', 
  "2" = 'Rod', 
  "3" = 'Rod', 
  "4" = 'Rod', 
  "5" = 'Rod', 
  "6" = 'Rod', 
  "7" = 'Rod', 
  "8" = 'Rod',
  "9" = 'Rod', 
  "14" = 'Rod', 
  "11" = 'Cone', 
  "21" = 'Cone', 
  "33" = 'Cone', 
  "28" = 'HC', 
  "31" = 'HC', 
  "34" = 'HC',
  "38" = 'HC', 
  "26" = 'doublets', 
  "10" = 'MG', 
  "35" = 'BC', 
  "12" = 'BC', 
  "22" = 'BC', 
  "15" = 'BC', 
  "16" = 'BC', 
  "17" = 'BC', 
  "18" = 'BC', 
  "19" = 'BC', 
  "25" = 'BC', 
  "37" = 'BC', 
  "29" = 'BC',
  "30" = 'Vasculature cells',
  "36" = 'RPE',
  "24" = 'doublets',
  "28" = 'doublets',  
  "32"= 'AC')  
head(Idents(myObject))


cell_values <- c("AC", "Rod", "Cone", "HC", "MG", "BC", "RPE", "Vasculature cells")
mySubset <- subset(myObject, idents = cell_values, invert = FALSE)


figure_name <- ""
figure_name <- paste(mysample, "_AnnotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(mySubset, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(mySubset, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()

myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(mySubset, file = myRDS)



