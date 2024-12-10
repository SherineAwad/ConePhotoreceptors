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
myRDS <- args[1]

split_string <- strsplit(myRDS, "_")[[1]]
mysample <- split_string[1]
print(mysample)

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"

table(myObject@active.ident)
table(myObject@meta.data[,'sample'])


myObject <- RenameIdents(
  object = myObject,
  "0" = 'Rod', 
  "1" = 'Rod', 
  "2" = 'Rod',
  "3" = 'Rod',
  "4" = 'Rod',
  "5" = 'Rod',
  "6" = 'Rod',
  "7" = 'Rod',
  "9" = 'Rod',
  "10" = 'Rod', 
  "15" = 'Rod',

  "32" = "Cones", 
  "14" = "Cones",
  "16" = "Cones",

  "8" = "MG", 
  "26" = "MG", 
  "25" = "MG",
  "38" = "MG", 

  "11" = 'AC', 
  "17" = 'AC',  
  "29" = 'AC', 
  "37" = 'AC', 
  "31" = 'AC', 
  
  "35" = 'RPE',

  "19" = "BC", 
  "36" = "BC",
  "20" = "BC",
  "24" = "BC",
  "21" = "BC",
  "12" = "BC",
  "22" = "BC",
  "27" = "BC",
  "34" = "BC",
  "13"= "BC",
  "18" = "BC",

  "30" = "HC",
  "23" = "HC",
  "33" = "HC", 

  "28" = "Vasculature cells")


figure_name <- ""
figure_name <- paste(mysample, "_AnnotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(myObject, file = myRDS)






