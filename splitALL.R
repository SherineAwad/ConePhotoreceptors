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


#After analysing and merging 
myRDS <- readRDS("photoreceptors.rds")
groups <- sample(c("15dayS1", "15dayS2", "30dayS1", "30dayS2"), size =27095 ,replace = TRUE)
names(groups) <- colnames(myObject)
myObject  <-  AddMetaData(object = myObject, metadata = groups, col.name = "sample")
saveRDS(myObject, file = "photoreceptors.rds")





