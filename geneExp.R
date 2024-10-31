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
write.csv(DEGs, "photoreceptors_conesDGES.csv") 

cones <- NormalizeData(object = cones, normalization.method = "LogNormalize", scale.factor = 10000)
AverageExpression <- AverageExpression(object = cones,assays = "RNA",  group.by ="sample")
write.csv(AverageExpression,"photoreceptors_conesAvgExp.csv") 


myRDS <- paste(mysample, "_cones.rds", sep="")
saveRDS(cones, file = myRDS)








