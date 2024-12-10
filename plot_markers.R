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


figure_name <- ""
figure_name <- paste(mysample, "_UMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()



figure_name <- "" 
figure_name <- paste(mysample, "_Markers.pdf", sep="") 
pdf(file =figure_name, width =12)
FeaturePlot(myObject, features = c("VWF", "ACKR1","IFI27","Tie1", "Kcnj8", "Pecam1"), order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("Cbln4","Bhlhe23","Cnga3", "Opn1sw", "Opn1mw","Onecut2","Abca8a","Aldh1a1","Kcnj10","Csf2rb","Csf2rb2","Cx3cr1","Mpeg1" ,"Ptprc","Tmem119"),order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("Isl2","Pou4f1","Pou4f2","Pou4f3","Rbpms", "Guca1b","Rom1","Cdk4","Pmel","Rpe65"), order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)  
FeaturePlot(myObject,  features = c("Otx2","HucCD","Crx","Nrl","Rhodopsin","Rxrg","Thrb","Arrestin","Brn3a/b","Atoh7","Isl1"), order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject,  features = c("Olig2","Fonx4","Pou2f2","Rbpj","Sox9","NFI","Lhx1","Calb1","Onecut1","Lhx2","Sox2","Chat","Vsx2","Vsx1","Cabp5","Scgn"),order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("Pax6","Tfap2a","Gad1/2","Prox1","Ascl1","Neurog2","Prkca","Scrt1","Scrt2","Bcl11a","Bcl11b"), order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01) 
FeaturePlot(myObject , features =  c( "Prdm1", "Meis2", "Lhx4", "Ankrd33b"), order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features =  c( "Dscam", "Otor", "Gnat2", "Gngt2"), order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features =  c( "Alk", "Rgr", "Slit2", "Lrat"),order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("Rbfox3", "Sebox", "Gad1", "Elavl3"),order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("Sox9", "Glul",  "Rlbp1", "Ascl1"), order=TRUE, reduction =  "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c( "Otx2", "Olig2", "Crx", "Neurog2"),order=TRUE, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("Rho", "Arr3", "Tfap2b", "Vsx1"), order=TRUE, reduction =  "umap", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("Insm1", "Prdm1", "GFP", "Elavl4"), order=TRUE, reduction =  "umap", cols = c("lightgrey", "red"), pt.size = 0.01)


dev.off()







