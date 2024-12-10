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


myGenes <- c("Galnt1", "Galnt10", "Galnt11", "Galnt12", "Galnt13", "Galnt14", "Galnt2", "Galnt3", "Galnt4", "Galnt5", "Galnt6", "Galnt7", "Galnt9", "Galnt16", "Galntl5", "Galntl6", "Galnt17") 
figure_name <- ""
figure_name <- paste(mysample, "N-Acetylgalactosaminyltransferases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red")) 
dev.off()

myGenes <- c("A4gnt", "B3gnt2", "B3gnt3", "B3gnt4", "B3gnt8", "Gcnt1", "Gcnt3", "Mgat1", "Mgat2", "Mgat3", "Mgat4a", "Mgat4b", "Mgat4c", "Mgat5", "Mgat5b", "Ogt", "Pomgnt1")

figure_name <- ""
figure_name <- paste(mysample, "N-Acetylglucosaminyltransferases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red")) 
dev.off()

myGenes <- c("B3glct", "B4galt1", "B4galt2", "B4galt3", "B4galt5", "C1galt1") 
figure_name <- paste(mysample, "Galactosyltransferases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Uggt1", "Uggt2")
figure_name <- paste(mysample, "Glucosyltransferases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Edem1", "Edem2", "Edem3", "Man1a", "Man1a2", "Man1b1", "Man1c1", "Man2a1", "Man2a2", "Man2b1", "Manba")
figure_name <- paste(mysample, "Mannosidases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("A4gnt", "B3glct", "B3gnt8", "B4galt5", "C1galt1", "C1galt1c1", "Galnt1", "Galnt10", "Galnt11", "Galnt12", "Galnt13", "Galnt14", "Galnt2", "Galnt3", "Galnt4", "Galnt5", "Galnt6", "Galnt7", "Galnt9", "Galnt16", "Galntl5", "Galntl6", "Gcnt1", "Gcnt3", "Ogt", "Pofut1", "Pofut2", "Pomgnt1", "Pomt1", "Pomt2", "St3gal1", "St3gal2", "St6galnac1", "St8sia3", "St8sia6", "Galnt17")

figure_name <- paste(mysample, "O-Linked Glycosylation.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Aga", "B3gnt2", "B3gnt3", "B3gnt8", "B4galt1", "B4galt2", "B4galt3", "Edem1", "Edem2", "Edem3", "Fuca1", "Fuca2", "Fut11", "Fut8", "Ganab", "Glb1", "Gnptab", "Gnptg", "Hexa", "Hexb"," Man1a", "Man1a2", "Man1b1", "Man1c1", "Man2a1", "Man2a2", "Man2b1", "Manba", "Mgat1", "Mgat2", "Mgat3", "Mgat4a", "Mgat4b", "Mgat4c", "Mgat5", "Mgat5b", "Mogs", "Nagpa", "Neu1", "Neu2", "Neu3", "Neu4", "Prkcsh", "St6gal1", "St8sia2", "St8sia3", "St8sia4", "St8sia6", "Uggt1", "Uggt2")

figure_name <- paste(mysample, "N-Linked Glycosylation.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Aga", "C1galt1c1")
figure_name <- paste(mysample, "OtherGlycosylation.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("B3gnt3", "B3gnt4", "B4galt1", "B4galt2", "B4galt3", "Glb1", "Hexa", "Hexb", "St3gal1", "St3gal2", "St8sia6", "Galnt17") 
figure_name <- paste(mysample, "GlycosphingolipidSynthesis.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Gnptab", "Gnptg", "Nagpa")
figure_name <- paste(mysample, "Mannose6PhosphateSynthesisCatabolism.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("St3gal1", "St3gal2", "St6gal1", "St6galnac1", "St8sia2", "St8sia3", "St8sia4", "St8sia6")
figure_name <- paste(mysample, "Sialyltransferases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Neu1", "Neu2", "Neu3", "Neu")
figure_name <- paste(mysample, "Sialidases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Fuca1", "Fuca2", "Fut11", "Fut8", "Pofut1", "Pofut2")
figure_name <- paste(mysample, "FucosidasesFucosyltransferases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Ganab", "Glb1", "Hexa", "Hexb", "Mogs", "Prkcsh")
figure_name <- paste(mysample, "GalactosidesGlucosidasesHexosaminidases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Pomt1", "Pomt2") 
figure_name <- paste(mysample, "Mannosyltransferases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Gfpt1", "Gnpnat1", "Pgm3", "Uap1", "Ogt", "Oga", "Hk1", "Hk2", "Gpl") 
figure_name <- paste(mysample, "Hexoaminepathway.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Fas","Parp1","Casp8","Casp9","Casp3","Ripk1","Ripk3","Atg5","Sqstm1","Alox12","Alox15","Acsl4","Lpcat3","Gch1","Hmgcr","Chac1","Ptgs2","Rpl8")
figure_name <- paste(mysample, "Apoptosis.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


myGenes <- c("Slc1a5","Slc7a5","Gls","Gls2","Glud1","Asns")
figure_name <- paste(mysample, "glutaminetransportcatabolism.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Got1","Got2","Bcat1","Bcat2","Me1","Me2")
figure_name <- paste(mysample, "transaminases.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Gpx4","Gpx1","Aifm2","Sod1","Sod2","Cat","G6PD","Pgd","Shmt1","Shmt2","Mthfd1","Gsr","Gclc","Slc7a11","Gss")
figure_name <- paste(mysample, "reoxbalance.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Idh3g","Ogdh","Sucla2","Sucgl1","Sucgl2","Sdha","Sdhb","Sdhc","Sdhd","Fh","Mdh1","Mdh1b","Mdh2","Cs","Acly","Pdha","Pdhb","Hk2","Pdk1","Pdk2","Pdk3","Pdk4","Pcx")
figure_name <- paste(mysample, "TCA.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

myGenes <- c("Tfrc","Slc40a1","Ftl1","Fth1","Ncoa4","Slc11a2","Slc39a14","Slc39a8","Hamp","Ireb2","Hif1α","Epas1","Hmox1","Heph","Cp")
figure_name <- paste(mysample, "Ironmetabolism.pdf", sep="")
pdf(file =figure_name, width=18)
DoHeatmap(object = myObject, features = myGenes, size = 2.5, group.by ="sample",  slot="data",group.colors = brewer.pal(9, "Blues")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


