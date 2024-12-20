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
library(heatmaply)

mysample =  "photoreceptors" 
grad_colors <- colorRampPalette(c("blue", "white", "red"))(100)
data <- read.csv("photoreceptors_coneAvgExp.csv", header = TRUE, row.names = 1)
df <- as.matrix(data)

myGenes <- c("Galnt1", "Galnt10", "Galnt11", "Galnt12", "Galnt13", "Galnt14", "Galnt2", "Galnt3", "Galnt4", "Galnt5", "Galnt6", "Galnt7", "Galnt9", "Galnt16", "Galntl5", "Galntl6", "Galnt17")
df_selected <- df[rownames(df) %in% myGenes, ]

figure_name <- ""
figure_name <- paste(mysample, "N-Acetylgalactosaminyltransferases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none", 
col = grad_colors, 
margins = c(7, 7),  
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),   
labCol = c("Ctrl", "KO"),
xlab = "Samples", 
ylab = "Genes",
trace = "none", 
key = TRUE, 
keysize = 1.2, 
density.info ="none",
dendrogram="none",
main= "N-Acetylgalactosaminyltransferases"
)
dev.off() 



myGenes <- c("A4gnt", "B3gnt2", "B3gnt3", "B3gnt4", "B3gnt8", "Gcnt1", "Gcnt3", "Mgat1", "Mgat2", "Mgat3", "Mgat4a", "Mgat4b", "Mgat4c", "Mgat5", "Mgat5b", "Ogt", "Pomgnt1")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "N-Acetylglucosaminyltransferases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "N-Acetylglucosaminyltransferases"
)
dev.off()


myGenes <- c("B3glct", "B4galt1", "B4galt2", "B4galt3", "B4galt5", "C1galt1") 
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Galactosyltransferases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Galactosyltransferases"
)
dev.off()


myGenes <- c("Uggt1", "Uggt2")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Glucosyltransferases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2, 
density.info ="none",
dendrogram="none",
main= "Glucosyltransferases"
)
dev.off()


myGenes <- c("Edem1", "Edem2", "Edem3", "Man1a", "Man1a2", "Man1b1", "Man1c1", "Man2a1", "Man2a2", "Man2b1", "Manba")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Mannosidases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
rowsep,
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Mannosidases"
)
dev.off()



myGenes <- c("A4gnt", "B3glct", "B3gnt8", "B4galt5", "C1galt1", "C1galt1c1", "Galnt1", "Galnt10", "Galnt11", "Galnt12", "Galnt13", "Galnt14", "Galnt2", "Galnt3", "Galnt4", "Galnt5", "Galnt6", "Galnt7", "Galnt9", "Galnt16", "Galntl5", "Galntl6", "Gcnt1", "Gcnt3", "Ogt", "Pofut1", "Pofut2", "Pomgnt1", "Pomt1", "Pomt2", "St3gal1", "St3gal2", "St6galnac1", "St8sia3", "St8sia6", "Galnt17")
figure_name <- ""
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "O-LinkedGlycosylation.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "O-LinkedGlycosylation"
)
dev.off()


myGenes <- c("Aga", "B3gnt2", "B3gnt3", "B3gnt8", "B4galt1", "B4galt2", "B4galt3", "Edem1", "Edem2", "Edem3", "Fuca1", "Fuca2", "Fut11", "Fut8", "Ganab", "Glb1", "Gnptab", "Gnptg", "Hexa", "Hexb"," Man1a", "Man1a2", "Man1b1", "Man1c1", "Man2a1", "Man2a2", "Man2b1", "Manba", "Mgat1", "Mgat2", "Mgat3", "Mgat4a", "Mgat4b", "Mgat4c", "Mgat5", "Mgat5b", "Mogs", "Nagpa", "Neu1", "Neu2", "Neu3", "Neu4", "Prkcsh", "St6gal1", "St8sia2", "St8sia3", "St8sia4", "St8sia6", "Uggt1", "Uggt2")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "N-LinkedGlycosylation.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "N-LinkedGlycosylation"
)
dev.off()




myGenes <- c("Aga", "C1galt1c1")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "OtherGlycosylation.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "OtherGlycosylation"
)
dev.off()



myGenes <- c("B3gnt3", "B3gnt4", "B4galt1", "B4galt2", "B4galt3", "Glb1", "Hexa", "Hexb", "St3gal1", "St3gal2", "St8sia6", "Galnt17") 
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "GlycosphingolipidSynthesis.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "GlycosphingolipidSynthesis"
)
dev.off()


myGenes <- c("Gnptab", "Gnptg", "Nagpa")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Mannose6PhosphateSynthesisCatabolism.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Mannose6PhosphateSynthesisCatabolism"
)
dev.off()

myGenes <- c("St3gal1", "St3gal2", "St6gal1", "St6galnac1", "St8sia2", "St8sia3", "St8sia4", "St8sia6")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Sialyltransferases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Sialyltransferases"
)
dev.off()

myGenes <- c("Neu1", "Neu2", "Neu3", "Neu")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Sialidases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Sialidases"
)
dev.off()

myGenes <- c("Fuca1", "Fuca2", "Fut11", "Fut8", "Pofut1", "Pofut2")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "FucosidasesFucosyltransferases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "FucosidasesFucosyltransferases"
)
dev.off()

myGenes <- c("Ganab", "Glb1", "Hexa", "Hexb", "Mogs", "Prkcsh")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "GalactosidesGlucosidasesHexosaminidases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "GalactosidesGlucosidasesHexosaminidases"
)
dev.off()

myGenes <- c("Pomt1", "Pomt2") 
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Mannosyltransferases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Mannosyltransferases"
)
dev.off()

myGenes <- c("Gfpt1", "Gnpnat1", "Pgm3", "Uap1", "Ogt", "Oga", "Hk1", "Hk2", "Gpl") 
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Hexoaminepathway.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Hexoaminepathway"
)
dev.off()

myGenes <- c("Fas","Parp1","Casp8","Casp9","Casp3","Ripk1","Ripk3","Atg5","Sqstm1","Alox12","Alox15","Acsl4","Lpcat3","Gch1","Hmgcr","Chac1","Ptgs2","Rpl8")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Apoptosis.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Apoptosis"
)
dev.off()

myGenes <- c("Slc1a5","Slc7a5","Gls","Gls2","Glud1","Asns")
df_selected <- df[rownames(df) %in% myGenes, ]

figure_name <- ""
figure_name <- paste(mysample, "glutaminetransportcatabolism.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "glutaminetransportcatabolism"
)
dev.off()


myGenes <- c("Got1","Got2","Bcat1","Bcat2","Me1","Me2")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "transaminases.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "transaminases"
)
dev.off()

myGenes <- c("Gpx4","Gpx1","Aifm2","Sod1","Sod2","Cat","G6PD","Pgd","Shmt1","Shmt2","Mthfd1","Gsr","Gclc","Slc7a11","Gss")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "reoxbalance.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "reoxbalance"
)
dev.off()

myGenes <- c("Idh3g","Ogdh","Sucla2","Sucgl1","Sucgl2","Sdha","Sdhb","Sdhc","Sdhd","Fh","Mdh1","Mdh1b","Mdh2","Cs","Acly","Pdha","Pdhb","Hk2","Pdk1","Pdk2","Pdk3","Pdk4","Pcx")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "TCA.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "TCA"
)
dev.off()

myGenes <- c("Tfrc","Slc40a1","Ftl1","Fth1","Ncoa4","Slc11a2","Slc39a14","Slc39a8","Hamp","Ireb2","Hif1α","Epas1","Hmox1","Heph","Cp")
df_selected <- df[rownames(df) %in% myGenes, ]
figure_name <- ""
figure_name <- paste(mysample, "Ironmetabolism.pdf", sep="")
pdf(file =figure_name)
gplots::heatmap.2(df_selected,
scale = "none",
col = grad_colors,
margins = c(7, 7),
cexRow = 0.8,
cexCol = 0.8,
labRow = rownames(df_selected),
labCol = c("Ctrl", "KO"),
xlab = "Samples",
ylab = "Genes",
trace = "none",
key = TRUE,
keysize = 1.2,
density.info ="none",
dendrogram="none",
main= "Ironmetabolism"
)
dev.off()
