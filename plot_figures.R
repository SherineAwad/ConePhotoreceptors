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
myObject <- readRDS(myRDS) 

DefaultAssay(myObject) <- "RNA"


myGenes <- c("Gls")
figure_name <- paste(mysample, "GLSFP.pdf", sep="")
pdf(file = figure_name, width = 12)
FeaturePlot(myObject , features =myGenes, split.by = "sample", reduction = "umap", pt.size = 0.4)
dev.off() 


figure_name <- paste(mysample, "GLSVP.pdf", sep="")
pdf(file = figure_name, width = 12)
VlnPlot(myObject , features =myGenes, group.by = "sample",cols =  brewer.pal(9, "RdYlBu"),pt.size=0.4)
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


myRDS <- paste(mysample, "_annotated.rds", sep="")
myObject <- readRDS(myRDS) 

table(myObject@meta.data$sample, myObject@active.ident)
#Rod Cones MG AC RPE BC HC Vasculature cells
df <- data.frame(CellType = c("Rod","Cones","MG","AC","RPE","BC","HC","Vasculature cells"),                
day15S1=c(4468,243,259,234,16,655,162,41),
day15S2=c(3385,176,234,176,15,598,153,35),
day30S1=c(6320,244,331,377,24,749,110,24),
day30S2=c(6246,211,341,351,27,789,63,38) )
library(tidyr)
df_long <- df %>% gather(key = "Condition", value = "Value", -CellType)
df_percent <- df_long %>%
  group_by(Condition) %>%
  mutate(Percent = Value / sum(Value) * 100)


stallion = c("AC"="#D51F26","Rod"="#003782","Cone"="#208A42","HC"="#820078","MGPC"="#F47D2B", "MG"="#FFA500","BC"="#8A9FD1","Vasculature cells" ="#E6C122","RPE" ="#4B4BF7")
pdf(file = "Cellratio.pdf", width=4, height=4, onefile=FALSE)
ggplot(df_percent, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" makes it a percent stacked barplot
  scale_y_continuous(labels = scales::percent) +    # Show y-axis as percentage
  labs(y = "Cell Ratio", x = "Sample", fill = "Cell Type") +
  theme_minimal() +
  ggtitle("Celltype Ratio") +
   scale_fill_manual(values = stallion) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1),
   axis.line = element_line(size = 1.5))


dev.off()


