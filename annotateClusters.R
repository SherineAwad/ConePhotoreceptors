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

myRDS <- paste(mysample, "_analysed.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

#head(myObject) 
#levels(myObject)
#head(Idents(myObject) 

table(myObject@active.ident)
table(myObject@meta.data[,'sample'])


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

#get no. of cells per sample and cells 
table(myObject@meta.data$sample, myObject@meta.data$cells)

df <- data.frame(CellType = c("AC", "Rod", "Cone", "HC", "MG", "BC", "Vasculature cells", "RPE"), day15S1=c(296,4467,243,114,199,641,41,16),
day15S2=c(213,3390,176,133,180,581,35,15),
day30S1=c(393,6315,244,95,198,748,24,24),
day30S2=c(324,6239,212,102,199,780,38,27) )
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



myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(mySubset, file = myRDS)



