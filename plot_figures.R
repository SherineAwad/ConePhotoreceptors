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
pdf(file =figure_name, width=8, height=5)
DotPlot(myObject, features = myGenes, col.min = -100, col.max =100, group.by="sample", dot.scale = 8) + RotatedAxis()
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
  mutate(Percent = Value / sum(Value) * 100, ymax = cumsum(Percent), ymin = c(0, head(ymax, n=-1) ) ) 

stallion = c("AC"="#D51F26","Rod"="#003782","Cones"="#208A42","HC"="#820078","MGPC"="#F47D2B", "MG"="#FFA500","BC"="#8A9FD1","Vasculature cells" ="#E6C122","RPE" ="#4B4BF7")
pdf(file = "Cellratio.pdf", width=4, height=4, onefile=FALSE)
custom_order <- c("Cones", "Rod", "AC", "BC", "HC", "RPE", "MG", "Vasculature cells")  
df_percent$CellType <- factor(df_percent$CellType, levels = custom_order)

ggplot(df_percent, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" makes it a percent stacked barplot
  scale_y_continuous(labels = scales::percent) +    # Show y-axis as percentage
  labs(y = "Cell Ratio", x = "Sample", fill = "Cell Type") +
  theme_minimal() +
  ggtitle("Celltype Ratio") +
  scale_fill_manual(values = stallion) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  axis.line = element_line(size = 1.5)  )
dev.off()


ggplot(df_percent, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity") +  # Create stacked bar plot
  geom_text(aes(y = (ymax + ymin) / 2,  # Position labels in the middle of each segment
                label = paste0(round(Percent, 1), "%")),  # Format labels as percentages
            color = "white", size = 1.5) +  # Text properties
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Add percentage scale on y-axis
  theme_minimal() + labs(y = "Cell Ratio", x = "Sample", fill = "Cell Type") +theme_minimal() +
  ggtitle("Celltype Ratio") +
  scale_fill_manual(values = stallion) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  axis.line = element_line(size = 1.5)  )
dev.off()



grad_colors <- colorRampPalette(c("blue", "white", "red"))(100)


ConesGenes = c("Arr3", "Gls", "Opn1sw", "Rxrg", "Gnat2", "Pde6h", "Cngb3")
data <- read.csv("photoreceptors_coneAvgExp.csv", header = TRUE, row.names = 1)
df <- as.matrix(data)
df_selected <- df[rownames(df) %in% ConesGenes , ]

figure_name <- ""
figure_name <- paste(mysample, "_conesAvgExpHeatmap.pdf", sep="")
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
main= "Cones Genes"
)
dev.off()



