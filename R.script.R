#!/usr/bin/env Rscript

#############################################################
# Author: Gabriela Chavez Calvillo
# Copyright (c) gatoo, `r paste(format(Sys.Date(), "%Y"))`
# Email:  gabvet.c@gmail.com
# 
# Date: `r paste(Sys.Date())`
#
# Script Name: RNA-seq analysis
#
# Script Description: This script was created based on "Differential gene and transcript 
# expression analysis of RNA-seq experiments with TopHat and Cufflinks" on Nature Communications
# the script can run after cufflinks analysis. 

#Install and load packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("cummeRbund")

#I use this script for Mouse transcriptomic analysis

library(cummeRbund)
setwd("/Users/joneslab/Documents/RNA-seq-Data_Analysis/HSV1-mouse")

cuff_data <- readCufflinks('diff_out_Male')
#cuff_data <- readCufflinks('diff_out_Fem')

csDensity(genes(cuff_data))
#csScatter(genes(cuff_data), 'Fem_dLAT', 'Fem_U1')
#csVolcano(genes(cuff_data), 'Fem_dLAT', 'Fem_U1')

#csScatter(genes(cuff_data), 'Male_dLAT', 'Male_U1')
csVolcano(genes(cuff_data), 'Male_dLAT', 'Male_U1')

#Write a table with significant expressed Genes
gene_diff_data <- diffData(genes(cuff_data))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))

#Retrive significant gene IDs (XLOC) with a pre-specified alpha
diffGeneIDs <- getSig(cuff_data,level="genes",alpha=0.05)

#Use returned identifiers to create a CuffGeneSet object with all relevant info for given genes
diffGenes<-getGenes(cuff_data,diffGeneIDs)

#gene_short_name values (and corresponding XLOC_* values) can be retrieved from the CuffGeneSet by using:
names<-featureNames(diffGenes)
row.names(names)=names$tracking_id
diffGenesNames<-as.matrix(names)
diffGenesNames<-diffGenesNames[,-1]

# get the data for the significant genes
diffGenesData<-diffData(diffGenes)
row.names(diffGenesData)=diffGenesData$gene_id
diffGenesData<-diffGenesData[,-1]

# merge the two matrices by row names
diffGenesOutput<-merge(diffGenesNames,diffGenesData,by="row.names")

write.table(diffGenesOutput, 'diff_genes_Male.txt', sep='\t',row.names = F, col.names = T, quote = F)

#gene_short_name values (and corresponding XLOC_* values) can be retrieved from the CuffGeneSet by using:
names<-featureNames(diffGenes)
row.names(names)=names$tracking_id
diffGenesNames<-as.matrix(names)
diffGenesNames<-diffGenesNames[,-1]

# get the data for the significant genes
diffGenesData<-diffData(diffGenes)
row.names(diffGenesData)=diffGenesData$gene_id
diffGenesData<-diffGenesData[,-1]

# merge the two matrices by row names
diffGenesOutput<-merge(diffGenesNames,diffGenesData,by="row.names")

tail(diffGenesOutput)

disp<-dispersionPlot(genes(cuff_data))
disp 
genes.scv<-fpkmSCVPlot(genes(cuff_data))
isoforms.scv<-fpkmSCVPlot(isoforms(cuff_data))

densRep<-csDensity(genes(cuff_data),replicates=T)
densRep
b<-csBoxplot(genes(cuff_data))
b

brep<-csBoxplot(genes(cuff_data),replicates=T)
brep

s<-csScatterMatrix(genes(cuff_data))
s

sp<-csScatter(genes(cuff_data),"Fem_dLAT","Fem_U1",smooth=T)
sp

v<-csVolcanoMatrix(genes(cuff_data))
v

replicates(cuff_data)

#h<-csHeatmap(diffGenesData,cluster='both', replicates=T)
#h

#######
#Building Venn diagrams
#VennDiagram
#install.packages("VennDiagram", "yarr", "RColorBrewer")
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("eulerr", "tidyverse")
#install.packages("tidyverse")
library(VennDiagram)
#install.packages("yarr")
library(yarrr)
library(RColorBrewer)
library(eulerr)
library(tidyverse)

myCol6 <- brewer.pal(6, "Pastel2")

FPKM_all <- read.csv("All-sorted-FPKM.csv", TRUE, ",")

# Prepare a palette of X number of colors with R colorbrewer:
#library(RColorBrewer)
#myCol <- brewer.pal(3, "Pastel2")
myCol3 <- brewer.pal(3, "Pastel2")
myCol4 <- brewer.pal(4, "Pastel2")

##Trying to fix it
#require(venneuler)
#install.packages("venneuler")
m <- (read.csv("All-sorted-FPKM.csv", TRUE, ","))
id <- FPKM_all$X
m2 <- data.matrix(m[,2:6])
rownames(m[,2:6]) -> id

#Plot the venn.diagram and save it as a png file in the working directory for males
venn.diagram(
  x= list(FPKM_all$Male.U1,FPKM_all$Male.WT,FPKM_all$Male.dLAT),
  category.names = c("U1", "WT", "dLAT"),
  filename = "male.png",
  output=TRUE,
  col = myCol3,
  fill= myCol3,
  #Adjusting the size of the font for the numbers
  cex = 1.7,
  fontfamily = "sans",
  #fontface = "bold",
  #Adjusting the size of the font for the category's labels
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)


venn.diagram(
  x= list(FPKM_all$Fem.U1,FPKM_all$Male.U1,FPKM_all$Fem.WT, FPKM_all$Male.WT),
  category.names = c("Fem-U1", "Male-U1","Fem-WT", "Male-WT" ),
  filename = "fem.png",
  output=TRUE,
  col = myCol4,
  fill= myCol4,
  #Adjusting the size of the font for the numbers
  cex = 1.7,
  fontfamily = "sans",
  #fontface = "bold",
  #Adjusting the size of the font for the category's labels
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)


set.seed(12) #for reproducible plot
data.frame(FPKM_all$Male.U1, FPKM_all$Male.WT)%>% #combine the vectors to a data frame
  mutate_at(1:2, as.logical) %>% #convert to logical where ZERO is FALSE and anythingelse is TRUE
  euler(shape = "ellipse", input = "disjoint") %>% #calculate euler object, plot as ellipse
  plot(quantities = T)

#Create a matrix an eliminate zeros to generate a heatmap
m <- (read.csv("All-sorted-FPKM.csv", TRUE, ","))
id <- FPKM_all$X
m2 <- data.matrix(m[,2:7])
rownames(m2[,1:6]) <- id
#Do we wanna remove the ZEROs??
m2[m2==0] <- NA
head(m2)

data2<-m2[complete.cases(m2),]
head(data2)
heatmap(data2)
#Summary: m is the data frame with all information, m2 is the same but in matrix with NA and data2 is without zeros

#Create the correlations
cormat <- round(cor(m2),2)
head(cormat)

#Create the correlation heatmap with ggplot2
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
#Get the lower and the upper triangles of the correlation

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#Usage:
upper_tri <- get_upper_tri(cormat)
upper_tri

#Finish correlation matrix heatmap
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#Reorder correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

#Add correlation coefficients
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

#To create an interactive heatmap 
#install.packages("heatmaply")
library(heatmaply)
library(RColorBrewer)
#Loading data
m <- read.csv("All-sorted-FPKM.csv", TRUE, ",")
id <- m$X                              # assign labels in column 1 to "rnames"
m2 <- data.matrix(m[,2:ncol(m)])       # transform column 2-7 into a matrix
rownames(m2) <- id                     # assign row names
#Normalize data
df <- normalize(m2)

heatmaply(head(df))

#for with different colors

x <- heatmaply(head(df))
y <- heatmaply(tail(df))

x <- heatmaply(
  head(df),
  colors = colorRampPalette(brewer.pal(3, "RdBu"))(256),
  k_col = 2, 
  k_row = 2
)

heatmaply(
  tail(df),
  colors = colorRampPalette(brewer.pal(3, "RdBu"))(256),
  k_col = 2, 
  k_row = 2
)

###To generate a file:
heatmaply(df[1:10,], file = "heatmaply_plot.html") 