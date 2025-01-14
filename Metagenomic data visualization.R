###########################################################################################################################################################################
# This is an R script
#Programmer: Wenqiao He
#Last update: December 2024
#Purpose: Visualization of Aedes mosquito virus metagenomic data analysis
###########################################################################################################################################################################

install.packages("eulerr")
install.packages("viridis")
install.packages("VennDiagram")
install.packages("gplots")
install.packages("tidyverse")
install.packages("data.table")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("stringr")
install.packages("tmap")
install.packages("tmaptools")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("writexl")
library(viridis)
library(eulerr)
library(VennDiagram)
library(readxl)
library(writexl)
library(gplots)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(FactoMineR)
library(factoextra)
library(stringr)
library(tmap)       
library(tmaptools)  
library(ggrepel)

### Set working directory###
setwd("/PATH/")

### Load data for Venn analysis at family level
family_summary <- read_excel("Venn_analysis_family_example_data.xlsx")
family_summary <- as.data.frame(family_summary)
rownames(family_summary) <- family_summary$Family
family_summary <- family_summary[,-9]
family_summary$K_total <- family_summary$K1+family_summary$K2+family_summary$K3+family_summary$K4+family_summary$K5
family_summary$M_total <- family_summary$M1
family_summary$V_total <- family_summary$V2+family_summary$V3
family_summary$total <- family_summary$M_total+family_summary$K_total+family_summary$V_total

family_summary[["K_total"]][family_summary[["K_total"]] != 0] <- 1
family_summary[["M_total"]][family_summary[["M_total"]] != 0] <- 1
family_summary[["V_total"]][family_summary[["V_total"]] != 0] <- 1

set_K <- family_summary$K_total == 1
set_M <- family_summary$M_total == 1
set_V <- family_summary$V_total == 1

venn(list(K = which(set_K), M = which(set_M), V = which(set_V)))

family_venn <- euler(c("Kimpese city" = 14, Malanga = 6, Viaza = 3, "Kimpese city&Malanga" = 5, "Kimpese city&Viaza" = 5, "Viaza&Malanga" = 1, "Kimpese city&Malanga&Viaza" = 17), counts = TRUE)
png("/results/family_venn.png", width = 800, height = 600)
plot(family_venn, 
     fills = list(fill = c("red", "lightgreen", "blue"), alpha = 0.5), 
     labels = list(font = 2),
     quantities = TRUE)
dev.off()
### Load data for Venn analysis at genus level

genus_summary <- read_excel("Venn_analysis_genus_example_data.xlsx")
genus_summary <- as.data.frame(genus_summary)
rownames(genus_summary) <- genus_summary$genus
genus_summary <- genus_summary[,-9]
genus_summary$K_total <- genus_summary$K1+genus_summary$K2+genus_summary$K3+genus_summary$K4+genus_summary$K5
genus_summary$M_total <- genus_summary$M1
genus_summary$V_total <- genus_summary$V2+genus_summary$V3
genus_summary$total <- genus_summary$M_total+genus_summary$K_total+genus_summary$V_total

genus_summary[["K_total"]][genus_summary[["K_total"]] != 0] <- 1
genus_summary[["M_total"]][genus_summary[["M_total"]] != 0] <- 1
genus_summary[["V_total"]][genus_summary[["V_total"]] != 0] <- 1
genus_summary$shared_genus <- genus_summary$K_total+genus_summary$M_total+genus_summary$V_total

set_K <- genus_summary$K_total == 1
set_M <- genus_summary$M_total == 1
set_V <- genus_summary$V_total == 1

venn(list(K = which(set_K), M = which(set_M), V = which(set_V)))

genus_venn <- euler(c("Kimpese city" = 26, Malanga = 12, Viaza = 5, "Kimpese city&Malanga" = 5, "Kimpese city&Viaza" = 9, "Viaza&Malanga" = 1, "Kimpese city&Malanga&Viaza" = 6), counts = TRUE)
png("/results/genus_venn.png", width = 800, height = 600)
plot(genus_venn, 
     fills = list(fill = c("red", "lightgreen", "blue"), alpha = 0.5), 
     labels = list(font = 2),
     quantities = TRUE)
dev.off()
### Load data for viral relative abundance analysis at family level (Plot top 5 families)
top5_family <- read_excel("relative_abundance_family_example_data.xlsx")
top5_family$Family_reorder  <- with(top5_family, reorder(top5_family$Family, top5_family$RelativeAbundance))
top5_family$label <- paste(top5_family$Family, round(top5_family$RelativeAbundance,digits=2), sep = " ")

family_bar<-ggplot(data=top5_family, aes(x=top5_family$SampleName, y=top5_family$RelativeAbundance,fill=top5_family$Family_reorder)) +
  geom_bar(stat="identity")+ 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_viridis(discrete=TRUE) +
  labs(x="Sample", y="Relative abundance", fill="Family", size=14)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),     
        axis.text.y = element_text(size = 10),      
        axis.title.x = element_text(size = 14),     
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size = 10),         
        legend.title = element_text(size = 12))   

family_bar

### Load data for viral relative abundance analysis at genus level (Plot top 5 genera)
top5_genus <- read_excel("relative_abundance_genus_example_data.xlsx")
top5_genus$Genus_reorder  <- with(top5_genus, reorder(top5_genus$Genus, top5_genus$RelativeAbundance))
top5_genus$label <- paste(top5_genus$Genus, round(top5_genus$RelativeAbundance,digits=2), sep = " ")

genus_bar<-ggplot(data=top5_genus, aes(x=top5_genus$SampleName, y=top5_genus$RelativeAbundance,fill=top5_genus$Genus_reorder)) +
  geom_bar(stat="identity")+ 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_viridis(discrete=TRUE) +
  labs(x="Sample", y="Relative abundance", fill="Genus", size=14)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),     
        axis.text.y = element_text(size = 10),      
        axis.title.x = element_text(size = 14),     
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size = 10),         
        legend.title = element_text(size = 12))   
genus_bar

### PCA analysis based on viral genus annotation

PCAdf <- read_excel("PCA_analysis_genus_example_data.xlsx")
PCAdf <- as.data.frame(PCAdf)
rownames(PCAdf) <- PCAdf[,1]
PCAdf <- PCAdf[,-1]
data_matrix <- as.matrix(PCAdf)
sum(is.na(data_matrix))
out.pca <- prcomp(data_matrix, scale. = TRUE)
attributes(out.pca)
plot(out.pca)
biplot(out.pca)
Fp <- out.pca$x
dim(Fp)

plot(Fp[,1],Fp[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis",main="CHD")#from the plot, we find some outliners with Fp[,1] >50 or <-50
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")

la <- apply(Fp,2,var)
la
fr <- la/sum(la)#calculate the percentage of the single factor can explain
fr
cu <- cumsum(fr)#calculate the sum of the percentage of the factors can explain
cu

collection_site <- factor(c("Kimpese city", "Kimpese city", "Kimpese city","Kimpese city","Kimpese city","Malanga", "Viaza","Viaza"))
pca_data <- as.data.frame(out.pca$x[, 1:2])  # Get the first 2 PCs
pca_data$CollectionSite <- collection_site  
pca_data$Sample <- rownames(pca_data)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = CollectionSite, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA anaysis",
       x = "PC1 (29.8%)", y = "PC2 (24.3%)",
       legend = "Sample collection sites") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()  # Remove the gridlines
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  # Horizontal line at y=0
  geom_vline(xintercept = 0, color = "black",linetype = "dashed") +
  scale_x_continuous(breaks = seq(-6, 10, by = 2)) +  # Set ticks on x-axis at an interval of 2
  scale_y_continuous(breaks = seq(-8,8, by = 2)) +
  geom_text_repel() + 
  scale_color_manual(values = c("Kimpese city" = "red", "Malanga" = "lightgreen", "Viaza" = "blue"))
pca_plot

### Plotting blast and krakenuniq results for huamn and animal viruses

annotation_result <- read_excel("BLAST_krakenuniq_results_example_data.xlsx")
annotation_plot <- ggplot(annotation_result, aes(x=Sample, y=Viruses, color=`Krakenuniq annonation`, fill=`Krakenuniq annonation`, shape=Blast))+ geom_point(size=12.5)+scale_shape_manual(values = c(15,16))+
  labs(x="Mosquito pool", y="", color="Krakenuniq annotation result",fill="Krakenuniq annotation result", shape="Blast result")+
  theme(axis.text.y = element_text(size = 30, face="italic"),
        axis.text.x = element_text(size = 30),# Specify family here
        axis.title = element_text(size = 30),  # Specify family here
        legend.text = element_text(size = 28),  # Specify family here
        legend.title = element_text(size = 30))+
  scale_color_manual(values = c("negative" = "white", "nt and viral databases" = "green", "viral database" = "purple", "nt database"= "blue"))+
  scale_fill_manual(values = c("negative" = "black", "nt and viral databases" = "green", "viral database" = "purple", "nt database"= "blue"))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=27))+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  guides(shape = guide_legend(override.aes = list(fill = "white")))
annotation_plot

### Plotting blood meal analysis results
Blood_meal_annotation <- read_excel("Blood_meal_analysis_example_data.xlsx")
Blood_meal_plot <- ggplot(Blood_meal_annotation, aes(x=Blood_meal_annotation$`Sample name`, y=Genus))+ 
  geom_point(aes(fill=as.factor(Blood_meal_annotation$Detection), color=as.factor(Blood_meal_annotation$Detection)), size=6, shape=15)+
  scale_fill_manual(values = c("0" = "white", "1" = "red"),labels = c("Negative", "Positive")) +
  scale_color_manual(values = c("0" = "black", "1" = "red"),labels = c("Negative", "Positive"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12,angle=15, face="italic"),axis.title = element_text(size = 14),legend.text = element_text(size = 12),legend.title = element_text(size = 13),legend.position="bottom")+
  labs(color = "Detection", fill = "Detection")+
  xlab("Mosquito pool")+
  theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(axis.text.x = element_text(angle = 45,hjust=1))

Blood_meal_plot
