###########################################################################################################################################################################
# This is an R script
#Programmer: Wenqiao He
#Last update: Aug 2025
#Purpose: Visualization of Aedes mosquito viral metagenomic data analysis (after filtering) and phylogenetic analysis
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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
install.packages("ggnewscale")
install.packages("scico")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("treeio")
install.packages("RColorBrewer")
library(RColorBrewer)
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
library(treeio)
library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(scico)

### Set working directory###
setwd("/PATH/")

### Load data for heatmap at family level
family_RPM <- read_excel("family_filtering.xlsx")
family_long <- family_RPM %>%
  pivot_longer(
    cols = -taxName,
    names_to = "Sample",
    values_to = "RPM"
  )%>%
  mutate(RPM_log = log10(RPM + 1))

ggplot(family_long, aes(x = Sample, y = taxName, fill = RPM_log)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.y = element_text(size=12, face = "italic"))


### Load data for heatmap at genus level
genus_RPM <- read_excel("genus_filtering.xlsx")
genus_long <- genus_fltering %>%
  pivot_longer(
    cols = -taxName,
    names_to = "Sample",
    values_to = "RPM"
  )%>%
  mutate(RPM_log = log10(RPM + 1))

ggplot(genus_long, aes(x = Sample, y = taxName, fill = RPM_log)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.y = element_text(size=12, face = "italic"))

### Load data for viral relative abundance analysis at family level 
family_RA <- family_RPM
family_RA[, -1] <- lapply(family_RA[, -1], function(x) as.numeric(x) / 1e6 * 100)
family_RA_long_table <- family_RA %>%
  pivot_longer(
    cols = -taxName,              # all columns except taxName
    names_to = "Sample",          # new column for sample names
    values_to = "Relative abundance"    # new column for the numeric value
  )

family_bar<-ggplot(data=family_RA_long_table, aes(x=Sample, y=family_RA_long_table$`Relative abundance`,fill=taxName)) +
  geom_bar(stat="identity")+ ####geom_text(aes(label = paste0(label, "%")), vjust = 0,angle=30, size=2) +
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_viridis(discrete=TRUE) +
  labs(x="Sample", y="Relative abundance", fill="Family", size=14)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),     
        axis.text.y = element_text(size = 12),      
        axis.title.x = element_text(size = 16),     
        axis.title.y = element_text(size = 16), 
        legend.text = element_text(size = 12),         
        legend.title = element_text(size = 14))   
family_bar

### Load data for viral relative abundance analysis at genus level 
genus_RA <- genus_RPM
genus_RA[, -1] <- lapply(genus_RA[, -1], function(x) as.numeric(x) / 1e6 * 100)

genus_RA_long_table <- genus_RA %>%
  pivot_longer(
    cols = -taxName,              # all columns except taxName
    names_to = "Sample",          # new column for sample names
    values_to = "Relative abundance"    # new column for the numeric value
  )

genus_bar<-ggplot(data=genus_RA_long_table, aes(x=Sample, y=genus_RA_long_table$`Relative abundance`,fill=taxName)) +
  geom_bar(stat="identity")+ ####geom_text(aes(label = paste0(label, "%")), vjust = 0,angle=30, size=2) +
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_viridis(discrete=TRUE) +
  labs(x="Sample", y="Relative abundance", fill="Genus", size=14)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),     
        axis.text.y = element_text(size = 12),      
        axis.title.x = element_text(size = 16),     
        axis.title.y = element_text(size = 16), 
        legend.text = element_text(size = 12),         
        legend.title = element_text(size = 14))   
genus_bar

### PCA analysis based on viral genus annotation
PCA_df <- read_excel("PCA_genus_filtering.xlsx")
PCAdf <- as.data.frame(PCA_df)
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

collection_site <- factor(c("Kiasungua", "Kiasungua", "Kiasungua","Kiasungua","Kiasungua","Malanga", "Viaza","Viaza"))
collection_site_kimpese_city <- factor(c("Kimpese city", "Kimpese city", "Kimpese city","Kimpese city","Kimpese city","Malanga", "Viaza","Viaza"))

pca_data <- as.data.frame(out.pca$x[, 1:2])  # Get the first 2 PCs
pca_data$CollectionSite <- collection_site  
pca_data$CollectionSite_kimpese_city <- collection_site_kimpese_city  
pca_data$Sample <- rownames(pca_data)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = CollectionSite_kimpese_city, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA anaysis",
       x = "PC1 (45.9%)", y = "PC2 (25.1%)", ##if use PCA_analysis_8pools.xlsx,  x = "PC1 (29.5%)", y = "PC2 (23.9%)"
       legend = "Sample collection sites") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        
        # Draw only axis lines at x=0 and y=0 using `geom_hline()` and `geom_vline()`
        panel.grid = element_blank()  # Remove the gridlines
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  # Horizontal line at y=0
  geom_vline(xintercept = 0, color = "black",linetype = "dashed") +
  scale_x_continuous(breaks = seq(-6, 10, by = 2)) +  # Set ticks on x-axis at an interval of 2
  scale_y_continuous(breaks = seq(-8,8, by = 2)) +
  geom_text_repel() + 
  theme(axis.text.x = element_text(size = 12),     
        axis.text.y = element_text(size = 12),      
        axis.title.x = element_text(size = 16),     
        axis.title.y = element_text(size = 16), 
        legend.text = element_text(size = 12),         
        legend.title = element_text(size = 14))  +
  scale_color_manual(name = "Collection Site", values = c("Kimpese city" = "#117733", "Malanga" = "#DDCC77", "Viaza" = "#332288"))

pca_plot
                         
# Top contributing variables to PC1
sort(abs(out.pca$rotation[,1]), decreasing = TRUE)[1:10]
# Top contributing variables to PC2
sort(abs(out.pca$rotation[,2]), decreasing = TRUE)[1:10]


### Plotting BLAST and KrakenUniq results for huamn and animal viruses
blast_result <- read_excel("BLAST_confirmation_filtering.xlsx")

annotation_plot <- ggplot(blast_result, aes(x = Sample, y = Viruses, fill = Annotation, color = Annotation)) + 
  geom_point(size = 12.5, shape = 22) +
  labs(x = "Mosquito pool", y = "Viruses", fill = "Annotation", color = "Annotation") +
  scale_color_manual(values = c(
    "KrakenUniq positive+blast positive" = "#D55E00", 
    "KrakenUniq positive+blast negative" = "#0072B2", 
    "KrakenUniq negative+blast negative" = "#000000"
  )) +
  scale_fill_manual(values = c(
    "KrakenUniq positive+blast positive" = "#D55E00", 
    "KrakenUniq positive+blast negative" = "#0072B2", 
    "KrakenUniq negative+blast negative" = "#FFFFFF"
  )) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 27)) +
  theme_bw() +  # base theme
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  guides(shape = guide_legend(override.aes = list(fill = "white")))

annotation_plot

###Pan-DENV reference selection, download, and phylogenetic tree visualization
##References (>10,000 bp) representing different genotypes were randomly selected from various geographic regions.  
##The file "all_ref.xlsx" contains all DENV reference sequences from the NCBI Virus database, downloaded on June 25, 2025. 
all_ref <- read_excel("all_ref.xlsx")
all_ref <- all_ref %>%
  mutate(Organism_Name = str_replace(Organism_Name, "dengue virus type I", "dengue virus type 1"))%>%
  mutate(Organism_Name = str_replace(Organism_Name, "dengue virus", "Dengue virus"))

result <- all_ref %>%
  group_by(Collection_year, Country, Organism_Name) %>%
  slice_max(order_by = Length, n = 1, with_ties = FALSE) %>%
  ungroup()

countries <- ne_countries(returnclass = "sf") %>%
  sf::st_drop_geometry() %>%
  select(name, region_un, subregion,continent = continent) %>%
  rename(Country = name) 

result_with_continent <- result%>%
  left_join(countries, by = "Country")
                   
##References were selected to ensure representation across diverse geographical regions and serotypes.
random_all_ref <- result_with_continent %>%
  filter(Length >= 10000) %>%                               # Remove reads with length < 10000
  group_by(continent, Organism_Name) %>%                    # Group by region and serotype
  slice_sample(n = 3) %>%                                   # Randomly select 3 reads per group
  ungroup()
                   
##Here is the final list of references included in the Pan-DENV phylogenetic analysis. Some references were removed due to high diversity compared to others of the same serotype. A Zika virus sequence was included as the outgroup.)
pan_denv_refs <- read_excel("DENV_all_serotypes_metadata.xlsx")
pan_denv_published_refs <- subset(pan_denv_refs,pan_denv_refs$This_study=="no") 

##Download references from NCBI 
acc_all <- pan_denv_published_refs$Accession
fasta_seqs <- lapply(acc_all, function(acc) {
  entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
})

combined_fasta <- unlist(fasta_seqs)
writeLines(combined_fasta, "Pan_DENV_all_refs.fasta")

##Load the phylogenetic tree for visualization
Pan_DENV_raxmlng <- read.tree("DENV_all_serotypes_raxml_ng.raxml.support")
Pan_DENV_raxmlng_tree <- ggtree(all_serotypes_raxmlng) %<+% pan_denv_refs +
  geom_tiplab(aes(color = This_study), size = 2.5) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=4.5, y=-1.2, width=0.2, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()
Pan_DENV_raxmlng_tree

continent_data <- as.data.frame(pan_denv_refs[, c("Continent","Accession")])
rownames(continent_data) <- continent_data$Accession
continent_data$Accession <- NULL

serotype_data <- as.data.frame(pan_denv_refs[, c("Organism_Name","Accession")])
rownames(serotype_data) <- serotype_data$Accession
serotype_data$Accession <- NULL

Pan_DENV_merge_1 <- gheatmap(Pan_DENV_raxmlng_tree, continent_data,
                    width = 0.1, offset = 1.15,
                    colnames = TRUE, font.size = 3,
                    colnames_position = "top",
                    colnames_offset_x = -0.15, colnames_offset_y = 0.5,
                    colnames_angle = 45, color = NA, hjust = 0) +
  scale_fill_brewer(palette = "Dark2", na.value = "gray90", name = "Continent")
Pan_DENV_merge_1

Pan_DENV_merge_1_1 <- Pan_DENV_merge_1 + new_scale_fill() 

Pan_DENV_merge_2 <- gheatmap(
  Pan_DENV_merge_1_1, serotype_data,
  width = 0.1, offset = 0.65,
  colnames = TRUE, font.size = 3,
  colnames_position = "top",
  colnames_offset_x = -0.15, colnames_offset_y = 0.5,
  colnames_angle = 45, color = NA, hjust = 0
) +
  scale_fill_brewer(palette = "Paired", name = "Serotype", na.value = "gray90")

Pan_DENV_merge_2

###DENV-2 reference selection, download, and phylogenetic tree visualization
## Add continent information
all_ref_with_continent <- all_ref%>%
  left_join(countries, by = "Country")
DENV2_ref <- all_ref_with_continent %>%
  filter(Length >= 10000) %>% # Remove reads with length < 10000
  subset(Organism_Name == "Dengue virus type 2")

##Select the longest one from each country and year
DENV2_longest <- DENV2_ref %>%
     group_by(Collection_year, Country) %>%
     slice_max(order_by = Length, n = 1, with_ties = FALSE) %>%
     ungroup()

##Select references from different geographical regions
DENV2_random_ref <- DENV2_longest %>%
  group_by( continent) %>%
  slice_sample(n=15) %>%
  ungroup()

##This is the final list of references included in the DENV-2 phylogenetic analysis. Some were excluded due to high genetic divergence. References with genotype information from a published study (PMID: 34188091) were included. A DENV-4 sequence was used as the outgroup.
denv2_refs <- read_excel("DENV2_metadata.xlsx")
denv2_published_refs <- subset(denv2_refs,denv2_refs$This_study=="no") 
##Download references from NCBI 
acc_denv2 <- denv2_published_refs$Accession
fasta_seqs_denv2 <- lapply(acc_denv2, function(acc) {
  entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
})

combined_fasta_denv2 <- unlist(fasta_seqs_denv2)
writeLines(combined_fasta_denv2, "DENV2_all_refs.fasta")

##Load phylogenetic tree for visualization
DENV2_raxmlng_tree <- read.tree("DENV2_raxml_ng.raxml.support")
DENV2_tree <- ggtree(DENV2_raxmlng_tree)%<+% denv2_refs +
  geom_tiplab(aes(color = This_study), size = 2.5) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=2, y=-1, width=0.1, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()

DENV2_tree

DENV2_continent_data <- as.data.frame(denv2_refs[, c("Continent","Accession")])
rownames(DENV2_continent_data) <- DENV2_continent_data$Accession
DENV2_continent_data$Accession <- NULL

DENV2_genotype_data <- as.data.frame(denv2_refs[, c("Genotype","Accession")])
rownames(DENV2_genotype_data) <- DENV2_genotype_data$Accession
DENV2_genotype_data$Accession <- NULL

DENV2_merge_1 <- gheatmap(DENV2_tree, DENV2_continent_data,
                          width = 0.1, offset = 0.35,
                          colnames = TRUE, font.size = 3,
                          colnames_position = "top",
                          colnames_offset_x = -0.15, colnames_offset_y = 0.5,
                          colnames_angle = 45, color = NA, hjust = 0) +
  scale_fill_brewer(
    palette = "Dark2",
    name = "Continent",
    limits = c("Africa", "Asia", "Oceania", "South America", "Europe", "North America"),
    na.value = "gray90"
  )
DENV2_merge_1
DENV2_merge_1_1 <- DENV2_merge_1 + new_scale_fill()

library(RColorBrewer)

DENV2_merge_2 <- gheatmap(
  DENV2_merge_1_1, DENV2_genotype_data,
  width = 0.1, offset = 0.1,
  colnames = TRUE, font.size = 1,
  colnames_position = "top",
  colnames_offset_x = 0.1, colnames_offset_y = 1,
  colnames_angle = 0, color = NA, hjust = 0
) +
  scale_fill_brewer(
    palette = "PiYG",   # or "Dark2"
    name = "Genotype",
    na.value = "gray90"
  )

DENV2_merge_2

###DENV-4 reference selection, download, and phylogenetic tree visualization
## Add continent information
DENV4_ref <- all_ref_with_continent %>%
  filter(Length >= 10000) %>% # Remove reads with length < 10000
  subset(Organism_Name == "Dengue virus type 4")

##Select the longest one from each country and year
DENV4_longest <- DENV4_ref %>%
  group_by(Collection_year, Country) %>%
  slice_max(order_by = Length, n = 1, with_ties = FALSE) %>%
  ungroup()

##Select references from different geographical regions
DENV4_random_ref <- DENV4_longest %>%
  group_by( continent) %>%
  slice_sample(n=20) %>%
  ungroup()

##This is the final list of references included in the DENV-4 phylogenetic analysis. Some were excluded due to high genetic divergence. References with genotype information from a published study (PMID: 37293986) were included. A DENV-2 sequence was used as the outgroup.
denv4_refs <- read_excel("DENV4_metadata.xlsx")
denv4_published_refs <- subset(denv4_refs,denv4_refs$This_study=="no") 
##Download references from NCBI 
acc_denv4 <- denv4_published_refs$Accession
fasta_seqs_denv4 <- lapply(acc_denv4, function(acc) {
  entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
})

combined_fasta_denv4 <- unlist(fasta_seqs_denv4)
writeLines(combined_fasta_denv4, "DENV4_all_refs.fasta")

##Load phylogenetic tree for visualization
DENV4_raxmlng_tree <- read.tree("DENV4_raxml_ng.raxml.support")
DENV4_tree <- ggtree(DENV4_raxmlng_tree)%<+% denv4_refs +
  geom_tiplab(aes(color = This_study), size = 2.5) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=1.75, y=-1, width=0.1, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()

DENV4_tree

DENV4_continent_data <- as.data.frame(denv4_refs[, c("continent","Accession")])
rownames(DENV4_continent_data) <- DENV4_continent_data$Accession
DENV4_continent_data$Accession <- NULL

DENV4_genotype_data <- as.data.frame(denv4_refs[, c("Genotype.x","Accession")])
rownames(DENV4_genotype_data) <- DENV4_genotype_data$Accession
DENV4_genotype_data$Accession <- NULL

DENV4_merge_1 <- gheatmap(DENV4_tree, DENV4_continent_data,
                          width = 0.1, offset = 0.35,
                          colnames = TRUE, font.size = 3,
                          colnames_position = "top",
                          colnames_offset_x = -0.15, colnames_offset_y = 0.5,
                          colnames_angle = 45, color = NA, hjust = 0) +
  scale_fill_brewer(
    palette = "Dark2",
    name = "Continent",
    limits = c("Africa", "Asia", "Oceania", "South America", "North America"),
    na.value = "gray90"
  )
DENV4_merge_1

DENV4_merge_1_1 <- DENV4_merge_1 + new_scale_fill()

DENV4_merge_2 <- gheatmap(
  DENV4_merge_1_1, DENV4_genotype_data,
  width = 0.1, offset = 0.1,
  colnames = TRUE, font.size = 6,
  colnames_position = "top",
  colnames_offset_x = 0.1, colnames_offset_y = 1,
  colnames_angle = 0, color = NA, hjust = 0
) +
  scale_fill_brewer(
    palette = "RdBu",   # or "Dark2"
    name = "Genotype",
    na.value = "gray90"
  )

DENV4_merge_2

###Bat faecal associated dicistrovirus 4 reference (dicistrovirus sequences from differnet hosts) download and phylogenetic tree visualization
bat_dic_refs <- read_excel("bat_dic_metadata.xlsx")
bat_dic_published_refs <- subset(bat_dic_refs,bat_dic_refs$This_study=="no") 

##Download references from NCBI 
acc_bat_dic <- bat_dic_published_refs$Accession
fasta_seqs_bat_dic <- lapply(acc_bat_dic, function(acc) {
  entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
})

combined_fasta_bat_dic <- unlist(fasta_seqs_bat_dic)
writeLines(combined_fasta_bat_dic, "bat_dic_all_refs.fasta")

##Load phylogenetic tree for visualization
bat_dic_metadata <- read_excel("bat_dic_metadata.xlsx")
bat_dic <- read.tree("Bat_dic_raxml_ng.raxml.support")
bat_dic_tree <- ggtree(bat_dic)%<+% bat_dic_metadata +
  geom_tiplab(aes(color = This_study), size = 5,, align=TRUE, linetype="dashed", linesize=0.5, offset=0.01) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=1.75, y=-1, width=0.1, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()
bat_dic_tree

bat_dic_host_data <- as.data.frame(bat_dic_metadata[, c("Host","Accession")])
rownames(bat_dic_host_data) <- bat_dic_host_data$Accession
bat_dic_host_data$Accession <- NULL

bat_dic_tree_merge_1 <- gheatmap(
  bat_dic_tree, bat_dic_host_data,
  width = 0.1, offset = 0.6,
  colnames = TRUE, font.size = 6,
  colnames_position = "top",
  colnames_offset_x = 0.1, colnames_offset_y = 1,
  colnames_angle = 0, color = NA, hjust = 0
) +
  scale_fill_brewer(
    palette = "Set2",   # or "Dark2"
    name = "Host/Source",
    na.value = "gray90"
  )

bat_dic_tree_merge_1

###Human blood-associated dicistrovirus reference (dicistrovirus sequences from different hosts) download and phylogenetic tree visualization
human_dic_refs <- read_excel("human_dic_metadata.xlsx")
human_dic_published_refs <- subset(human_dic_refs,bat_dic_refs$This_study=="no") 
##Download references from NCBI 
acc_human_dic <- human_dic_published_refs$Accession
fasta_seqs_human_dic <- lapply(acc_human_dic, function(acc) {
  entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
})

combined_fasta_human_dic <- unlist(fasta_seqs_human_dic)
writeLines(combined_fasta_human_dic, "human_dic_all_refs.fasta")

##Load metadata for phylogenetic tree visualization. Three separate trees were generated for Human blood-associated dicistrovirus, as reads mapped to different regions of the genome.
human_dic_metadata <- read_excel("human_dic_metadate.xlsx")
human_dic_K5_1 <- read.tree("Human_dic_K5_1_raxml_ng.raxml.support")
human_dic_K5_1_tree <- ggtree(human_dic_K5_1)%<+% human_dic_metadata+
  geom_tiplab(aes(color = This_study), size = 5, align=TRUE, linetype="dashed", linesize=0.5, offset=0.01) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=1.25, y=-0.5, width=0.1, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()
human_dic_K5_1_tree

human_dic_host_data <- as.data.frame(human_dic_metadata[, c("Host","Accession")])
rownames(human_dic_host_data) <- human_dic_host_data$Accession
human_dic_host_data$Accession <- NULL

human_dic_K5_1_tree_merge_1 <- gheatmap(
  human_dic_K5_1_tree, human_dic_host_data,
  width = 0.1, offset = 0.6,
  colnames = TRUE, font.size = 6,
  colnames_position = "top",
  colnames_offset_x = 0.1, colnames_offset_y = 1,
  colnames_angle = 0, color = NA, hjust = 0
) +
  scale_fill_brewer(
    palette = "Set2",   # or "Dark2"
    name = "Host/Source",
    na.value = "gray90"
  )
human_dic_K5_1_tree_merge_1

human_dic_K5_2 <- read.tree("Human_dic_K5_2_raxml_ng.raxml.support")
human_dic_K5_2_tree <- ggtree(human_dic_K5_2)%<+% human_dic_metadata+
  geom_tiplab(aes(color = This_study), size = 5, align=TRUE, linetype="dashed", linesize=0.5, offset=0.01) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=1.25, y=-0.5, width=0.1, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()
human_dic_K5_2_tree

human_dic_K5_2_tree_merge_1 <- gheatmap(
  human_dic_K5_2_tree, human_dic_host_data,
  width = 0.1, offset = 0.6,
  colnames = TRUE, font.size = 6,
  colnames_position = "top",
  colnames_offset_x = 0.1, colnames_offset_y = 1,
  colnames_angle = 0, color = NA, hjust = 0
) +
  scale_fill_brewer(
    palette = "Set2",   # or "Dark2"
    name = "Host/Source",
    na.value = "gray90"
  )
human_dic_K5_2_tree_merge_1

human_dic_K4 <- read.tree("Human_dic_K4_raxml_ng.raxml.support")
human_dic_K4_tree <- ggtree(human_dic_K4)%<+% human_dic_metadata+
  geom_tiplab(aes(color = This_study), size = 5, align=TRUE, linetype="dashed", linesize=0.5, offset=0.01) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=1.25, y=-0.5, width=0.1, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()
human_dic_K4_tree

human_dic_K4_tree_merge_1 <- gheatmap(
  human_dic_K4_tree, human_dic_host_data,
  width = 0.1, offset = 0.6,
  colnames = TRUE, font.size = 6,
  colnames_position = "top",
  colnames_offset_x = 0.1, colnames_offset_y = 1,
  colnames_angle = 0, color = NA, hjust = 0
) +
  scale_fill_brewer(
    palette = "Set2",   # or "Dark2"
    name = "Host/Source",
    na.value = "gray90"
  )
human_dic_K4_tree_merge_1


###Human pegivirus reference (pegivirus sequences from different hosts) download and phylogenetic tree visualization
pegivirus_refs <- read_excel("Pegivirus_metadata.xlsx")
pegivirus_published_refs <- subset(pegivirus_refs,pegivirus_refs$This_study=="no") 

##Download references from NCBI 
acc_pegivirus <- pegivirus_published_refs$Accession
fasta_seqs_pegivirus <- lapply(acc_pegivirus, function(acc) {
  entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
})

combined_fasta_pegivirus <- unlist(fasta_seqs_pegivirus)
writeLines(combined_fasta_pegivirus, "pegivirus_all_refs.fasta")

##Load metadata for phylogenetic tree visualization
pegivirus_metadata <- read_excel("pegivirus_metadata.xlsx")
pegivirus_host_data <- as.data.frame(pegivirus_metadata[, c("Host","Accession")])
rownames(pegivirus_host_data) <- pegivirus_host_data$Accession
pegivirus_host_data$Accession <- NULL

Pegivirus <- read.tree("Pegivrus_raxml_ng.raxml.support")
Pegivirus_tree <- ggtree(Pegivirus)%<+% pegivirus_metadata +
  geom_tiplab(aes(color = This_study), size = 5, align=TRUE, linetype="dashed", linesize=0.5, offset=0.01) +
  scale_color_manual(
    values = c("yes" = "red", "no" = "black")
  ) +
  geom_treescale(x=1.75, y=-0.5, width=0.1, fontsize = 3)+
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  theme_tree()
Pegivirus_tree

Pegivirus_tree_merge_1 <- gheatmap(
  Pegivirus_tree, pegivirus_host_data,
  width = 0.1, offset = 0.6,
  colnames = TRUE, font.size = 6,
  colnames_position = "top",
  colnames_offset_x = 0.1, colnames_offset_y = 1,
  colnames_angle = 0, color = NA, hjust = 0
) +
  scale_fill_brewer(
    palette = "Set2",   # or "Dark2"
    name = "Host",
    na.value = "gray90"
  )
Pegivirus_tree_merge_1



