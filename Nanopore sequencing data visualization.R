###########################################################################################################################################################################
# This is an R script
#Programmer: Wenqiao He
#Last update: July 2025
#Purpose: Visualization of Aedes mosquito nanopore targeted sequencing data analysis
###########################################################################################################################################################################

library(writexl)
library(gplots)
library(dplyr)
library(ggplot2)

### Plotting blood meal analysis results
Blood_meal_annotation <- read_excel("/example data/Blood_meal_analysis_example_data.xlsx")

Blood_meal_annotation_plot <- ggplot(Blood_meal_annotation, aes(x = Blood_meal_annotation$`Sample name`, y = Genus)) + 
  geom_point(
    aes(
      fill = factor(Detection), 
      color = factor(Detection)
    ), 
    size = 12, shape = 22  # shape 22 = filled square
  ) +
  scale_color_manual(values = c(
    "1" = "red", 
    "0" = "black"
  )) +
  scale_fill_manual(values = c(
    "1" = "red", 
    "0" = "white"
  )) +
  labs(
    x = "Mosquito pool", 
    color = "Detection", 
    fill = "Detection"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

Blood_meal_annotation_plot

