###########################################################################################################################################################################
# This is an R script
#Programmer: Wenqiao He
#Last update: December 2024
#Purpose: create maps of sample collection sites within the Democratic Republic of the Congo (DRC)
###########################################################################################################################################################################

install.packages("sf", type = "source", 
                 configure.args = "--with-gdal-config=/usr/bin/gdal-config --with-proj-lib=/usr/lib")
install.packages("tmap")
install.packages("sp")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("ggplot2")

library(ggplot2)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(tmap)

###Set working directory###
setwd("/Path/")

# Load world data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Get the data of African and the Democratic Republic of the Congo (DRC)
africa <- world[world$continent == "Africa", ]  # Extract only African countries
drc <- africa[africa$admin == "Democratic Republic of the Congo", ]  # Filter for DRC

# Plot Africa and highlight DRC
ggplot(data = africa) +
  geom_sf(fill = "lightgray", color = "black") +  # Base map of Africa
  geom_sf(data = drc, fill = "red", color = "black") +  # Highlight DRC
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(filename = "results/Africa_map.png", width = 8, height = 6)

# Get the data of DRC, load administrative boundaries (.shp file) from https://gadm.org/
congo <- ne_countries(scale = "medium", country = "Democratic Republic of the Congo", returnclass = "sf")
city_boundaries <- st_read("gadm41_COD_1.shp")

# Base map of DRC with city borders
Congo_map <- ggplot(data = congo) +
  geom_sf(fill = "lightgrey", color = "white")+
  geom_sf(data = city_boundaries, fill = NA, color = "white", size = 0.5)

# Highlight Kongo-central
KongoCentral <- subset(city_boundaries, city_boundaries$NAME_1 == "Kongo-Central")

Congo_map <- Congo_map+
  geom_sf(data = KongoCentral, fill = "red", color = "white")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(Congo_map)

ggsave(filename = "results/DRC_map.png", plot = Congo_map, width = 8, height = 6)

## Get Kongo-Central base mape and get boundaries from https://data.humdata.org/ (use RDC_Aires de santé.shp file), highlight Kimpese

KongoCentral_map <- ggplot(data = KongoCentral) +
  geom_sf(fill = "lightgrey", color = "white")+
  geom_sf(data = KongoCentral, fill = NA, color = "white", size = 0.5)

boundaries <- st_read("RDC_Aires de santé.shp")

# Get Kimpese data and highlight Kimpese

Kimpese <- boundaries[boundaries$ZS == "Kimpese", ]
Kimpese_map <- KongoCentral_map+
  geom_sf(data = Kimpese, fill = "red", color = "white")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(Kimpese_map)

ggsave(filename = "results/KCP_map.png", plot = Kimpese_map, width = 8, height = 6)

# Filter the data for Ceco, Malanga, and Viaza
Ceco <- Kimpese[Kimpese$AS_ == "CECO",]
Malanga <- Kimpese[Kimpese$AS_ == "MALANGA", ]
Viaza <- Kimpese[Kimpese$AS_ == "VIAZA", ]


# Highlight Ceco, Malanga, and Viaza

Kimpese_base <- tm_shape(Kimpese) +
  tm_borders() +
  tm_layout(legend.position = c("left", "bottom"),scale = 0.8) +
  tm_scale_bar(position = c("left", "bottom"), width = 0.2, text.size = 0.7) 

Ceco_highlight <- tm_shape(Ceco) +
  tm_fill("red", alpha = 0.3) +  # Example: Highlight color and transparency
  tm_borders(col = "black") 

Viaza_highlight <- tm_shape(Viaza) +
  tm_fill("blue", alpha = 0.3) +  # Example: Highlight color and transparency
  tm_borders(col = "black") 

Malanga_highlight <- tm_shape(Malanga) +
  tm_fill("lightgreen", alpha = 0.3) +  # Example: Highlight color and transparency
  tm_borders(col = "black")

# Show Kimpese city with a red dot, search the longitude and latitude of Kimpese City in Google
longitude <- 14.4288
latitude <- -5.5618

Kimpese_city <- data.frame(
  name = "Kimpese City",
  lon = longitude,
  lat = latitude
)

point_sf <- st_as_sf(Kimpese_city, coords = c("lon", "lat"), crs = 4326)

Kimpese_map <- Kimpese_base+Ceco_highlight+Viaza_highlight+Malanga_highlight+tm_shape(point_sf) +
  tm_dots(col = "red", size = 1) 

# Display the map
print(Kimpese_map)

tmap_save(Kimpese_map, filename = "results/Kimpese_map.png", width = 800, height = 600)
