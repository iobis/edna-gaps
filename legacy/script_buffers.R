library(dplyr)
library(sf)
library(glue)
library(robis)
library(stringi)
library(mapview)

shapefile <- "https://github.com/iobis/mwhs-shapes/raw/master/output/marine_world_heritage.gpkg"

fish_classes <- c("Actinopteri", "Cladistii", "Coelacanthi", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti", "Teleostei")
turtle_orders <- c("Testudines")
mammal_classes <- c("Mammalia")

# read shapefile

sf_use_s2(FALSE)

centroids <- read_sf(shapefile) %>%
  mutate(name_simplified = gsub("_+", "_", gsub("[^[:alnum:]]", "_", tolower(stri_trans_general(name, "latin-ascii"))))) %>%
  group_by(name_simplified) %>%
  summarize() %>%
  st_centroid() 

buffers_terra <- terra::buffer(terra::vect(sites_centroids), 1000000)
# plot(buffers)
# terra::writeVector(buffers, filename = "sites.shp", overwrite = TRUE)
buffers_sf <- sf::st_as_sf(buffers_terra)
# mapview(buffers_sf, legend = FALSE)
buffers <- buffers_sf %>%
  mutate(wkt = st_as_text(geometry)) %>%
  st_drop_geometry()

