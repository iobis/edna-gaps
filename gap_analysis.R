library(h3)
library(mapview)
library(arrow)
library(DBI)
library(sf)
library(terra)
library(dplyr)
library(stringr)
library(furrr)
library(tidyr)
library(readr)

# TODO: check if abundant species are truly absent from reference db

# config

h3_res <- 3
h3_ring <- 3
export_file <- "~/Desktop/temp/obis_20231025.parquet"
fish_classes <- c("Actinopteri", "Cladistii", "Coelacanthi", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti", "Teleostei")
turtle_orders <- c("Testudines")
mammal_classes <- c("Mammalia")
reference_databases <- list(
  "12s_mimammal" = "reference_databases/12S/202311/12S_mammal_ncbi_1_50000_pcr_pga_taxa.tsv",
  "12s_mifish" = "reference_databases/12S/202311/12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa.tsv",
  "12s_teleo" = "reference_databases/12S/202311/12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa.tsv",
  "coi" = "reference_databases/COI_ncbi/COI_ncbi_1_50000_pcr_pga_taxon.tsv",
  "16s" = "reference_databases/16S/202311/16S_ncbi_euk_1_50000_pga_taxa.tsv"
)

# species list by H3 cell

row_to_geo <- function(row) {
  geo_to_h3(c(row$decimalLatitude, row$decimalLongitude), h3_res)
}

reader <- open_dataset(export_file) %>%
  select(decimalLongitude, decimalLatitude, class, order, species, speciesid) %>%
  as_record_batch_reader()

if (file.exists("temp.db")) {
  file.remove("temp.db")
}
con <- dbConnect(RSQLite::SQLite(), "temp.db")

while (TRUE) {
  batch <- reader$read_next_batch()
  if (is.null(batch)) {
    break
  }
  batch %>%
    as.data.frame() %>%
    mutate(h3 = row_to_geo(.)) %>%
    group_by(h3, class, order, species, speciesid) %>%
    summarize(records = n()) %>%
    dbWriteTable(con, "occurrence", ., append = TRUE)
}

species_by_h3 <- dbGetQuery(con, "select h3, class, `order`, species, speciesid, sum(records) as records from occurrence group by h3, class, `order`, species, speciesid") %>%
  filter(!is.na(speciesid) & !is.na(class) & !is.na(order)) %>%
  mutate(
    group = case_when(
      class %in% fish_classes ~ "fish",
      order %in% turtle_orders ~ "turtle",
      class %in% mammal_classes ~ "mammal"
    )
  ) %>%
  filter(!is.na(group) & species != "Homo sapiens") %>%
  select(-class, -order)

# add presence in reference database

reference_species <- read_csv("reference_databases/reference_species.csv.gz") %>%
  filter(!is.na(speciesid)) %>%
  distinct(speciesid, marker) %>%
  mutate(has_marker = TRUE) %>%
  pivot_wider(id_cols = speciesid, names_from = marker, values_from = has_marker, names_prefix = "marker_") %>%
  mutate_at(vars(starts_with("marker")), ~replace_na(., FALSE)) %>%
  mutate(speciesid = as.numeric(speciesid))

species_by_h3 <- species_by_h3 %>%
  left_join(reference_species, by = "speciesid") %>%
  mutate_at(vars(starts_with("marker")), ~replace_na(., FALSE))

write.csv(species_by_h3, file = gzfile("species_lists/species_lists.csv.gz"), row.names = FALSE, na = "")

# SECTION BELOW IS EXPERIMENTAL

# H3 neighborhood approach demo

hex <- geo_to_h3(c(51, 2), h3_res)
h3_to_geo_boundary_sf(hex) %>% mapview()

ring <- k_ring(hex, h3_ring)
ring <- ring[ring != "0"]
ring_sf <- h3_to_geo_boundary_sf(ring) %>%
  st_wrap_dateline()
ring_sf %>% mapview(legend = FALSE)

# species list by H3 neighborhood

h3_stats <- function(hex) {
  ring <- k_ring(hex, h3_ring)
  ring <- ring[ring != "0"]
  ring_species <- species_by_h3 %>%
    filter(h3 %in% ring) %>%
    select(-h3) %>%
    distinct()
  ring_stats <- ring_species %>%
    select(-speciesid) %>%
    group_by(group) %>%
    summarize(
      species = n(),
      across(everything(), \(x) sum(x, na.rm = TRUE))
    ) %>%
    mutate_at(vars(starts_with("marker_")), ~./species) %>%
    pivot_wider(names_from = group, values_from = c(species, starts_with("marker_")), names_glue = "{group}_{.value}")
  return(ring_stats)
}

h3_indexes <- function(h3_res) {
  pol <- st_sfc(st_polygon(list(rbind(c(0, -90), c(0, 90), c(180, 90), c(180, -90), c(0, -90)))), crs = 4326)
  pol2 <- st_sfc(st_polygon(list(rbind(c(-180, -90), c(-180, 90), c(0, 90), c(0, -90), c(-180, -90)))), crs = 4326)
  indexes <- h3jsr::polygon_to_cells(pol, res = h3_res, simple = TRUE)
  indexes2 <- h3jsr::polygon_to_cells(pol2, res = h3_res, simple = TRUE)
  unique(c(indexes[[1]], indexes2[[1]]))
}

all_indexes <- h3_indexes(h3_res)
stats_list <- purrr::map(all_indexes, h3_stats, .progress = TRUE)
names(stats_list) <- all_indexes
stats <- bind_rows(stats_list, .id = "h3") %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

# stats

sf_use_s2(FALSE)

stats_sf <- h3_to_geo_boundary_sf(stats$h3) %>%
  left_join(stats, by = c("h3_index" = "h3")) %>%
  st_wrap_dateline() %>%
  st_make_valid() %>%
  st_difference(land)

# visualize

mapview(stats_sf, zcol = "fish_marker_coi")

