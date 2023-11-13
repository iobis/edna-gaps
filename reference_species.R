library(h3)
library(mapview)
library(arrow)
library(DBI)
library(sf)
library(terra)
library(dplyr)
library(stringr)
library(furrr)

# config
# TODO: check if this taxonomy works as expected

fish_classes <- c("Actinopteri", "Cladistii", "Coelacanthi", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti", "Teleostei", "Chondrichthyes")
turtle_orders <- c("Testudines")
mammal_classes <- c("Mammalia")

reference_databases <- list(
  "12s_mimammal" = "reference_databases/12S/202311/12S_mammal_ncbi_1_50000_pcr_pga_taxa.tsv",
  "12s_mifish" = "reference_databases/12S/202311/12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa.tsv",
  "12s_teleo" = "reference_databases/12S/202311/12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa.tsv",
  "coi" = "reference_databases/COI_ncbi/COI_ncbi_1_50000_pcr_pga_taxon.tsv",
  "16s" = "reference_databases/16S/202311/16S_ncbi_euk_1_50000_pga_taxa.tsv"
)

# reference databases

get_marker_species_list <- function(marker) {
  ncol <- max(count.fields(reference_databases[[marker]], sep = "\t"), na.rm = TRUE)
  reference_database <- read.table(reference_databases[[marker]], sep = "\t", header = FALSE, col.names = 1:ncol, fill = TRUE, quote = "\"", na.strings = c("", "nan"))[,1:9] %>%
    setNames(c("seqid", "taxid", "kingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
    filter(species != "Homo_sapiens") %>%
    distinct(taxid, kingdom, phylum, class, order, family, genus, species) %>%
    mutate(species = str_replace_all(species, "_", " ")) %>%
    mutate(
      group = case_when(
        class %in% fish_classes ~ "fish",
        order %in% turtle_orders ~ "turtle",
        class %in% mammal_classes ~ "mammal"
      )
    ) %>%
    mutate(marker = marker)
  return(reference_database)
}

full_species_lists <- purrr::map(names(reference_databases), get_marker_species_list) %>%
  bind_rows()

species_lists <-  full_species_lists %>%
  filter(!is.na(group))

# name matching

species_names <- species_lists %>%
  filter(
    !str_detect(species, " sp.") &
    !str_detect(species, " aff.") &
    !str_detect(species, " cf.")
  ) %>%
  distinct(species) %>%
  pull(species)

name_batches <- split(species_names, as.integer((seq_along(species_names) - 1) / 50))
plan(multisession, workers = 10)
matched_batches <- future_map(name_batches, function(.x) {
  res <- worrms::wm_records_names(.x, marine_only = FALSE)
  names(res) <- .x
  return(bind_rows(res, .id = "input"))
})
matched_species <- bind_rows(matched_batches) %>%
  select(species = input, speciesid = valid_AphiaID)

# output

reference_species <- full_species_lists %>%
  left_join(matched_species, by = "species")

write.csv(reference_species, file = gzfile("reference_databases/reference_species.csv.gz"), row.names = FALSE, na = "")
