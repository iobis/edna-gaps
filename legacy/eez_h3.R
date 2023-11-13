library(arrow)
library(h3)
library(mapview)

eez <- read_parquet("eez_h3_res5.parquet")
eez_3 <- unique(unlist(purrr::map(eez$h3_index, h3_to_parent, 3)))
h3_to_geo_boundary_sf(eez_3) %>% mapview()
