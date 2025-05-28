filter_traits <- function(grid_sp, try_wide, sp_cols) {
  library(tidyr)
  library(dplyr)
  library(tibble)

  # Step 1: Transpose selected columns in grid_sp to long format
  grid_long <- grid_sp %>%
    pivot_longer(
      cols = all_of(sp_cols),   # Use all_of() for column selection
      names_to = "species",     # New column for species names
      values_to = "value"       # New column for values
    ) %>%
    filter(value > 0)           # Keep only non-zero values

  # Step 2: Ensure 'species' column exists in try_wide
  try_wide <- try_wide %>%
    rownames_to_column("species")  # Convert rownames to species column if needed

  # Step 3: Match species in grid_long with those in try_wide
  sp_traits <- grid_long %>%
    distinct(species) %>%
    inner_join(try_wide, by = "species")  # Match species with traits

  # Step 4: Filter grid_long for matched species
  matched_species <- sp_traits$species

  site_sp <- grid_long %>%
    filter(species %in% matched_species) %>%       # Keep only matched species
    pivot_wider(
      names_from = species,                        # Reshape species to columns
      values_from = value,                         # Values for species counts
      values_fill = 0                              # Fill missing values with 0
    )

  # Step 5: Create site_xy table with filtered grid_id, x, y
  site_xy <- grid_sp %>%
    filter(grid_id %in% site_sp$grid_id) %>%       # Match filtered grid IDs
    select(grid_id, x, y)                          # Select grid_id, x, y

  # Return list of 3 tables
  return(list(sp_traits = sp_traits, site_sp = site_sp, site_xy = site_xy))
}
