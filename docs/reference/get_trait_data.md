# Scrape and Analyze Wikipedia & Trait Data for a Species

Given a binomial species name, this function retrieves optional metadata
from Wikipedia (taxonomic summary, taxonomy, image, and a color palette)
and joins relevant plant/trait data from a TRY-style or user-provided
trait table. Fuzzy matching is used for both TRY and local tables to
handle minor spelling or naming mismatches.

## Usage

``` r
get_trait_data(
  species,
  remove_bg = FALSE,
  do_palette = TRUE,
  do_taxonomy = TRUE,
  do_summary = TRUE,
  do_image = TRUE,
  bg_thresh = 80,
  green_delta = 20,
  n_palette = 5,
  preview = FALSE,
  save_folder = NULL,
  use_try = FALSE,
  try_data = NULL,
  trait_species_col = "AccSpeciesName",
  local_trait_df = NULL,
  local_species_col = "species",
  max_dist = 1
)
```

## Arguments

- species:

  Character. Species name (binomial, e.g. `"Acacia karroo"`).

- remove_bg:

  Logical. If `TRUE`, call remove.bg via `remove_bg_and_save()` to
  remove the background from the Wikipedia infobox image and use the
  processed PNG for preview/palette. Default: `FALSE`.

- do_palette, do_taxonomy, do_summary, do_image:

  Logical. Control which metadata to scrape. All default to `TRUE`.

- bg_thresh:

  Integer. Deprecated/ignored when `remove_bg = TRUE`. Kept for
  signature compatibility. Default: `80`.

- green_delta:

  Integer. Deprecated/ignored when `remove_bg = TRUE`. Kept for
  signature compatibility. Default: `20`.

- n_palette:

  Integer. Number of colors to extract for the palette. Default: `5`.

- preview:

  Logical. Print the processed image (magick) in the console. Default:
  `FALSE`.

- save_folder:

  Character or `NULL`. If non-`NULL`, write the PNG used for
  palette/preview into this folder. When `remove_bg = TRUE`, the
  background-removed PNG is written here; otherwise the original image
  is written. If `NULL` and `remove_bg = TRUE`, a temporary directory is
  used.

- use_try:

  Logical. If `TRUE`, join plant traits using a TRY-format
  database/table. Default: `FALSE`.

- try_data:

  Character (path) or `data.frame`. Path to a TRY file, or a data frame
  containing TRY-style trait data (must include `trait_species_col`,
  `TraitName`, and `OrigValueStr`).

- trait_species_col:

  Name of the species column in the TRY trait table. Default:
  `"AccSpeciesName"`.

- local_trait_df:

  Optional. `data.frame` of local trait data (any wide table). All
  columns except the species column are returned.

- local_species_col:

  Name of the species column in the local trait table. Default:
  `"species"`.

- max_dist:

  Numeric. Maximum distance for fuzzy matching (Jaro–Winkler via
  [`fuzzyjoin::stringdist_left_join`](https://rdrr.io/pkg/fuzzyjoin/man/stringdist_join.html)).
  Default: `1`.

## Value

A tibble (one row) with columns: `species`, optional metadata
(`summary`, taxonomy ranks, `img_url`, `palette`, `image_file`), and all
available trait columns found via TRY/local joins. `image_file` contains
the normalized path to the PNG used for palette/preview (or `NA` if
none).

## Details

When `remove_bg = TRUE`, the infobox image background is removed using
the remove.bg API via an internal helper (`remove_bg_and_save()`), the
resulting PNG is re-read with **magick**, and the palette is extracted
from that processed image. Set the environment variable
`REMOVE_BG_API_KEY` to a valid remove.bg API key before calling.

- **Wikipedia:** summary via REST API; taxonomy parsed from the infobox.

- **Image:** first infobox image is used; when `remove_bg = TRUE` the
  function calls the remove.bg API. Set
  `Sys.setenv(REMOVE_BG_API_KEY = "…")`.

- **Palette:** simple k-means on non-transparent pixels of the
  (processed) PNG.

- **Traits (TRY):** wide table produced from `TraitName` and numeric
  `OrigValueStr`.

- **Traits (local):** returns all columns except the species column.

- **Dependencies:** `dplyr`, `purrr`, `tibble`, `fuzzyjoin`, `rvest`,
  `httr`, `stringr`, `jsonlite`, `magick`, `abind`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using TRY table
get_trait_data("Acacia karroo",
               use_try = TRUE,
               try_data = try_traits,
               trait_species_col = "SpeciesName")

# Using a local trait table
get_trait_data("Acraea horta", local_trait_df = traits, local_species_col = "species")

# Metadata only, with background removal and saving to a folder
Sys.setenv(REMOVE_BG_API_KEY = "<your-removebg-key>")
get_trait_data("Acacia karroo", use_try = FALSE, remove_bg = TRUE, save_folder = "out/")
} # }
```
