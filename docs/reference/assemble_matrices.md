# Assemble community matrices for invasion-fitness workflows

`assemble_matrices()` creates the standard inputs used throughout the
invasion-fitness pipeline from either a single long-format table or a
set of pre-assembled tables. It returns aligned objects including:

- `site_df`: site by x,y coordinates

- `env_df`: site by environment numeric matrix

- `comm_res`: site by resident abundance matrix

- `pa_res`: site by resident presence-absence matrix

- `traits_res`: resident by trait data frame (mixed types allowed)

## Usage

``` r
assemble_matrices(
  long_df = NULL,
  site_df = NULL,
  env_df = NULL,
  comm_res = NULL,
  traits_res = NULL,
  comm_long = c("auto", "long", "wide"),
  site_col = "site",
  x_col = "x",
  y_col = "y",
  species_col = "species",
  count_col = "count",
  env_cols = NULL,
  env_prefix = "^env",
  trait_cols = NULL,
  trait_prefix = "^trait",
  drop_empty_sites = TRUE,
  drop_empty_species = TRUE,
  return_diversity = TRUE,
  make_plots = FALSE
)
```

## Arguments

- long_df:

  Optional long-format data frame with at least the columns `site`, `x`,
  `y`, `species`, and `count`.

- site_df:

  Optional data frame with columns `site`, `x`, and `y`, used when
  `long_df` is `NULL`.

- env_df:

  Optional site by environment numeric data frame or matrix. Row names
  must correspond to site identifiers.

- comm_res:

  Optional site by resident numeric matrix or data frame (wide), or a
  long-format table with columns `site`, `species`, and `count`.

- traits_res:

  Optional resident by trait data frame with row names corresponding to
  species identifiers.

- comm_long:

  Character string indicating how to interpret a separately supplied
  `comm_res`; one of `"auto"`, `"long"`, or `"wide"`.

- site_col, x_col, y_col, species_col, count_col:

  Column names used when building from `long_df`.

- env_cols:

  Character vector of environment column names to extract.

- env_prefix:

  Regular expression used to select environment columns.

- trait_cols:

  Character vector of trait column names to extract.

- trait_prefix:

  Regular expression used to select trait columns.

- drop_empty_sites:

  Logical. If `TRUE`, drop sites with zero total abundance.

- drop_empty_species:

  Logical. If `TRUE`, drop resident species with zero total abundance
  across all sites.

- return_diversity:

  Logical. If `TRUE`, compute and return site-level diversity summaries.

- make_plots:

  Logical. If `TRUE`, generate quick diagnostic plots.

## Value

A named list of aligned community objects with components:

- site_df:

  Site-by-coordinate data frame (rows = sites).

- env_df:

  Site-by-environment numeric data frame.

- comm_res:

  Site-by-resident abundance matrix.

- pa_res:

  Site-by-resident presence-absence matrix.

- traits_res:

  Resident-by-trait data frame (mixed types allowed).

- diversity:

  Optional site-level diversity summaries (present when
  `return_diversity = TRUE`).

- sites, residents:

  Character vectors of aligned site and resident identifiers.

- n_sites, n_env, n_residents, n_traits:

  Integer dimension counts.

- plots:

  List of diagnostic plots (empty unless `make_plots = TRUE`).

## Examples

``` r
if (FALSE) { # \dontrun{
site_env_spp <- read.csv(
  system.file("extdata", "site_env_spp_simulated.csv.gz", package = "invasimapr")
)
names(site_env_spp)[names(site_env_spp) == "site_id"] <- "site"
core <- assemble_matrices(long_df = site_env_spp)
str(core, 1)
} # }
```
