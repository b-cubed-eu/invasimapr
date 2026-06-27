# Prepare inputs (assemble and align core tables)

Thin wrapper around
[assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)
that standardises and stores assembled inputs in a
[new_invasimapr_fit](https://b-cubed-eu.github.io/invasimapr/reference/utils_internal.md)
object for downstream modelling. Use this function to run the
input-assembly step once and pass a single structured container through
subsequent pipelines.

## Usage

``` r
prepare_inputs(..., make_plots = FALSE)
```

## Arguments

- ...:

  Arguments passed on to
  [`assemble_matrices`](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)

  `long_df`

  :   Optional long-format data frame with at least the columns `site`,
      `x`, `y`, `species`, and `count`.

  `site_df`

  :   Optional data frame with columns `site`, `x`, and `y`, used when
      `long_df` is `NULL`.

  `env_df`

  :   Optional site by environment numeric data frame or matrix. Row
      names must correspond to site identifiers.

  `comm_res`

  :   Optional site by resident numeric matrix or data frame (wide), or
      a long-format table with columns `site`, `species`, and `count`.

  `traits_res`

  :   Optional resident by trait data frame with row names corresponding
      to species identifiers.

  `comm_long`

  :   Character string indicating how to interpret a separately supplied
      `comm_res`; one of `"auto"`, `"long"`, or `"wide"`.

  `site_col,x_col,y_col,species_col,count_col`

  :   Column names used when building from `long_df`.

  `env_cols`

  :   Character vector of environment column names to extract.

  `env_prefix`

  :   Regular expression used to select environment columns.

  `trait_cols`

  :   Character vector of trait column names to extract.

  `trait_prefix`

  :   Regular expression used to select trait columns.

  `drop_empty_sites`

  :   Logical. If `TRUE`, drop sites with zero total abundance.

  `drop_empty_species`

  :   Logical. If `TRUE`, drop resident species with zero total
      abundance across all sites.

  `return_diversity`

  :   Logical. If `TRUE`, compute and return site-level diversity
      summaries.

- make_plots:

  Logical; if `TRUE`, propagate `make_plots` to
  [assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)
  so that diagnostic plots are produced during assembly. Defaults to
  `FALSE`.

## Value

An object of class `invasimapr_fit` with components:

- inputs:

  The full list returned by
  [assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md).

- meta:

  A list with elements `n_sites`, `n_residents`, and `n_traits`.

## Details

**What this does**

1.  Calls
    [assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)
    to build core site-by-species and site-by-trait matrices, perform
    consistency checks, and harmonise keys and factor levels.

2.  Wraps the returned list into a lightweight S3 container of class
    `invasimapr_fit`, with a dedicated `inputs` slot and a small `meta`
    block containing basic counts used by later steps.

**Object layout**

    invasimapr_fit
      inputs : list   (output from assemble_matrices)
      meta   : list   (n_sites, n_residents, n_traits)

**Notes**

- No mutation of the assembled inputs occurs here; the function only
  packages and annotates them.

- To inspect the assembled data, use `fit$inputs` on the returned
  object.

## See also

[assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md),
[new_invasimapr_fit](https://b-cubed-eu.github.io/invasimapr/reference/utils_internal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- prepare_inputs(
  sites      = site_df,
  residents  = resident_df,
  invaders   = invader_df,
  traits     = trait_df,
  make_plots = TRUE
)

str(fit, 1)
str(fit$inputs, 1)
fit$meta
} # }
```
