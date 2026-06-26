# Simulate hypothetical invader trait profiles from a resident trait pool

Generate \\n\_{\mathrm{inv}}\\ rows of trait values for hypothetical
invaders by resampling the empirical distribution of resident traits. By
default, traits are sampled independently by column, creating novel
trait combinations. Alternatively, entire rows can be resampled to
preserve the covariance structure of the resident trait space.

Row names of the returned object are set to the invader identifiers.

## Usage

``` r
simulate_invaders(
  resident_traits,
  n_inv = 10,
  species_col = "species",
  trait_cols = NULL,
  mode = c("columnwise", "rowwise"),
  numeric_method = c("bootstrap", "normal", "uniform"),
  keep_bounds = TRUE,
  inv_prefix = "inv",
  keep_species_column = TRUE,
  seed = NULL
)
```

## Arguments

- resident_traits:

  A data frame containing either a species identifier column (specified
  by `species_col`) or species identifiers as row names, plus one or
  more trait columns.

- n_inv:

  Integer; number of invaders to simulate (default `10`).

- species_col:

  `NULL` or character; name of the species identifier column in
  `resident_traits`. If `NULL`, species identifiers are taken from row
  names.

- trait_cols:

  `NULL` or character vector specifying which trait columns to use.
  Defaults to all columns except `species_col` when present.

- mode:

  Character; either `"columnwise"` or `"rowwise"`.

- numeric_method:

  Character; for numeric traits in columnwise mode, one of
  `"bootstrap"`, `"normal"`, or `"uniform"`.

- keep_bounds:

  Logical; if `TRUE`, normal or uniform draws are constrained to the
  observed minimum and maximum values.

- inv_prefix:

  Character; prefix used to construct invader identifiers.

- keep_species_column:

  Logical; if `TRUE` and `species_col` is not `NULL`, retain the species
  identifier column after setting row names.

- seed:

  `NULL` or integer; optional random seed for reproducibility.

## Value

A data frame of simulated invader traits. Row names correspond to
invader identifiers. If `species_col` is not `NULL` and
`keep_species_column = TRUE`, the species identifier column contains the
same identifiers.

## Details

**Species identifiers** Species identifiers can be supplied via a
dedicated column specified by `species_col`, or taken from the row names
of `resident_traits` when `species_col = NULL`. Newly simulated invaders
receive fresh, unique identifiers constructed from `inv_prefix`, which
become the row names of the returned data frame. When `species_col` is
not `NULL`, the same identifiers are also stored in that column unless
`keep_species_column = FALSE`.

**Trait selection** If `trait_cols` is `NULL`, all columns except
`species_col` (when present) are treated as traits. Otherwise, only the
intersection of `trait_cols` and existing column names is used.

**Sampling modes**

- `mode = "columnwise"`: Each trait is generated independently. Numeric
  traits can be drawn by bootstrap sampling, from a normal distribution
  parameterised by the empirical mean and standard deviation, or from a
  uniform distribution on the observed range. Factor and character
  traits are sampled from observed values.

- `mode = "rowwise"`: Entire rows are resampled with replacement from
  the resident trait table, preserving joint structure across traits.

**Identifier collisions** If proposed invader identifiers collide with
existing resident identifiers, they are made unique using
[make.unique](https://rdrr.io/r/base/make.unique.html).

## See also

[sample](https://rdrr.io/r/base/sample.html),
[make.unique](https://rdrr.io/r/base/make.unique.html),
[rnorm](https://rdrr.io/r/stats/Normal.html),
[runif](https://rdrr.io/r/stats/Uniform.html),
[sd](https://rdrr.io/r/stats/sd.html)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
residents <- data.frame(
  species  = paste0("sp", 1:5),
  height   = c(10.2, 9.8, 11.1, 10.5, 9.9),
  SLA      = c(15.0, 15.2, 14.7, 15.5, 15.1),
  lifeform = factor(c("tree", "shrub", "shrub", "tree", "herb"))
)

inv <- simulate_invaders(
  residents,
  n_inv = 4,
  species_col = "species",
  mode = "columnwise",
  numeric_method = "bootstrap"
)
head(inv)
} # }
```
