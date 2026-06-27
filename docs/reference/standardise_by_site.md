# Standardize a site-by-species matrix by site (row-wise z) with optional robust SD

For each site (row), subtract the row-mean and divide by a row-scale.
The scale can be classical SD or a robust max(MAD, SD) with a small
floor.

## Usage

``` r
standardise_by_site(X, robust = FALSE)
```

## Arguments

- X:

  Numeric matrix (sites × species).

- robust:

  Logical. If TRUE, uses max(MAD, SD, 1e-8) per row; else classical SD
  with a safe 1 fallback.

## Value

List with `z` (z-scored matrix), `mu` (row means), `sd` (row scales).

## Examples

``` r
if (FALSE) { # \dontrun{
std = standardise_by_site(C_js, robust = TRUE)
C_js_z = std$z
} # }
```
