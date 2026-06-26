# Build residents-by-sites long table and fit the residents-only GLMM

Creates a long data set with one row per (site, resident) pair, attaches
standardized environment and resident traits, then fits a GLMM for
residents. Also returns fixed-effect predictions on the link scale as a
site×resident matrix (r_js) and the response-scale mean (mu_js).

## Usage

``` r
prep_resident_glmm(
  comm_res,
  env_df_z,
  traits_res_glmm,
  fml,
  family = glmmTMB::tweedie(link = "log"),
  seed = 123,
  fit_model = TRUE
)
```

## Arguments

- comm_res:

  Numeric matrix (sites × residents) of observed counts or abundance.
  Row names are site IDs; column names are resident species IDs.

- env_df_z:

  Data frame of standardized site-level environment. Row names are site
  IDs and must match `rownames(comm_res)`.

- traits_res_glmm:

  Data frame of standardized resident traits used in the GLMM. Row names
  are resident IDs and must match `colnames(comm_res)`.

- fml:

  A model formula for
  [`glmmTMB::glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)
  (e.g., from
  [`build_model_formula()`](https://b-cubed-eu.github.io/invasimapr/reference/build_model_formula.md)).

- family:

  A `glmmTMB` family object. Default is `glmmTMB::tweedie(link="log")`.

- seed:

  Integer. Random seed.

- fit_model:

  logical; if FALSE, return dat_r/fml without fitting.

## Value

list(dat_r, fit_r, grid_res, r_js, mu_js, sites, res_ids)

## Examples

``` r
if (FALSE) { # \dontrun{
# Called inside model_residents(); see that wrapper for the typical path.
res <- prep_resident_glmm(comm_res, env_df_z, traits_res_glmm, fml)
} # }
```
