# Build standardized invader predictors r_is_z, C_is_z, S_is_z (expects env/traits frames already projected to model columns)

Build standardized invader predictors r_is_z, C_is_z, S_is_z (expects
env/traits frames already projected to model columns)

## Usage

``` r
build_invader_predictors(
  fit_r,
  env_df_model,
  traits_inv_model,
  sites,
  inv_ids,
  r_mu_s,
  r_sd_s,
  W_site,
  gower_all,
  res_ids,
  sigma_alpha,
  C_mu_s,
  C_sd_s,
  S_s_z,
  verbose = TRUE
)
```

## Arguments

- fit_r:

  Fitted residents-only `glmmTMB` model.

- env_df_model:

  Data frame of environment predictors used by the model (rows = sites).
  Either `ENV_PC*` (if PCA was used) or standardized `env_df_z`.

- traits_inv_model:

  Data frame of invader traits used by the model (rows = invaders).
  Includes `TR_PC*` (if PCA was used) and any raw factor traits present
  in the model formula.

- sites:

  Character vector of site identifiers (must match row names of
  `env_df_model` and `W_site`).

- inv_ids:

  Character vector of invader identifiers (must match row names of
  `traits_inv_model`).

- r_mu_s:

  Numeric vector of site-wise means used to standardize abiotic
  suitability (`r_js`).

- r_sd_s:

  Numeric vector of site-wise standard deviations used to standardize
  abiotic suitability (`r_js`).

- W_site:

  Site × resident weighting matrix (e.g. relative abundance or
  presence–absence weights).

- gower_all:

  Square matrix or `dist` object of Gower distances among all residents
  and invaders.

- res_ids:

  Character vector of resident species identifiers.

- sigma_alpha:

  Numeric kernel bandwidth for crowding effects. If `NULL` or
  non-positive, a data-driven default is used.

- C_mu_s:

  Numeric vector of site-wise means used to standardize crowding
  (`C_is`).

- C_sd_s:

  Numeric vector of site-wise standard deviations used to standardize
  crowding (`C_is`).

- S_s_z:

  Numeric vector of standardized site saturation values.

- verbose:

  Logical; if `TRUE`, emit diagnostic messages.

## Value

A named list with components:

- r_is_z:

  Site-by-invader standardised abiotic-suitability matrix.

- C_is_z:

  Site-by-invader standardised niche-crowding matrix.

- S_is_z:

  Site-by-invader standardised resident-competition matrix.

- df:

  Long-format tidy data frame of the three standardised predictors.

- kernels:

  List with the invader-resident trait kernel `K_inv_res` and distance
  matrix `D_inv_res`.

- inputs_used:

  List recording `sites`, `inv_ids`, `res_ids` and the kernel bandwidth
  `sigma_alpha`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Called internally by predict_invaders(); see that function for usage.
inv <- build_invader_predictors(fit_r, env_df_model, traits_inv_model,
                                sites, inv_ids)
} # }
```
