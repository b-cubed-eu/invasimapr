# Build standardized invader predictors (r_is_z, C_is_z, S_is_z)

Uses resident-side moments and the residents-only model to construct
invader-level predictors on the **site-standardised** scale:
\\r^{(z)}\_{is}\\, \\C^{(z)}\_{is}\\, \\S^{(z)}\_{is}\\. Handles PCA
projection with “dummy” factor expansion so that invader covariates
exactly match the **training design** used in
[`model_residents()`](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.md)
/ the stored resident GLMM.

## Usage

``` r
predict_invaders(fit, traits_inv)
```

## Arguments

- fit:

  `invasimapr_fit` after
  [`prepare_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.md)
  /
  [`learn_sensitivities()`](https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.md)
  and resident model fitting (i.e., it contains resident GLMM artefacts,
  PCA objects/metadata, and resident moments such as `r_mu_s`, `r_sd_s`,
  `C_mu_s`, `C_sd_s`).

- traits_inv:

  **data.frame** of *raw* invader trait values; columns must match the
  resident trait names used at training (numerics and/or factors).

## Value

The input `fit` with an updated `fit$invaders` list containing:

- `traits_inv_raw`:

  Raw (unmodified) invader trait table (as supplied).

- `traits_inv_glmm`:

  Design-scale trait data passed to the resident model (post
  dummy-expansion, standardisation, PCA, and factor alignment).

- `r_is_z`, `C_is_z`, `S_is_z`:

  Site-standardised matrices for invaders.

- `df`:

  Tidy table used/returned by
  [`build_invader_predictors()`](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.md)
  for joins.

## Details

**Pipeline**

1.  **Environment**: if an environment PCA was used at training, rebuild
    the environment design matrix with the saved factor map and project
    using the stored rotation; otherwise reuse the stored `env_df_z`.

2.  **Traits**: coerce invader raw trait factors to training levels
    (redirect unseen levels to `_other_` if present); rebuild the dummy
    matrix and standardise with stored means/SDs; project via trait PCA
    if present; then append any raw factor terms that remained in the
    fixed-effects formula.

3.  **Predictors**: call
    [[`build_invader_predictors()`](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.md)](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.html)
    to compute \\r^{(z)}\_{is}\\, \\C^{(z)}\_{is}\\, \\S^{(z)}\_{is}\\
    using resident moments, crowding kernels, and similarity structures
    stored in `fit`.

**Assumptions & safeguards**

- Requires a trained resident model stored in `fit$residents$fit_r` and
  PCA metadata in `fit$model$*_pca*` if PCA was used.

- Uses resident moments (`r_mu_s`, `r_sd_s`, `C_mu_s`, `C_sd_s`) and
  similarity/crowding information (`W_site`, `gower`).

- Site ordering is taken from `fit$meta$sites`; inputs are conformed to
  it.

## Rationale

Many downstream steps (e.g. invasion fitness, establishment probability)
assume invader predictors live on the **same centred/scaled basis** as
the resident training data. This function guarantees alignment by: (i)
coercing raw trait factor levels to training levels (with `_other_`
fallbacks), (ii) rebuilding the design matrix with the training dummy
map, (iii) applying the stored centring/scaling, and (iv) projecting
through the stored PCA rotations before calling
[`build_invader_predictors()`](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.md).

## See also

- Input assembly:
  [[`assemble_matrices()`](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.html)

- Predictor builder:
  [[`build_invader_predictors()`](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.md)](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.html)

- Residents GLMM:
  [[`fit_auxiliary_residents_glmm()`](https://b-cubed-eu.github.io/invasimapr/reference/fit_auxiliary_residents_glmm.md)](https://b-cubed-eu.github.io/invasimapr/reference/fit_auxiliary_residents_glmm.html)

- Sensitivities:
  [[`learn_sensitivities()`](https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.md)](https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- prepare_inputs(sites = site_df, residents = resident_df,
                      invaders = invader_df, traits = trait_df)
fit <- learn_sensitivities(fit)

# New invader trait table (raw scale, same columns as residents' traits)
inv_traits <- data.frame(height = c(1.3, 0.9), SLA = c(12, 18),
                         life_form = c("shrub","forb"),
                         row.names = c("spA","spB"))

fit <- predict_invaders(fit, inv_traits)
str(fit$invaders, 1)
dim(fit$invaders$r_is_z)  # |sites| × |invaders|
} # }
```
