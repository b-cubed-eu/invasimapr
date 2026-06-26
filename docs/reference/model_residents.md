# Resident GLMM and construction of standardized predictors

`model_residents()` fits a residents-only generalized linear mixed model
and constructs site-standardized predictors used for learning
sensitivities and projecting invaders. Specifically, the function:

1.  builds a fixed-effects design from environment and resident traits,

2.  optionally expands factors and applies PCA compression with stored
    centering and scaling metadata,

3.  fits a `glmmTMB` model to resident abundance,

4.  predicts resident linear predictors \\r\_{js}\\ and applies row-wise
    z-standardisation to obtain \\r^{(z)}\_{js}\\,

5.  standardises resident crowding \\C\_{js}\\ and computes site
    saturation summaries \\S_s\\, \\S^{(z)}\_s\\, and \\S^{(z)}\_{js}\\.

## Usage

``` r
model_residents(
  fit,
  family = glmmTMB::tweedie(link = "log"),
  include_env_trait_interactions = TRUE,
  saturation_mode = c("evenness_deficit", "opportunity_penalty", "modelled_dominance"),
  robust_r = TRUE,
  robust_c = TRUE,
  fit_model = TRUE,
  max_dense_gb = 8,
  reduce_strategy = c("auto", "none", "no_interactions", "pca"),
  pca_env_k = 5,
  pca_trait_k = 5,
  verbose = TRUE
)
```

## Arguments

- fit:

  An object of class `invasimapr_fit` produced by
  [`prepare_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.md).
  Must contain resident inputs, environment covariates, and raw crowding
  matrices.

- family:

  A `glmmTMB` family. Default is `glmmTMB::tweedie(link = "log")`.

- include_env_trait_interactions:

  Logical; whether to include environment by trait interactions when not
  reduced.

- saturation_mode:

  Character string controlling site saturation computation. Passed to
  [`compute_site_saturation()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_site_saturation.md).

- robust_r, robust_c:

  Logical; whether to use robust row-wise z-standardisation for
  \\r\_{js}\\ and \\C\_{js}\\.

- fit_model:

  Logical; if `FALSE`, perform preflight only and return without fitting
  the model.

- max_dense_gb:

  Numeric; maximum allowed dense fixed-effects size (GB) in automatic
  reduction mode.

- reduce_strategy:

  One of `"auto"`, `"none"`, `"no_interactions"`, or `"pca"`.

- pca_env_k, pca_trait_k:

  Integers giving the number of retained principal components for
  environment and trait blocks.

- verbose:

  Logical; if `TRUE`, emit progress messages.

## Value

The input `fit` object with updated components:

- `model`:

  Standardised inputs, PCA objects, and metadata required for invader
  projection.

- `residents`:

  Resident GLMM fit, prediction grid, raw and standardised predictor
  matrices, and site saturation summaries.

- `model$pca_used`:

  Flags indicating whether PCA was applied and the retained
  dimensionality.

## Details

Model residents and build standardized predictors

**Design construction**. Environment and trait tables are standardised
prior to model fitting. Character variables are coerced to factors to
ensure stable dummy expansion. Under PCA reduction, each block is
one-hot encoded, column-standardised, and passed to
[`prcomp()`](https://rdrr.io/r/stats/prcomp.html). All centering,
scaling, and rotation metadata are stored to enable reproducible
projection of invaders.

**Preflight memory checks**. In automatic reduction mode, the dense
fixed-effects design size is estimated before fitting. If the design
exceeds `max_dense_gb`, interactions are removed and/or PCA compression
is applied until the constraint is satisfied.

**Z-standardisation**. Row-wise z-standardisation stores site-level
means and standard deviations for later use in invasion fitness and
probability mapping.

## Reduction strategies

The argument `reduce_strategy` controls fixed-effects complexity.

- `"auto"`: estimate dense design size and iteratively reduce complexity
  until the memory budget is met.

- `"no_interactions"`: drop all environment by trait interactions.

- `"pca"`: dummy-expand, standardise, and apply
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) to environment and
  trait blocks, retaining the first components.

- `"none"`: retain the full requested fixed-effects structure.

## See also

- [`build_model_formula()`](https://b-cubed-eu.github.io/invasimapr/reference/build_model_formula.md)

- [`prep_resident_glmm()`](https://b-cubed-eu.github.io/invasimapr/reference/prep_resident_glmm.md)

- [`compute_site_saturation()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_site_saturation.md)

- [`prepare_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.md)

- [`predict_invaders()`](https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- prepare_inputs(sites = site_df, residents = resident_df,
                      invaders = invader_df, traits = trait_df)
fit <- prepare_trait_space(fit, traits_inv)
fit <- model_residents(fit, reduce_strategy = "auto", verbose = TRUE)

# Inspect z-standardised predictors for residents
dim(fit$residents$r_js_z); dim(fit$residents$C_js_z)
fit$residents$messages
} # }
```
