# Learn sensitivities (alpha_i, beta_i, theta_i/gamma_i) and optional site-varying alpha_is, Gamma_is

Fits an auxiliary GLMM on resident data to estimate invader-level
sensitivities to crowding and saturation (\\\alpha_i\\, \\\beta_i\\),
and abiotic conversion slopes (\\\theta_i\\ or \\\gamma_i\\). When
supported by the data, optional site-level random slopes yield
site-varying adjustments (\\\alpha\_{is}\\, \\\Gamma\_{is}\\).

Results are written into the `fit$sensitivities` slot of an
[new_invasimapr_fit](https://b-cubed-eu.github.io/invasimapr/reference/utils_internal.md)
object for downstream invasion-fitness and establishment calculations.

## Usage

``` r
learn_sensitivities(fit, use_site_random_slopes = TRUE, lrt = TRUE)
```

## Arguments

- fit:

  An object produced by
  [prepare_inputs](https://b-cubed-eu.github.io/invasimapr/reference/prepare_inputs.md)
  and
  [assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md),
  containing resident predictor matrices (`r_js_z`, `C_js_z`, `S_js_z`),
  trait-space structures (`Q_res`, `Q_inv`), and the resident community
  layout (`inputs$comm_res`).

- use_site_random_slopes:

  Logical; if `TRUE`, the auxiliary model includes site-level random
  slopes for abiotic and crowding terms, enabling estimation of
  site-varying \\\alpha\_{is}\\ and \\\Gamma\_{is}\\ when supported by
  the data. Defaults to `TRUE`.

- lrt:

  Logical; if `TRUE`, compute Wald or likelihood-ratio tests for key
  contrasts (e.g., trait-varying versus global slopes). Defaults to
  `TRUE`.

## Value

The input `fit` object (invisibly), with an updated `fit$sensitivities`
list.

## Details

Fits an auxiliary GLMM on **resident** data to estimate invader-level
sensitivities to crowding and saturation (alpha_i, beta_i), and abiotic
conversion slopes (theta_i or gamma_i), with optional **site-varying**
random slopes that yield per-site adjustments (alpha_is, Gamma_is).
Results are written into the `fit$sensitivities` slot of an
[`invasimapr_fit`](https://b-cubed-eu.github.io/invasimapr/reference/invasimapr_fit.md)
object for downstream invasion-fitness and establishment steps.

**Workflow**

1.  Fit an auxiliary GLMM on resident responses using
    [fit_auxiliary_residents_glmm](https://b-cubed-eu.github.io/invasimapr/reference/fit_auxiliary_residents_glmm.md),
    optionally including site-level random slopes for \\r_z\\ and
    \\C_z\\.

2.  Convert GLMM coefficients to sensitivities (\\\alpha_i\\,
    \\\beta_i\\, \\\theta_i\\ or \\\gamma_i\\) using
    [derive_sensitivities](https://b-cubed-eu.github.io/invasimapr/reference/derive_sensitivities.md),
    returning signed and unsigned variants plus inference summaries.

3.  When supported, extract site-varying effects (\\\alpha\_{is}\\,
    \\\Gamma\_{is}\\) via
    [site_varying_alpha_beta_gamma](https://b-cubed-eu.github.io/invasimapr/reference/site_varying_alpha_beta_gamma.md).

The resulting components are stored in `fit$sensitivities`, including:

- global and trait-varying sensitivities;

- inference diagnostics and clamping summaries;

- optional site-varying matrices and compact decomposition tables.

## See also

[prepare_inputs](https://b-cubed-eu.github.io/invasimapr/reference/prepare_inputs.md),
[assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md),
[fit_auxiliary_residents_glmm](https://b-cubed-eu.github.io/invasimapr/reference/fit_auxiliary_residents_glmm.md),
[derive_sensitivities](https://b-cubed-eu.github.io/invasimapr/reference/derive_sensitivities.md),
[site_varying_alpha_beta_gamma](https://b-cubed-eu.github.io/invasimapr/reference/site_varying_alpha_beta_gamma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- prepare_inputs(
  sites = site_df,
  residents = resident_df,
  invaders = invader_df,
  traits = trait_df
)

fit <- learn_sensitivities(fit, use_site_random_slopes = TRUE)
names(fit$sensitivities)
} # }
```
