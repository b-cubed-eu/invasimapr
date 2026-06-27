# Site-varying alpha and gamma with trait-dependent beta

Constructs site-by-invader crowding penalties and density-dependence
scalars by combining trait-derived slopes with site-level random
effects.

Site-level random slopes for \\C_z\\ and \\r_z\\ are added to the
trait-dependent systems to produce site-varying crowding penalties
\\\alpha\_{is}\\ and density-dependence scalars \\\Gamma\_{is}\\.
Saturation effects \\\beta_i\\ are trait-varying only; no site-varying
\\\beta\_{is}\\ is constructed because \\S_z\\ is site-only.

The function also exposes signed slopes and clamping diagnostics
returned by
[`derive_sensitivities()`](https://b-cubed-eu.github.io/invasimapr/reference/derive_sensitivities.md).

## Usage

``` r
site_varying_alpha_beta_gamma(
  fit_coeffs,
  Q_inv,
  sites,
  inv_ids,
  lrt = TRUE,
  quiet = FALSE
)
```

## Arguments

- fit_coeffs:

  Fitted residents-only GLMM (e.g. `glmmTMB`) used to extract fixed
  effects and site-level random slopes.

- Q_inv:

  Data frame of invader trait scores with columns `tr1` and `tr2`; row
  names must correspond to `inv_ids`.

- sites:

  Character vector of site identifiers (rows of the output matrices).

- inv_ids:

  Character vector of invader identifiers (columns of the output
  matrices).

- lrt:

  Logical; passed to
  [`derive_sensitivities()`](https://b-cubed-eu.github.io/invasimapr/reference/derive_sensitivities.md)
  to control likelihood-ratio testing for trait effects.

- quiet:

  Logical; if `TRUE`, suppress warnings about missing or negligible
  site-level random effects.

## Value

A list with the following elements:

- `alpha_is`: matrix of site-by-invader crowding penalties

- `Gamma_is`: matrix of site-by-invader density-dependence scalars

- `beta_i`, `beta_signed_i`: invader-level saturation effects

- `alpha_signed_i`: signed crowding slopes prior to clamping

- `theta_i`, `gamma_i`: invader-level abiotic sensitivities

- `slope_C_i`: trait-derived crowding slopes

- `delta_C_s`, `delta_r_s`: site-level random effects

- `notes`: character vector of diagnostics

- `df`: tidy long-format data frame

## Examples

``` r
if (FALSE) { # \dontrun{
# Called inside learn_sensitivities() when site-varying effects are requested.
svg <- site_varying_alpha_beta_gamma(fit_coeffs, Q_inv, sites, inv_ids)
} # }
```
