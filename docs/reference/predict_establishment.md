# Compute invasion fitness (\\\lambda\\) and optional establishment probability (\\P\\)

Wraps
[[`compute_invasion_fitness()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.html)
to generate invader-by-site invasion fitness \\\lambda\_{is}\\ under
model Options **A-E**, and (optionally) maps \\\lambda\\ to
probabilities via
[[`compute_establishment_probability()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.html).
Supports calibration of \\\kappa\\ so that **resident** mean \\\lambda\\
is approximately zero, and can overlay an `sf` boundary on map-like
plots.

## Usage

``` r
predict_establishment(
  fit,
  option = c("A", "B", "C", "D", "E"),
  calibrate_kappa = FALSE,
  prob_method = c(NULL, "probit", "logit", "hard"),
  prob_args = list(),
  method = NULL,
  prob_scale = NULL,
  boundary_sf = NULL,
  boundary_params = list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3)
)
```

## Arguments

- fit:

  An `invasimapr_fit` from
  [[`predict_invaders()`](https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.md)](https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.html),
  containing `invaders$r_is_z/C_is_z/S_is_z`, resident summaries, and
  learned sensitivities.

- option:

  Character, one of `c("A","B","C","D","E")` selecting the
  invasion-fitness specification:

  A

  :   Baseline: \\\lambda = r^{(z)} - \alpha C^{(z)} - \beta S^{(z)}\\.

  B

  :   Global abiotic scaling \\\theta_0 \cdot r^{(z)}\\.

  C

  :   Trait-varying abiotic scaling \\\theta_i \cdot r^{(z)}\\.

  D

  :   Site-varying abiotic and crowding \\\Gamma\_{is}, \alpha\_{is}\\.

  E

  :   Signed saturation effect (facilitation allowed via signed
      \\\beta_i\\).

- calibrate_kappa:

  Logical; if `TRUE`, set \\\kappa\\ so mean **resident** \\\lambda
  \approx 0\\ (scale alignment for communication/comparison).

- prob_method:

  (legacy) `NULL` or one of `c("probit","logit","hard")`; preserved for
  backward compatibility.

- prob_args:

  (legacy) List of arguments passed to
  [[`compute_establishment_probability()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.html)
  (e.g., `sigma`, `tau`, `predictive`, `sigma_mat`, `site_df`,
  `return_long`, `make_plots`).

- method:

  Alias of `prob_method` (preferred in user code).

- prob_scale:

  Alias of `prob_args` (preferred in user code).

- boundary_sf:

  Optional **sf** object to overlay on map-like probability plots.

- boundary_params:

  Named list for
  [`ggplot2::geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html)
  aesthetics; default
  `list(inherit.aes=FALSE, fill=NA, color="black", size=0.3)`.

## Value

The input `fit` with:

- `$fitness`:

  List returned by
  [`compute_invasion_fitness()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.md)
  (including `lambda_is`, `lambda_long`, `option`, and any plots/maps).

- `$prob`:

  If requested, list returned by
  [`compute_establishment_probability()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md)
  with probability tables and plots (with boundary overlay applied when
  provided).

## Details

**What this wrapper does**

1.  Chooses the appropriate sensitivity inputs for the requested
    `option` (e.g., `theta_i` for C, `Gamma_is`/`alpha_is` for D).

2.  Calls ... with site-standardised inputs \\(r^{(z)}, C^{(z)},
    S^{(z)})\\ held in `fit`.

3.  (Optional) Converts \\\lambda\\ to probability \\P\\ via probit,
    logit, or a hard threshold, forwarding `prob_*` settings.

4.  (Optional) If `boundary_sf` is supplied and plots are available,
    overlays the boundary on map-like plots using `geom_sf()`.

**Calibration (`calibrate_kappa = TRUE`)** Aligns invader-resident
scales by shifting \\\lambda\\ so that the **resident** mean is ~0,
using resident moments and trait-plane slopes stored in `fit`.

## See also

- Fitness:
  [[`compute_invasion_fitness()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.html)

- Probabilities:
  [[`compute_establishment_probability()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.html)

- Sensitivities:
  [[`learn_sensitivities()`](https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.md)](https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.html)

- Invader predictors:
  [[`predict_invaders()`](https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.md)](https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit = fit |>
  predict_invaders(traits_inv) |>
  predict_establishment(option = "C", calibrate_kappa = TRUE,
                        method = "probit", prob_scale = list(sigma = 1))

# Overlay a boundary on probability maps
fit = predict_establishment(fit, option = "D", method = "logit",
                            prob_scale = list(tau = 1),
                            boundary_sf = rsa_boundary)
} # }
```
