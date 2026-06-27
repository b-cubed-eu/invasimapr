# Compute invasion fitness \\\lambda\_{is}\\ for multiple model options

`compute_invasion_fitness()` evaluates the invasion-fitness surface
\\\lambda\_{is}\\ over sites (\\s\\) and invaders (\\i\\) using three
standardized predictors:

- Abiotic suitability \\r^{(z)}\_{is}\\ (alignment between invader
  traits and the local environment),

- Niche crowding \\C^{(z)}\_{is}\\ (overlap with resident trait space,
  weighted by composition),

- Resident competition \\S^{(z)}\_{is}\\ (site-level saturation).

Five model variants are supported:

- Option A: \\\gamma = 1\\

- Option B: global slope \\\theta_0\\

- Option C: invader-specific slopes \\\theta_i\\

- Option D: site-varying slopes \\\Gamma\_{is}\\ and penalties
  \\\alpha\_{is}\\

- Option E: signed saturation effect

Optionally, an offset \\\kappa\\ can be calibrated so that the mean
resident invasion fitness is approximately zero when resident
standardized matrices are supplied.

## Usage

``` r
compute_invasion_fitness(
  r_is_z,
  C_is_z,
  S_is_z,
  option = c("A", "B", "C", "D", "E"),
  alpha_i = NULL,
  beta_i = NULL,
  theta0 = 1,
  theta_i = NULL,
  Gamma_is = NULL,
  alpha_is = NULL,
  beta_signed_i = NULL,
  calibrate_kappa = FALSE,
  r_js_z = NULL,
  C_js_z = NULL,
  S_js_z = NULL,
  Q_res = NULL,
  a0 = NULL,
  a1 = NULL,
  a2 = NULL,
  b0 = NULL,
  b1 = NULL,
  b2 = NULL,
  site_df = NULL,
  return_long = TRUE,
  label = NULL,
  standardise_inputs = FALSE
)
```

## Arguments

- r_is_z:

  Matrix of standardized abiotic suitability with dimensions \\S\\ by
  \\I\\.

- C_is_z:

  Matrix of standardized niche crowding with dimensions \\S\\ by \\I\\.

- S_is_z:

  Matrix of standardized site saturation with dimensions \\S\\ by \\I\\.

- option:

  Character string specifying the model option. One of `"A"`, `"B"`,
  `"C"`, `"D"`, or `"E"`.

- alpha_i:

  Optional named vector of invader-level crowding sensitivities.

- beta_i:

  Optional named vector of invader-level saturation sensitivities.

- theta0:

  Global abiotic slope used in options A, B, and E.

- theta_i:

  Optional invader-specific abiotic slopes used in option C.

- Gamma_is:

  Optional site by invader matrix of abiotic slopes used in option D.

- alpha_is:

  Optional site by invader matrix of crowding penalties used in option
  D.

- beta_signed_i:

  Optional signed saturation sensitivities used in option E.

- calibrate_kappa:

  Logical; if `TRUE`, compute a calibration offset using resident data.

- r_js_z, C_js_z, S_js_z:

  Optional resident standardized matrices required when
  `calibrate_kappa = TRUE`.

- Q_res:

  Optional data frame of resident trait-plane scores.

- a0, a1, a2, b0, b1, b2:

  Optional numeric coefficients used to derive resident analog slopes
  when calibrating \\\kappa\\.

- site_df:

  Optional site metadata table with columns `site`, `x`, and `y`.

- return_long:

  Logical; if `TRUE`, return a tidy long-format table.

- label:

  Optional character label attached to the output.

- standardise_inputs:

  Logical; if `TRUE`, re-standardise `r_is_z`, `C_is_z` and `S_is_z`
  onto a common global z-scale before computing fitness, so that no
  single term dominates. This is off by default because the `_z`
  arguments are expected to already be standardised; enable it only when
  passing raw (unstandardised) predictors. Errors if a predictor has
  zero or non-finite variance.

## Value

A list containing:

- `lambda_is`: invasion fitness matrix

- `GI`: abiotic slope used

- `AI`: crowding penalties used

- `BI`: saturation penalties used

- `kappa`: calibration offset

- `option`: option label

- `lambda_long`: tidy table (optional)

## Details

Compute invasion fitness with multiple model options

The invasion fitness is computed as:

- Option A: \\\lambda\_{is} = r^{(z)}\_{is} - \alpha_i C^{(z)}\_{is} -
  \beta_i S^{(z)}\_{is} + \kappa\\

- Option B: \\\lambda\_{is} = \theta_0 r^{(z)}\_{is} - \alpha_i
  C^{(z)}\_{is} - \beta_i S^{(z)}\_{is} + \kappa\\

- Option C: \\\lambda\_{is} = \theta_i r^{(z)}\_{is} - \alpha_i
  C^{(z)}\_{is} - \beta_i S^{(z)}\_{is} + \kappa\\

- Option D: \\\lambda\_{is} = \Gamma\_{is} r^{(z)}\_{is} - \alpha\_{is}
  C^{(z)}\_{is} - \beta_i S^{(z)}\_{is} + \kappa\\

- Option E: \\\lambda\_{is} = \theta_0 r^{(z)}\_{is} - \alpha_i
  C^{(z)}\_{is} + \beta_i^{(signed)} S^{(z)}\_{is} + \kappa\\

When `calibrate_kappa = TRUE`, the offset is chosen so that the mean
resident invasion fitness equals zero.

## Examples

``` r
## ---------------------------------------------------------------
## Minimal reproducible example (toy dimensions, no real ecology)
## ---------------------------------------------------------------

## S = number of sites, I = number of invaders
S = 5
I = 3
set.seed(1)

## r_is_z : intrinsic growth component
## rows = sites (s), columns = invaders (i)
r_is_z = matrix(
  rnorm(S * I),
  nrow = S,
  ncol = I,
  dimnames = list(paste0("s", 1:S), paste0("i", 1:I))
)

## C_is_z : crowding / competition component (same shape as r_is_z)
C_is_z = matrix(
  rnorm(S * I),
  nrow = S,
  ncol = I,
  dimnames = dimnames(r_is_z)
)

## S_is_z : site-only environmental gradient
## generated per-site, then broadcast across invaders
S_is_z = matrix(
  rep(scale(rnorm(S)), each = I),
  nrow = S,
  ncol = I,
  dimnames = dimnames(r_is_z)
)

## Invader-specific baseline parameters
alpha_i = setNames(runif(I, 0.2, 1.0), colnames(r_is_z))  # baseline crowding sensitivity
beta_i  = setNames(runif(I, 0.1, 0.5), colnames(r_is_z))  # strength of S effect

## ---------------------------------------------------------------
## Option A: fixed gamma = 1, no S effect (kappa = 0)
## ---------------------------------------------------------------
outA = compute_invasion_fitness(
  r_is_z, C_is_z, S_is_z,
  option      = "A",
  alpha_i     = alpha_i,
  beta_i      = beta_i,
  theta0      = 1,
  return_long = FALSE
)

## ---------------------------------------------------------------
## Option B: gamma shared across invaders (gamma = theta_0)
## ---------------------------------------------------------------
outB = compute_invasion_fitness(
  r_is_z, C_is_z, S_is_z,
  option      = "B",
  alpha_i     = alpha_i,
  beta_i      = beta_i,
  theta0      = 0.8,
  return_long = FALSE
)

## ---------------------------------------------------------------
## Option C: invader-specific gamma_i
## ---------------------------------------------------------------
theta_i = setNames(
  runif(I, 0.5, 1.2),
  colnames(r_is_z)
)

outC = compute_invasion_fitness(
  r_is_z, C_is_z, S_is_z,
  option      = "C",
  alpha_i     = alpha_i,
  beta_i      = beta_i,
  theta_i     = theta_i,
  return_long = FALSE
)

## ---------------------------------------------------------------
## Option D: fully site-varying Gamma_is and alpha_is
## ---------------------------------------------------------------

## gamma_is : site x invader matrix of density-dependence scalars
## constructed by repeating gamma_i across sites
Gamma_is = matrix(
  rep(theta_i, each = nrow(r_is_z)),
  nrow = nrow(r_is_z),
  ncol = ncol(r_is_z),
  dimnames = dimnames(r_is_z)
)

## alpha_is : site x invader crowding sensitivity
## start from invader-specific alpha_i and allow small site-level variation
alpha_is = matrix(
  NA_real_,
  nrow = nrow(r_is_z),
  ncol = ncol(r_is_z),
  dimnames = dimnames(r_is_z)
)

for (j in seq_len(ncol(r_is_z))) {
  alpha_is[, j] = alpha_i[j]
}

## Add site-level noise and enforce positivity.
## pmax() drops matrix attributes, so we re-wrap explicitly.
alpha_is = {
  tmp = pmax(
    0,
    alpha_is + matrix(
      rnorm(length(alpha_is), 0, 0.1),
      nrow = nrow(r_is_z),
      ncol = ncol(r_is_z)
    )
  )
  matrix(
    tmp,
    nrow = nrow(r_is_z),
    ncol = ncol(r_is_z),
    dimnames = dimnames(r_is_z)
  )
}

outD = compute_invasion_fitness(
  r_is_z, C_is_z, S_is_z,
  option      = "D",
  alpha_is    = alpha_is,
  beta_i      = beta_i,
  Gamma_is    = Gamma_is,
  return_long = FALSE
)

## ---------------------------------------------------------------
## Option E: signed S effect (positive or negative beta_i)
## ---------------------------------------------------------------
beta_signed_i = setNames(
  rnorm(I, 0, 0.3),
  colnames(r_is_z)
)

outE = compute_invasion_fitness(
  r_is_z, C_is_z, S_is_z,
  option          = "E",
  alpha_i         = alpha_i,
  beta_signed_i   = beta_signed_i,
  theta0          = 1,
  return_long     = FALSE
)
```
