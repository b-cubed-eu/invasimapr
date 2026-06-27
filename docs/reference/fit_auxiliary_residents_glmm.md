# Auxiliary GLMM for trait-varying and site-varying slopes

`fit_auxiliary_residents_glmm()` assembles a long table of resident site
by species cells with standardized predictors \\r^{(z)}\_{js}\\,
\\C^{(z)}\_{js}\\, and \\S^{(z)}\_{js}\\, together with two-dimensional
trait scores (`tr1`, `tr2`). It then fits a Gaussian generalized linear
mixed model to estimate how slopes on abiotic suitability, niche
crowding, and site saturation vary across the trait plane.

Optional site-level random slopes can be included for `r_z` and `C_z`.
Random slopes for `S_z` are intentionally excluded because `S_z` is
site-only and has no within-site variation across residents.

This function is intentionally forgiving: non-matrix inputs are coerced
using [`as.matrix()`](https://rdrr.io/r/base/matrix.html), missing
dimension names are repaired when possible, and inputs are aligned to
`comm_res`.

## Usage

``` r
fit_auxiliary_residents_glmm(
  comm_res,
  r_js_z,
  C_js_z,
  S_js_z,
  Q_res,
  use_site_random_slopes = TRUE,
  family = gaussian(),
  control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1e+05, eval.max = 1e+05)),
  na_action = c("drop", "error"),
  verbose = TRUE
)
```

## Arguments

- comm_res:

  Numeric matrix of site by resident abundances. Row names are site
  identifiers and column names are resident identifiers.

- r_js_z:

  Numeric matrix of standardized abiotic suitability with the same
  dimensions and names as `comm_res`.

- C_js_z:

  Numeric matrix of standardized niche crowding with the same dimensions
  and names as `comm_res`.

- S_js_z:

  Numeric matrix of standardized site saturation broadcast across
  residents, with the same dimensions and names as `comm_res`.

- Q_res:

  Data frame of resident trait scores. Must contain numeric columns
  `tr1` and `tr2`, with row names matching `colnames(comm_res)`.

- use_site_random_slopes:

  Logical; if `TRUE`, include random slopes for `r_z` and `C_z` at the
  site level.

- family:

  GLM family for the conditional model. Default is
  [`gaussian()`](https://rdrr.io/r/stats/family.html) applied to
  `log1p(abundance)`.

- control:

  Control object passed to
  [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html).

- na_action:

  How to handle rows with missing values. Either `"drop"` (default) or
  `"error"`.

- verbose:

  Logical; if `TRUE`, print a short model summary line.

## Value

A list with components:

- `fit`: the fitted `glmmTMB` model

- `data`: the long-format data frame used for fitting

- `formula`: the model formula

- `args`: resolved arguments used in the fit

## Details

Fit an auxiliary residents-only GLMM on standardized predictors

The fixed-effects component of the model is \$\$\log(1 +
\mathrm{abundance}) \sim (r_z + C_z + S_z) \times (tr1 + tr2),\$\$ with
random intercepts for species and site always included.

When `use_site_random_slopes = TRUE`, random slopes \\(0 + r_z \|
site)\\ and \\(0 + C_z \| site)\\ are added. Random slopes for `S_z` are
omitted because `S_z` is constant within each site and therefore not
identifiable.

## Role in invasion fitness

This auxiliary model estimates trait-dependent slope systems that can be
mapped to invader-level parameters such as \\\theta_i\\, \\\alpha_i\\,
and \\\beta_i\\, and optionally to site-varying counterparts via random
slopes. These parameters enter directly into the linear invasion-fitness
decomposition \$\$\lambda\_{is} = \Gamma\_{is} r^{(z)}\_{is} -
\alpha\_{is} C^{(z)}\_{is} - \beta_i S^{(z)}\_{is},\$\$ linking trait
position to abiotic suitability, niche crowding, and resident
competition.

## Examples

``` r
if (FALSE) { # \dontrun{
# Toy dimensions
set.seed(1)
S = 8  # sites
J = 6  # residents
sites   = paste0("s", 1:S)
res_ids = paste0("sp", 1:J)

# Site x resident abundance
comm_res = matrix(rpois(S*J, lambda = 2), S, J,
                   dimnames = list(sites, res_ids))

# Standardized predictors (row-wise z by site for r and C; site-only for S)
r_raw = matrix(rnorm(S*J), S, J, dimnames = dimnames(comm_res))
C_raw = matrix(rnorm(S*J), S, J, dimnames = dimnames(comm_res))
S_raw = matrix(rep(scale(rnorm(S))[,1], each = J), S, J,
                dimnames = dimnames(comm_res))
# Simple per-site z for demo
r_js_z = t(scale(t(r_raw))); r_js_z[!is.finite(r_js_z)] = 0
C_js_z = t(scale(t(C_raw))); C_js_z[!is.finite(C_js_z)] = 0
S_js_z = S_raw

# Resident trait-plane scores (PCoA on Gower in real workflow)
Q_res = data.frame(tr1 = rnorm(J), tr2 = rnorm(J))
rownames(Q_res) = res_ids

aux = fit_auxiliary_residents_glmm(
  comm_res   = comm_res,
  r_js_z     = r_js_z,
  C_js_z     = C_js_z,
  S_js_z     = S_js_z,
  Q_res      = Q_res,
  use_site_random_slopes = TRUE
)

summary(aux$fit)
head(aux$data)
aux$formula
} # }
```
