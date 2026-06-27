# Probabilistic establishment from invasion fitness

Invasion fitness \\\lambda\_{is}\\ integrates trait-space geometry
(distances, overlaps, convex hulls, centroids) with abiotic suitability
(alignment of invader traits to the local environment), niche crowding
(overlap with resident trait space weighted by composition), and
resident competition (site saturation).

`compute_establishment_probability()` maps \\\lambda\_{is}\\ to
probabilities of establishment using a unified interface:

- **Probit**: \\P = \Phi(\lambda / \sigma)\\, where \\\sigma\\ is a
  scalar, the residual standard deviation from a fitted auxiliary model,
  or a cell-wise predictive standard deviation.

- **Logistic**: \\P = \mathrm{logit}^{-1}(\lambda / \tau)\\, where
  \\\tau\\ is a scale parameter.

- **Hard rule**: \\P = I(\lambda \> 0)\\, yielding a binary map.

If `lambda_is` is not supplied, the function builds it from standardized
predictors using \\\lambda\_{is} = \gamma r^{(z)}\_{is} - \alpha
C^{(z)}\_{is} - \beta S^{(z)}\_{is} + \kappa\\.

## Usage

``` r
compute_establishment_probability(
  lambda_is = NULL,
  r_is_z = NULL,
  C_is_z = NULL,
  S_is_z = NULL,
  gamma = 1,
  alpha = NULL,
  beta = NULL,
  kappa = 0,
  method = c("probit", "logit", "hard"),
  sigma = NULL,
  tau = 1,
  fit = NULL,
  predictive = FALSE,
  sigma_mat = NULL,
  use_vcov = FALSE,
  Q_inv = NULL,
  site_df = NULL,
  return_long = TRUE,
  make_plots = TRUE,
  option_label = NULL
)
```

## Arguments

- lambda_is:

  Optional matrix of invasion fitness values with dimensions \\S\\ by
  \\I\\ (rows are sites, columns are invaders). If `NULL`, fitness is
  computed from the supplied components.

- r_is_z, C_is_z, S_is_z:

  Optional matrices of standardized abiotic suitability, niche crowding,
  and saturation. Used only when `lambda_is = NULL`.

- gamma:

  Optional vector of length \\I\\ or matrix of dimension \\S\\ by \\I\\
  giving the abiotic slope.

- alpha:

  Optional vector or matrix of crowding penalties. By convention, these
  values are constrained to be non-negative.

- beta:

  Optional vector of saturation penalties. Signed values allow
  facilitation.

- kappa:

  Optional scalar offset added to invasion fitness (default 0).

- method:

  Character string; one of `"probit"`, `"logit"`, or `"hard"`.

- sigma:

  Numeric scalar standard deviation for the probit transform.

- tau:

  Numeric scalar scale parameter for the logistic transform.

- fit:

  Optional fitted model object used to obtain a residual standard
  deviation when `method = "probit"`.

- predictive:

  Logical; if `TRUE`, a predictive standard deviation matrix is used for
  the probit transform.

- sigma_mat:

  Optional matrix of predictive standard deviations with the same
  dimensions as `lambda_is`.

- use_vcov:

  Logical; if `TRUE`, compute predictive standard deviations from the
  variance-covariance matrix of `fit`.

- Q_inv:

  Optional data frame of invader trait scores with columns `tr1` and
  `tr2`.

- site_df:

  Optional site metadata table with columns `site`, `x`, and `y`.

- return_long:

  Logical; if `TRUE`, include a long-format table in the output.

- make_plots:

  Logical; if `TRUE`, return diagnostic plots.

- option_label:

  Optional label attached to the output.

## Value

A list with components:

- `p_is`: matrix of establishment probabilities

- `lambda_is`: invasion fitness matrix

- `sigma_used`: standard deviation used by the probit transform

- `method`: transformation method

- `option_label`: label used for summaries

- `prob_long`: long-format table (optional)

- `plots`: list of plots (optional)

## Details

Convert invasion fitness to probabilistic establishment (probit, logit,
hard)

For the probit method, probabilities are computed as \\P\_{is} =
\Phi(\lambda\_{is} / \sigma)\\. A scalar or cell-wise standard deviation
may be used.

For the logistic method, probabilities are computed as \\P\_{is} =
\mathrm{logit}^{-1}(\lambda\_{is} / \tau)\\.

The hard rule returns a binary indicator equal to one when
\\\lambda\_{is} \> 0\\.

## Examples

``` r
## Minimal example (toy shapes)
set.seed(1)
S = 6; I = 4
sites   = paste0("s", 1:S)
inv     = paste0("i", 1:I)
r_is_z  = matrix(rnorm(S*I), S, I, dimnames=list(sites, inv))
C_is_z  = matrix(rnorm(S*I), S, I, dimnames=dimnames(r_is_z))
S_is_z  = matrix(rep(scale(rnorm(S)), each=I), S, I, dimnames=dimnames(r_is_z))
gamma   = setNames(runif(I, 0.5, 1.2), inv)
alpha   = setNames(runif(I, 0.2, 1.0), inv)
beta    = setNames(runif(I, 0.1, 0.6), inv)

# Build lambda internally, then get logistic probabilities
out_logit = compute_establishment_probability(
  r_is_z=r_is_z, C_is_z=C_is_z, S_is_z=S_is_z,
  gamma=gamma, alpha=alpha, beta=beta,
  method="logit", tau=1, return_long=FALSE, make_plots=FALSE
)
str(out_logit$p_is)
#>  num [1:6, 1:4] 0.274 0.57 0.302 0.93 0.509 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:6] "s1" "s2" "s3" "s4" ...
#>   ..$ : chr [1:4] "i1" "i2" "i3" "i4"

# Probit with a scalar sigma
out_probit = compute_establishment_probability(
  r_is_z=r_is_z, C_is_z=C_is_z, S_is_z=S_is_z,
  gamma=gamma, alpha=alpha, beta=beta,
  method="probit", sigma=1
)
# View site-mean probability map (requires ggplot2)
if (requireNamespace("ggplot2", quietly=TRUE)) print(out_probit$plots$site_mean)
#> NULL

# Hard rule (lambda>0)
out_hard = compute_establishment_probability(
  r_is_z=r_is_z, C_is_z=C_is_z, S_is_z=S_is_z,
  gamma=gamma, alpha=alpha, beta=beta,
  method="hard", return_long=TRUE, make_plots=FALSE
)
table(out_hard$prob_long$val)  # 0/1 counts
#> 
#>  0  1 
#> 10 14 
```
