# Flexible formula constructor for residents-only trait–environment models

`build_model_formula()` assembles a right-hand side (RHS) for a GLMM (or
LM/GLM) from environment terms, trait terms, and (optionally) all
pairwise environment-by-trait interactions. It also appends
random-effects structures, such as `(1 | site) + (1 | species)` and
optional zero-correlation random slopes like `(0 + r_z || site)`.

You can pass the terms directly as character vectors, or let the
function derive them from `env_df` and/or `trait_df` column names.

## Usage

``` r
build_model_formula(
  response = "abundance",
  env_terms = NULL,
  trait_terms = NULL,
  env_df = NULL,
  trait_df = NULL,
  include_intercept = TRUE,
  include_env_main = TRUE,
  include_trait_main = TRUE,
  include_env_trait_interactions = TRUE,
  extra_fixed = NULL,
  random_intercepts = c("site", "species"),
  random_slopes = NULL,
  backend = c("glmmTMB", "lme4"),
  verbose = FALSE
)
```

## Arguments

- response:

  Character scalar. Name of the response on the LHS (default
  `"abundance"`).

- env_terms:

  Optional character vector of environment term names to include as
  fixed effects. If `NULL`, they can be derived from `env_df`.

- trait_terms:

  Optional character vector of trait term names to include as fixed
  effects. If `NULL`, they can be derived from `trait_df`.

- env_df:

  Optional data frame containing environment predictors; used only to
  infer `env_terms` when `env_terms` is `NULL`.

- trait_df:

  Optional data frame containing trait predictors; used only to infer
  `trait_terms` when `trait_terms` is `NULL`.

- include_intercept:

  Logical. Include the fixed-effect intercept? If `FALSE`, the intercept
  is removed via `0` (equivalent to `-1`). Default `TRUE`.

- include_env_main:

  Logical. Include environment main effects? Default `TRUE`.

- include_trait_main:

  Logical. Include trait main effects? Default `TRUE`.

- include_env_trait_interactions:

  Logical. Include all pairwise environment-by-trait interactions?
  Implemented as `(E1 + E2 + ...):(T1 + T2 + ...)`. Default `TRUE`.

- extra_fixed:

  Optional character vector of additional fixed-effect terms to append
  verbatim (e.g., `"poly(temp,2)"`, `"I(pH^2)"`).

- random_intercepts:

  Character vector of grouping factors for random intercepts (e.g.,
  `c("site","species")`). Use `NULL` to omit. Default
  `c("site","species")`.

- random_slopes:

  Named list of the form `list(site = c("r_z","C_z"), species = "r_z")`
  to add zero-correlation slopes `(0 + term || group)`. Use `NULL`
  (default) for none.

- backend:

  Character flag used only for messaging; both **lme4** and **glmmTMB**
  accept the same syntax here. Default `c("glmmTMB","lme4")`.

- verbose:

  Logical. If `TRUE`, prints the assembled formula string.

## Value

An object of class `formula`, e.g.:
`abundance ~ env1 + env2 + tr1 + tr2 + (env1 + env2):(tr1 + tr2) + (1 | site) + (1 | species) + (0 + r_z || site)`

## Details

Build a GLMM-ready model formula from trait and environment terms

## Examples

``` r
if (FALSE) { # \dontrun{
# Toy data
set.seed(1)
env_df_z        = data.frame(env1 = rnorm(10), env2 = rnorm(10))
traits_res_glmm = data.frame(tr1  = rnorm(5),  tr2  = rnorm(5))

fml = build_model_formula(
  response   = "abundance",
  env_df     = env_df_z,
  trait_df   = traits_res_glmm,
  random_intercepts = c("site","species"),
  random_slopes     = list(site = c("r_z","C_z"))
)
fml

# Fit a GLMM (example; requires your long residents×sites table `dat_r`)
# library(glmmTMB)
# fit = glmmTMB::glmmTMB(fml, family = glmmTMB::tweedie(link="log"), data = dat_r)
# summary(fit)
} # }
```
