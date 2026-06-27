# Standardise model inputs (no leakage) for residents and invaders

Column-wise z-scores **environment** and **resident trait** numerics,
then scales **invader trait** numerics **using resident moments only**
(to avoid information leakage). Invader factor/character columns are
coerced to the resident levels; unseen levels become `NA`. Optionally
drops invader-only columns so the resident/invader trait schemas match.

## Usage

``` r
standardise_model_inputs(
  env_df = NULL,
  traits_res,
  traits_inv = NULL,
  drop_extra_invader_cols = FALSE,
  verbose = TRUE
)
```

## Arguments

- env_df:

  Optional `data.frame` (sites × environment). Numeric columns are
  z-scored; non-numeric are preserved.

- traits_res:

  `data.frame` (residents × traits). Mixed types allowed; numeric
  columns are z-scored.

- traits_inv:

  Optional `data.frame` (invaders × traits). Must contain at least the
  trait columns present in `traits_res`. Numeric columns are scaled
  using **resident** means/SDs; factors are coerced to resident levels.

- drop_extra_invader_cols:

  Logical; if `TRUE`, invader-only columns are dropped (not used
  downstream). If `FALSE`, they are still dropped for alignment but
  flagged in `info$notes`.

- verbose:

  Logical; print messages about what was standardised/coerced.

## Value

A named `list` with components:

- env_df_z:

  Environment table with numeric columns z-scored (or `NULL`).

- traits_res_glmm:

  Resident trait table with numeric columns z-scored.

- traits_inv_glmm:

  Invader trait table scaled to resident moments and factor levels
  matched (or `NULL`).

- moments:

  `list(env_means, env_sds, trait_means_res, trait_sds_res)` used for
  scaling.

- info:

  `list(notes=character())` with human-readable notes.

## Details

**What gets standardised and how**

- **Environment (`env_df`)**: numeric columns are z-scored (mean 0, sd
  1); non-numeric columns are kept as-is. Zero variance is guarded by
  setting sd=1.

- **Resident traits (`traits_res`)**: numeric columns are z-scored;
  mixed types allowed—non-numeric columns are kept.

- **Invader traits (`traits_inv`)**: numeric columns are scaled using
  the **resident** trait means/SDs only (never computed from invaders).
  Factor/ character columns are coerced to resident levels; unseen
  levels become `NA`. Extra invader columns are dropped (with a note).

**Returned objects**

- `env_df_z`: environment table with numeric columns standardised (or
  `NULL`)

- `traits_res_glmm`: resident traits with numeric columns standardised

- `traits_inv_glmm`: invader traits, scaled **like residents** + factor
  levels matched (or `NULL`)

- `moments`: resident/reference moments used for scaling (`env_*`,
  `trait_*`)

- `info$notes`: human-readable notes on coercions/drops

**Where this is used in the workflow**

- Called explicitly prior to GLMM fitting and when harmonising invaders,
  and implicitly by wrappers such as
  [`prepare_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.md)
  (if available).

## Invariants and guards

- Column names and row names are preserved.

- Zero-variance numeric columns use `sd=1` so z-scores stay defined.

- Invader trait numerics are always scaled by **resident** moments (no
  leakage).

- Invader extra columns are dropped for alignment; missing required
  columns error.

## See also

[`prepare_inputs()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_inputs.md),
[`prepare_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.md),
[`build_invader_predictors()`](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.md),
[`compute_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.md),
[`compute_resident_crowding()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.md)

## Examples

``` r
# Minimal reproducible example ----------------------------------------------
set.seed(1)
env_df = data.frame(elev = rnorm(5), temp = rnorm(5), zone = factor(sample(c("A","B"), 5, TRUE)))
rownames(env_df) = paste0("s", 1:5)

traits_res = data.frame(
  size = rlnorm(4), leaf = factor(c("broad","needle","broad","needle")),
  stringsAsFactors = FALSE
)
rownames(traits_res) = paste0("sp", 1:4)

traits_inv = data.frame(
  size = c(10, 1), leaf = factor(c("broad","unknown"))  # 'unknown' -> NA after coercion
)
rownames(traits_inv) = c("inv1","inv2")

std = standardise_model_inputs(env_df, traits_res, traits_inv, verbose = FALSE)
str(std, 1)
#> List of 5
#>  $ env_df_z       :'data.frame': 5 obs. of  3 variables:
#>  $ traits_res_glmm:'data.frame': 4 obs. of  2 variables:
#>  $ traits_inv_glmm:'data.frame': 2 obs. of  2 variables:
#>  $ moments        :List of 4
#>  $ info           :List of 1
head(std$traits_res_glmm)
#>           size   leaf
#> sp1 -0.3871555  broad
#> sp2 -0.4124477 needle
#> sp3 -0.6865489  broad
#> sp4  1.4861521 needle
head(std$traits_inv_glmm)   # note: 'leaf' for 'unknown' becomes NA
#>            size  leaf
#> inv1 31.8128564 broad
#> inv2  0.4876308  <NA>

# Workflow sketch ------------------------------------------------------------
# fit = prepare_inputs(long_df = longDF, ...)         # gives fit$inputs$env_df, $traits_res
# inv  = simulate_invaders(fit$inputs$traits_res, n_inv = 10)
# std  = standardise_model_inputs(fit$inputs$env_df, fit$inputs$traits_res, inv)
# std$traits_res_glmm; std$traits_inv_glmm              # pass to GLMM / trait-space steps
```
