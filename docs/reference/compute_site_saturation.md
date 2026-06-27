# Compute site-only saturation \\S\_{s}\\ and global z-score \\S\_{s}^{(z)}\\

Implements three site-only saturation definitions and returns (i) the
raw site statistic \\S_s\\, (ii) its global mean/SD, (iii) the global
z-scored site statistic \\S_s^{(z)}\\, and (iv) a **broadcast**
site-resident matrix \\S\_{js}^{(z)}\\ for downstream modelling.

## Usage

``` r
compute_site_saturation(mode, comm_res, res_ids, mu_js = NULL)
```

## Arguments

- mode:

  One of
  `c("opportunity_penalty","modelled_dominance","evenness_deficit")`.

- comm_res:

  Site x resident matrix (or frame) of observed counts. Will be coerced
  to a numeric matrix; non-numeric entries become `NA`.

- res_ids:

  Character vector of resident IDs (column names to consider). Columns
  not found in `comm_res` are ignored with a warning.

- mu_js:

  Optional site x resident matrix (or frame) of expected abundances,
  required for `mode = "modelled_dominance"`. Will be coerced to numeric
  matrix.

## Value

A list with components:

- `S_s`:

  Numeric vector (sites) of saturation values.

- `S_mu`:

  Global mean of `S_s`.

- `S_sd`:

  Global SD of `S_s` (set to 1 if 0/NA).

- `S_s_z`:

  Global z-scores of `S_s`.

- `S_js_z`:

  Site x resident matrix broadcasting `S_s_z` across `res_ids`.

## Definitions

Let \\X\_{js}\\ be observed resident counts and \\\mu\_{js}\\ expected
counts.

- `"opportunity_penalty"`:

  Abundance penalty: \\S_s = \log(1 + \sum_j X\_{js})\\.

- `"modelled_dominance"`:

  Model-based crowding: \\S_s = \log(1 + \sum_j \mu\_{js})\\. Requires
  `mu_js`.

- `"evenness_deficit"`:

  Pielou's evenness deficit: \\J_s = H_s / \log S_s^{(\\)}\\ with
  Shannon \\H_s\\ and richness \\S_s^{(\\)}\\; then \\S_s = 1 - J_s\\.

After computing \\S_s\\, a global z-score is applied: \\S_s^{(z)} =
(S_s - \bar S) / \mathrm{sd}(S_s)\\, and this is broadcast across
residents to form \\S\_{js}^{(z)}\\ (same value for all \\j\\ within a
site \\s\\).

## See also

[Package reference:
compute_site_saturation()](https://b-cubed-eu.github.io/invasimapr/reference/compute_site_saturation.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Evenness deficit from observed community
sat = compute_site_saturation("evenness_deficit", comm_res, colnames(comm_res))
head(sat$S_s_z)

# Modelled dominance using expected abundances mu_js
sat2 = compute_site_saturation("modelled_dominance", comm_res, colnames(comm_res), mu_js = mu_js)
dim(sat2$S_js_z)
} # }
```
