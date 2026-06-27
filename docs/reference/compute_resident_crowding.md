# Resident crowding \\C\_{js}\\ via Hellinger composition and Gower–Gaussian trait kernel

Given a site × resident community matrix and a resident traits table,
this function:

1.  normalizes site composition with Hellinger transformation
    (\\W\_{s\cdot}\\),

2.  computes Gower distances among residents (mixed types supported),

3.  converts distances to a Gaussian similarity kernel \\K\\ with
    bandwidth \\\sigma\_\alpha\\,

4.  returns the resident crowding matrix \\C = W \times K^\top\\ (sites
    × residents),

5.  optionally row-z-scores \\C\\ for within-site contrasts,

6.  optionally draws a heatmap of \\C\\ (row-z) and a geographic map of
    site-level summaries.

## Usage

``` r
compute_resident_crowding(
  comm_res,
  traits_res,
  site_df = NULL,
  overlay_sf = NULL,
  method_comp = "hellinger",
  distance_metric = "gower",
  sigma_alpha = NULL,
  row_z = TRUE,
  quantile_probs = 0.95,
  do_heatmap = TRUE,
  show_names = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  heatmap_palette = viridisLite::viridis,
  do_map = TRUE
)
```

## Arguments

- comm_res:

  A numeric matrix or data.frame of **sites × resident species**
  (non-negative abundances or transformed counts). Row names must be
  site IDs.

- traits_res:

  A data.frame of **resident species × traits**. Row names must be
  resident IDs that match the column names of `comm_res`. Mixed types
  allowed.

- site_df:

  Optional data.frame with columns `site`, `x`, `y` for mapping. If
  provided and `do_map=TRUE`, a site-level map is produced.

- overlay_sf:

  Optional `sf` object to outline on the map (e.g., study boundary).

- method_comp:

  Character method for
  [`vegan::decostand`](https://vegandevs.github.io/vegan/reference/decostand.html).
  Default `"hellinger"`.

- distance_metric:

  Character distance name for
  [`cluster::daisy`](https://rdrr.io/pkg/cluster/man/daisy.html).
  Default `"gower"`. You usually want `"gower"` for mixed traits.

- sigma_alpha:

  Optional numeric bandwidth for the Gaussian kernel. If `NULL`
  (default), uses the median of positive pairwise Gower distances. If
  that is missing or zero, falls back to `0.5`.

- row_z:

  Logical; if `TRUE` (default) return a row-z-scored version `C_js_z`
  for within-site contrasts.

- quantile_probs:

  Numeric vector of probabilities to summarize site-level crowding
  quantiles (e.g., `c(0.95)`). Use `NULL` to skip quantiles. Default
  `0.95`.

- do_heatmap:

  Logical; if `TRUE` (default) draw a clustered heatmap of `C_js_z`
  using pheatmap if available.

- show_names:

  Logical; show row/column names on the heatmap. Default `TRUE`.

- cluster_rows:

  Logical; cluster sites on the heatmap. Default `TRUE`.

- cluster_cols:

  Logical; cluster residents on the heatmap. Default `TRUE`.

- heatmap_palette:

  Function returning colors (e.g.,
  [`viridisLite::viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html)).
  Default
  [`viridisLite::viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html).

- do_map:

  Logical; if `TRUE`, draw a ggplot map of mean crowding per site when
  `site_df` is provided. Default `TRUE`.

## Value

A named list with:

- `W_site` — Hellinger composition matrix (sites × residents).

- `D_res` — resident × resident Gower distance matrix.

- `sigma_alpha` — Gaussian bandwidth used.

- `K_res_res` — Gaussian similarity kernel with zero diagonal.

- `C_js` — resident crowding matrix (sites × residents).

- `C_js_z` — row-z-scored version of `C_js` (if `row_z=TRUE`).

- `site_summary` — data.frame with per-site summaries (mean and
  requested quantiles).

- `heatmap` — `pheatmap` object or `NULL` if not drawn.

- `map` — `ggplot` object or `NULL` if not drawn.

## Details

Compute resident crowding from community composition and trait
similarity

The Hellinger transformation places compositions in a Euclidean-friendly
space. Gower distance handles mixed trait types. Distances are converted
to similarities with \\K\_{ij}=\exp(-D\_{ij}^2/(2\sigma\_\alpha^2))\\,
and the diagonal is set to 0 to remove self-crowding. Crowding for
resident \\j\\ at site \\s\\ is the weighted sum of trait-similar
resident composition at that site: \\C\_{sj} = \sum\_{r} W\_{sr}
K\_{jr}\\.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fake data: 6 sites × 5 residents
set.seed(42)
comm_res = matrix(rpois(30, lambda = 5), nrow = 6,
                   dimnames = list(paste0("s", 1:6), paste0("sp", 1:5)))
traits_res = data.frame(
  size = rnorm(5, 10, 2),
  diet = factor(sample(c("A","B"), 5, TRUE)),
  row.names = paste0("sp", 1:5)
)
site_df = data.frame(site = paste0("s", 1:6),
                      x = rep(1:3, each=2),
                      y = rep(1:2, 3))

out = compute_resident_crowding(
  comm_res   = comm_res,
  traits_res = traits_res,
  site_df    = site_df,
  quantile_probs = 0.95,
  do_heatmap = TRUE,
  do_map     = TRUE
)

str(out, max.level = 1)
out$heatmap  # pheatmap object
out$map      # ggplot map of mean C per site
} # }
```
