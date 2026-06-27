# Compute trait-space centrality (robust Mahalanobis) and hull status, with visuals

Given 2D trait coordinates for **residents** and **invaders** (typically
PCoA scores), this function:

- estimates a robust centre and covariance for residents (MCD via
  [`MASS::cov.rob`](https://rdrr.io/pkg/MASS/man/cov.rob.html), fallback
  to classical),

- computes Mahalanobis distances for residents and invaders,

- converts distances to **centrality** using the resident distance CDF
  (`centrality = 1 - F(d)`; 1 = core, 0 = peripheral),

- determines whether invaders sit **inside** the resident convex hull
  (familiar) or **outside** (novel),

- returns a tidy table and three ggplot figures.

## Usage

``` r
compute_centrality_hull(
  Q_res,
  Q_inv,
  ellipse_level = 0.5,
  point_size = 2.8,
  alpha = 0.95,
  stroke = 0.7,
  rank_by = c("centrality", "distance"),
  peripheral_first = TRUE,
  palette = "viridis"
)
```

## Arguments

- Q_res:

  Data frame with resident coordinates; must contain columns `tr1`,
  `tr2`. Row names are treated as resident IDs.

- Q_inv:

  Data frame with invader coordinates; must contain `tr1`, `tr2`. Row
  names are treated as invader IDs. Can be 0-row.

- ellipse_level:

  Numeric in (0,1); confidence level for the resident normal ellipse in
  the trait-plane plot. Default `0.50`.

- point_size:

  Numeric; point size in the trait-plane plot. Default `2.8`.

- alpha:

  Numeric in (0,1\]; point alpha. Default `0.95`.

- stroke:

  Numeric; outline stroke width for points. Default `0.7`.

- rank_by:

  Character; one of `"centrality"` or `"distance"`. Controls the invader
  ranking plot. Default `"centrality"`.

- peripheral_first:

  Logical; if `TRUE` and `rank_by = "centrality"`, sort by increasing
  centrality (peripheral first). If ranking by distance, `TRUE` sorts by
  decreasing distance (far first). Default `TRUE`.

- palette:

  Centrality fill palette; passed to `scale_fill_viridis_c`. Default
  `"viridis"`.

## Value

A list with:

- `df` — tidy table with columns: `id`, `grp` (`"resident"|"invader"`),
  `tr1`, `tr2`, `d_md` (Mahalanobis), `d_eu` (Euclidean), `centrality`
  (0–1), `in_hull` (logical for invaders; residents = TRUE).

- `center`, `cov` — robust centre and covariance used.

- `hull_df` — closed ring (tr1,tr2) of resident convex hull, or `NULL`
  if \<3 residents.

- `plots` — list with: `p_trait` (trait-plane scatter), `p_dist`
  (distance distributions), `p_rank` (invader ranking with hull flags;
  `NULL` if no invaders).

## Details

Centrality and convex-hull membership in trait space

Mahalanobis distances may fail if the covariance is nearly singular; a
tiny ridge is added internally if needed. Hull membership is computed
with sp when available, otherwise with sf if available; if neither is
installed and residents \< 3, hull membership is returned as `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume you already computed PCoA coordinates:
# Q_res, Q_inv each with columns tr1, tr2 and rownames as IDs.
out = compute_centrality_hull(Q_res, Q_inv)
head(out$df)
out$plots$p_trait
out$plots$p_dist
out$plots$p_rank
} # }
```
