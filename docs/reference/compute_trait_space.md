# Shared trait-space construction (Gower + PCoA), resident hull, centroid, and density plot

Builds a unified trait space for resident and invader species by: (1)
computing a Gower distance on the stacked trait table, (2) projecting to
2D via classical MDS/PCoA
([`stats::cmdscale`](https://rdrr.io/r/stats/cmdscale.html)), (3)
deriving the resident convex hull (realised niche region), (4)
estimating a kernel density over the 2D space, and (5) optionally
rendering a base R `filled.contour` map that overlays the hull, all
points (residents black, invaders red), and the cloud centroid (white
square with black outline). Optionally it also computes a functional
dendrogram (`hclust`) and a pretty dendrogram plot
([`factoextra::fviz_dend`](https://rdrr.io/pkg/factoextra/man/fviz_dend.html))
when available.

## Usage

``` r
compute_trait_space(
  traits_res,
  traits_inv,
  k = 2,
  kde_n = 100,
  pad_prop = 0.1,
  main_title = "Trait space density with convex hull and centroid",
  legend_line = paste("Hull = realised niche; white square = centroid;",
    "black dots = residents; red dots = invaders"),
  cex_main = 1,
  cex_sub = 0.72,
  cex_lab = 0.85,
  cex_axis = 0.75,
  highlight_level = 0.5,
  do_plot = TRUE,
  do_dend = FALSE
)
```

## Arguments

- traits_res:

  A data.frame of resident species (rows) by traits (columns). Mixed
  types are supported (numeric, integer, factor, ordered) via Gower
  distance. Row names are used as resident IDs; if absent, sequential
  IDs are created.

- traits_inv:

  A data.frame of invader species (rows) by traits (columns). Must share
  the same trait columns (names and types) as `traits_res`. Row names
  are used as invader IDs; if absent, sequential IDs are created. If you
  have no invaders, pass a 0-row data.frame with matching columns.

- k:

  Integer embedding dimension for PCoA; default 2.

- kde_n:

  Integer grid size for 2D kernel density estimation (per axis). Default
  `100`.

- pad_prop:

  Numeric padding proportion added to the plotting range on each axis.
  Default `0.10` (10%).

- main_title:

  Character main title for the contour plot.

- legend_line:

  Character subtitle used as an inline legend in the title area.

- cex_main, cex_sub, cex_lab, cex_axis:

  Numeric text scaling for main title, subtitle, axis labels, and axis
  tick labels, respectively. Defaults `1, 0.72, 0.85, 0.75`.

- highlight_level:

  Numeric in (0,1\]; draws a bold contour at this proportion of the
  maximum density (e.g., `0.5`=50% of max). Set `NA` to skip.

- do_plot:

  Logical; if `TRUE` renders the `filled.contour` figure. Default
  `TRUE`.

- do_dend:

  Logical; if `TRUE`, computes `hclust` on Gower and attempts to produce
  a dendrogram plot with
  [`factoextra::fviz_dend`](https://rdrr.io/pkg/factoextra/man/fviz_dend.html)
  when available. Default `FALSE`.

## Value

A list with:

- `gower` - Gower distance matrix (all species).

- `scores` - Data frame of 2D PCoA scores (`Q_all`) with columns
  `tr1,tr2`.

- `Q_res`, `Q_inv` - Subsets of `scores` for residents and invaders.

- `hull_res` - Data frame (tr1,tr2) of the resident convex hull ring
  (closed), or `NULL` if \< 3 residents.

- `centroid` - Numeric vector of the overall (res+inv) centroid
  `c(tr1, tr2)`.

- `density` - List `list(x, y, z)` from
  [`MASS::kde2d`](https://rdrr.io/pkg/MASS/man/kde2d.html).

- `colors` - Named vector with color per species ID (residents black,
  invaders red).

- `hc` - `hclust` object (if `do_dend=TRUE`), else `NULL`.

- `dend_plot` - A ggplot object from
  [`factoextra::fviz_dend`](https://rdrr.io/pkg/factoextra/man/fviz_dend.html)
  when available, else `NULL`.

- `xlim`, `ylim` - Numeric length-2 ranges used for plotting with
  padding.

## Details

Compute and plot shared trait space for residents and invaders

This function assumes that the trait columns in `traits_inv` are
compatible with those in `traits_res`. Gower distance
([`cluster::daisy`](https://rdrr.io/pkg/cluster/man/daisy.html))
supports mixed data types and handles factors/ordered factors
appropriately. The PCoA (`cmdscale`) is computed on the Gower distances
and the first two axes are returned by default for visualisation. The
resident convex hull is computed on the resident scores only and is
drawn as the "realised niche region".

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal reproducible example with numeric traits
set.seed(1)
traits_res = data.frame(
  trait_1 = rnorm(20),
  trait_2 = rnorm(20, 2),
  row.names = paste0("sp", 1:20)
)
traits_inv = data.frame(
  trait_1 = rnorm(5, 0.5),
  trait_2 = rnorm(5, 2.5),
  row.names = paste0("inv", 1:5)
)

out = compute_trait_space(
  traits_res = traits_res,
  traits_inv = traits_inv,
  do_plot    = TRUE,
  do_dend    = TRUE,
  main_title = "Trait space density with hull and centroid",
  legend_line = "Hull = realised niche region;
  white square = centroid; black = residents; red = invaders"
)

# Access returned objects
head(out$scores)
out$centroid
if (!is.null(out$dend_plot)) print(out$dend_plot)
} # }
```
