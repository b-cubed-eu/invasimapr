# Summarise invasiveness and invasibility (tables, maps, rankings)

Thin wrapper around
[summarise_invasiveness_invasibility](https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.md)
that extracts inputs from a
[new_invasimapr_fit](https://b-cubed-eu.github.io/invasimapr/reference/utils_internal.md)
container, forwards optional arguments, and optionally overlays a
boundary layer on the site map. This function is intended as a reporting
step after computing invasion fitness and/or establishment
probabilities.

## Usage

``` r
summarise_results(
  fit,
  boundary_sf = NULL,
  boundary_params = list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3),
  traits_inv = NULL,
  ...
)
```

## Arguments

- fit:

  An object produced by the invasimapr workflow, containing `fitness`
  and/or `prob` components created by upstream steps.

- boundary_sf:

  Optional object of class `sf` providing a boundary layer to overlay on
  the site map.

- boundary_params:

  Named list of aesthetics passed to
  [`ggplot2::geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html).
  Defaults to
  `list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3)`.

- traits_inv:

  Optional override of the invader trait table used for invader-ranked
  summaries. If `NULL`, traits are taken from the container.

- ...:

  Additional arguments forwarded to
  [summarise_invasiveness_invasibility](https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.md).

## Value

The input `fit` object with an added or updated `fit$summary` list.
Invisibly returns `fit` to support piping.

## Details

**Inputs used from the container**

- `fit$fitness$lambda_is`: site-by-invader invasion fitness matrix
  (optional).

- `fit$prob$p_is`: site-by-invader establishment probability matrix
  (optional).

- `fit$inputs$site_df`: site metadata with coordinates `x` and `y`
  (optional).

- `fit$inputs$comm_res`: resident community matrix for context summaries
  (optional).

- invader trait tables stored in the container, used for invader
  rankings.

At least one of `lambda_is` or `p_is` must be present; otherwise an
error is thrown.

**Behaviour**

1.  Extract invasion fitness and/or establishment probability matrices
    from the container.

2.  Select site metadata and an effective invader trait table when
    available.

3.  Call
    [summarise_invasiveness_invasibility](https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.md)
    to construct tidy summary tables and plots (site maps, rankings,
    heatmaps).

4.  If a boundary layer is supplied and plots are available, overlay it
    on the site map using `geom_sf()`.

**Output layout** The returned container gains a `summary` element
mirroring the structure returned by
[summarise_invasiveness_invasibility](https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.md),
typically containing summary tables, plots, and the arguments used.

## See also

[compute_invasion_fitness](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.md),
[compute_establishment_probability](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md),
[summarise_invasiveness_invasibility](https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.md),
[assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- prepare_inputs(
  sites = site_df,
  residents = resident_df,
  invaders = invader_df,
  traits = trait_df
)

fit <- learn_sensitivities(fit)
fit <- predict_invaders(fit, traits_inv = invader_traits)

fit <- summarise_results(fit)
fit$summary$plots$site_map
} # }
```
