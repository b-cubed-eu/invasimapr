# Prepare trait space and resident crowding (no site-z)

Orchestrates the *traits → space → centrality/hull → crowding* pipeline
for residents and invaders. Optionally standardises trait inputs,
constructs a joint trait space, computes convex hull / centroid /
centrality, and derives **raw** resident crowding \\C\_{js}\\ using a
Gaussian kernel (no per-site standardisation here). Row-wise z-scoring
is deferred to
[[`model_residents()`](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.md)](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.html).
Plot objects are *always* created and returned; they are only displayed
when `show_plots = TRUE`.

## Usage

``` r
prepare_trait_space(
  fit,
  traits_inv,
  crowding_sigma = NULL,
  show_plots = TRUE,
  do_standardise = TRUE,
  row_z = FALSE
)
```

## Arguments

- fit:

  An `invasimapr_fit` produced by
  [[`prepare_inputs()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_inputs.md)/[`assemble_matrices()`](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.html).
  Must contain `fit$inputs$traits_res`, `fit$inputs$comm_res`, and
  (optionally) `fit$inputs$site_df`.

- traits_inv:

  Data frame (or matrix) of **raw** invader traits (rows = invaders;
  columns aligned to resident trait columns).

- crowding_sigma:

  Optional numeric bandwidth for the resident crowding kernel. If
  `NULL`, a default/optimised value is chosen inside
  [[`compute_resident_crowding()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.html).

- show_plots:

  Logical. If `TRUE`, display plots as they are created. If `FALSE`,
  plots are still created/returned but not shown. Default `TRUE`.

- do_standardise:

  Logical. If `TRUE` and
  [[`standardise_model_inputs()`](https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.md)](https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.html)
  exists, standardise resident *and* invader trait inputs (environment
  handled later). Default `TRUE`.

- row_z:

  Logical. Whether to perform **row-wise** (site) z-scoring inside
  [`compute_resident_crowding()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.md).
  Leave `FALSE` here to keep **raw** `C_js`; row-z is typically applied
  in
  [[`model_residents()`](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.md)](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.html).

## Value

The input `invasimapr_fit` with updated components:

- `$traits`:

  List with `gower`, `Q_res`, `Q_inv`, `hull`, `centroid`, `density`,
  stored plots (`plots_ts`, `plots_ch`), centrality table, and hull
  vertices.

- `$crowding`:

  List with `W_site`, `D_res`, `sigma_alpha`, `K_res_res`, and **raw**
  `C_js` (unless `row_z = TRUE`).

- `$meta`:

  Lightweight cache of `residents`, `sites`, `invaders`, and
  `n_invaders`.

## Details

**Pipeline**

1.  *(Optional)* Standardise `traits_res/traits_inv` via
    [[`standardise_model_inputs()`](https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.md)](https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.html)
    when available; otherwise pass through raw traits.

2.  Build a **joint trait space** with
    [[`compute_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.html)
    (returns Gower distances, ordination scores for residents `Q_res`
    and invaders `Q_inv`, a density surface, dendrograms, etc.). Plots
    are requested but may be hidden depending on `show_plots`.

3.  Compute **centrality & hull** (resident convex hull, centroid,
    centrality scores) using
    [[`compute_centrality_hull()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_centrality_hull.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_centrality_hull.html).

4.  Compute **resident crowding** with
    [[`compute_resident_crowding()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.html),
    returning kernel weights `K_res_res`, distances `D_res`, site
    kernels `W_site`, chosen `sigma_alpha`, and **raw** `C_js` (unless
    `row_z = TRUE` is requested).

**Display control** To ensure reproducible plot objects without
cluttering the console, the function temporarily opens a null graphics
device when `show_plots = FALSE`.

**Assumptions & safeguards**

- Requires resident community and trait tables in `fit$inputs`.

- If standardisation fails or is unavailable, the function proceeds with
  raw traits.

- Site-wise z-scoring is intentionally **not** applied here by default.

## See also

- Trait space:
  [[`compute_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.html)

- Centrality & hull:
  [[`compute_centrality_hull()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_centrality_hull.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_centrality_hull.html)

- Resident crowding:
  [[`compute_resident_crowding()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.md)](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.html)

- Standardisation:
  [[`standardise_model_inputs()`](https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.md)](https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.html)

- Downstream modelling:
  [[`model_residents()`](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.md)](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute all plots but do not display them
fit2 <- prepare_trait_space(fit, traits_inv, show_plots = FALSE)

# Display plots during construction
fit3 <- prepare_trait_space(fit, traits_inv, show_plots = TRUE)

# Use a fixed kernel bandwidth and request row-z inside the crowding step
fit4 <- prepare_trait_space(fit, traits_inv, crowding_sigma = 0.35, row_z = TRUE)
} # }
```
