# Summaries of invasion fitness: species invasiveness, trait effects, and site invasibility

Invasion fitness \\\lambda\_{is}\\ integrates structure in trait space
(distances, overlaps, convex hull, cloud centroid) with abiotic
suitability (alignment to environment), niche crowding (overlap with
residents weighted by composition), and resident competition (site
saturation). `summarise_invasiveness_invasibility()` collapses the
site-by-invader surface into species-, trait-, and site-level summaries
using either a probabilistic measure \\P(\lambda \> 0)\\ or a hard rule
\\I\\\lambda \> 0\\\\.

Species invasiveness (per invader) summarises the breadth of
establishment across sites: \$\$V_i = \|S\|^{-1} \sum_s I\\\lambda\_{is}
\> 0\\\$\$ or \$\$\tilde V_i = \|S\|^{-1} \sum_s P(\lambda \> 0 \mid i,
s).\$\$

Site invasibility (per site) quantifies openness to newcomers: \$\$V_s =
\|I\|^{-1} \sum_i I\\\lambda\_{is} \> 0\\\$\$ or \$\$\tilde V_s =
\|I\|^{-1} \sum_i P(\lambda \> 0 \mid i, s).\$\$

Trait invasiveness scores which invader traits explain variation in
invasiveness, using standardized slopes for continuous traits and ANOVA
\\R^2\\ for categorical traits.

## Usage

``` r
summarise_invasiveness_invasibility(
  lambda_is = NULL,
  p_is = NULL,
  use_probabilistic = FALSE,
  prob_threshold = 0.5,
  site_df = NULL,
  traits_inv = NULL,
  comm_res = NULL,
  return_long = TRUE,
  make_plots = TRUE,
  label = NULL
)
```

## Arguments

- lambda_is:

  Matrix of invasion fitness with dimensions sites by invaders. Row and
  column names are required.

- p_is:

  Optional matrix of establishment probabilities with the same
  dimensions as `lambda_is`.

- use_probabilistic:

  Logical. If `TRUE`, summaries use `p_is`. If `FALSE`, summaries use
  the hard rule \\I\\\lambda \> 0\\\\.

- prob_threshold:

  Numeric between 0 and 1. If probabilistic summaries are used and a
  binary view is requested, values with `p_is >= prob_threshold` are
  treated as 1.

- site_df:

  Optional data frame with columns `site`, `x`, and `y` for spatial
  mapping.

- traits_inv:

  Optional data frame of invader traits; row names must match invader
  IDs. Numeric and factor traits are supported.

- comm_res:

  Optional site-by-resident matrix used to compute resident richness per
  site.

- return_long:

  Logical. If `TRUE`, returns a tidy long table of the site-by-invader
  surface.

- make_plots:

  Logical. If `TRUE`, returns ggplot objects.

- label:

  Optional character label added to plot titles and metadata.

## Value

A list with the following components:

- `species`: data frame of invader-level invasiveness metrics

- `site`: data frame of site-level invasibility metrics

- `trait_effects`: data frame of trait effect sizes

- `establish_long`: tidy long-format table of the working surface

- `plots`: list of ggplot objects (or `NULL`)

- `meta`: list describing the summary mode and threshold

## Details

Summarize invasion fitness into species, trait, and site metrics (with
plots)

## Examples

``` r
set.seed(42)
S = 8; I = 5
sites = paste0("s", 1:S)
inv   = paste0("i", 1:I)

# Fake fitness and probabilities
lambda_is = matrix(rnorm(S*I, sd=1), S, I, dimnames=list(sites, inv))
p_is      = pnorm(lambda_is)  # crude probit just for the example

# Minimal site coordinates for plotting (optional)
site_df = data.frame(site = sites,
                      x = rep(1:4, each=2)[1:S],
                      y = rep(1:2, times=4)[1:S])

# Minimal trait table for invaders (one numeric, one factor)
traits_inv = data.frame(trait_size = runif(I, 0, 1),
                         trait_type = factor(sample(c("A","B"), I, TRUE)),
                         row.names = inv)

# 1) Probabilistic summaries (use p_is)
outP = summarise_invasiveness_invasibility(
  lambda_is = lambda_is,
  p_is      = p_is,
  use_probabilistic = TRUE,
  site_df   = site_df,
  traits_inv = traits_inv,
  make_plots = FALSE
)
names(outP)
#> [1] "species"        "site"           "trait_effects"  "establish_long"
#> [5] "plots"          "meta"          
head(outP$species)
#> # A tibble: 5 × 5
#>   invader   V_i n_sites V_i_hard n_hard
#>   <chr>   <dbl>   <int>    <dbl>  <dbl>
#> 1 i1      0.636       8    0.625      5
#> 2 i2      0.625       8    0.5        4
#> 3 i3      0.380       8    0.25       2
#> 4 i4      0.514       8    0.5        4
#> 5 i5      0.349       8    0.375      3

# 2) Hard rule summaries (use I{lambda>0})
outH = summarise_invasiveness_invasibility(
  lambda_is = lambda_is,
  use_probabilistic = FALSE,
  site_df = site_df,
  traits_inv = traits_inv,
  make_plots = FALSE
)
head(outH$site)
#> # A tibble: 6 × 9
#>   site    V_s n_inv total_expected     x     y V_s_hard n_est n_fail
#>   <chr> <dbl> <int>          <dbl> <int> <int>    <dbl> <dbl>  <dbl>
#> 1 s1      0.8     5              4     1     1      0.8     4      1
#> 2 s2      0       5              0     1     2      0       0      5
#> 3 s3      0.6     5              3     2     1      0.6     3      2
#> 4 s4      0.6     5              3     2     2      0.6     3      2
#> 5 s5      0.4     5              2     3     1      0.4     2      3
#> 6 s6      0       5              0     3     2      0       0      5
```
