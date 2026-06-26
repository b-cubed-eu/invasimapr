# Getting started with invasimapr

## Overview

invasimapr is an R package to assess species and trait invasiveness and
site invasibility by visualising trait dispersion and computing invasion
fitness.

It implements a network invasibility framework that decomposes
**invasion fitness** – the per-capita growth rate of a rare invader
entering a resident community at equilibrium – into three mechanistic
components:

- **Abiotic suitability** (\\r^{(z)}\_{is}\\): how well the environment
  at site \\s\\ suits invader \\i\\.
- **Niche crowding** (\\C^{(z)}\_{is}\\): how much functional overlap
  invader \\i\\ has with the resident community at site \\s\\.
- **Resident competition** (\\S^{(z)}\_{is}\\): how saturated the
  resident community at site \\s\\ already is.

When abiotic benefits outweigh competitive penalties, an invader attains
positive invasion fitness and a higher probability of establishment.

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("b-cubed-eu/invasimapr")
```

## The workflow

The core workflow consists of eight wrapper functions, each returning an
updated `invasimapr_fit` container. The example below uses the demo data
shipped with the package.

``` r
library(invasimapr)

# Demo data: one long table with sites, coordinates, species counts,
# environment (env*) and trait (trait_*) columns.
site_env_spp <- read.csv(
  system.file("extdata", "site_env_spp_simulated.csv.gz", package = "invasimapr")
)
long_df <- site_env_spp
names(long_df)[names(long_df) == "site_id"] <- "site"
chr <- vapply(long_df, is.character, logical(1))
long_df[chr] <- lapply(long_df[chr], as.factor)

# 1. Prepare inputs (assemble and align core matrices)
fit <- prepare_inputs(
  long_df      = long_df,
  site_col     = "site",
  env_prefix   = "^env",
  trait_prefix = "^trait"
)

# 2. Simulate hypothetical invader profiles from the resident trait pool
traits_inv <- simulate_invaders(
  resident_traits = fit$inputs$traits_res,
  n_inv = 10, mode = "columnwise"
)

# 3. Trait space, convex hull and niche crowding
fit <- prepare_trait_space(fit, traits_inv = traits_inv)

# 4. Model resident abundances (GLMM) and derive standardised predictors
fit <- model_residents(fit)

# 5. Learn trait-dependent sensitivities
fit <- learn_sensitivities(fit)

# 6. Project invaders into the resident model space
fit <- predict_invaders(fit, traits_inv = traits_inv)

# 7. Compute invasion fitness and establishment probability
fit <- predict_establishment(fit, option = "C", prob_method = "probit")

# 8. Summarise invasiveness and invasibility (tables, maps, rankings)
fit <- summarise_results(fit)

fit
```

## Where to next

- The [tutorial
  articles](https://b-cubed-eu.github.io/invasimapr/articles/) walk
  through each step in detail with real data, including trait-space
  construction, clustering and risk scenarios, and the invasion-fitness
  synthesis.
- The [function
  reference](https://b-cubed-eu.github.io/invasimapr/reference/index.html)
  documents every wrapper and core function.

## Citing invasimapr

``` r
citation("invasimapr")
```
