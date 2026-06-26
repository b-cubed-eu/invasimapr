# Changelog

## invasimapr 0.2.0

Alignment with the [B-Cubed software development
guide](https://docs.b-cubed.eu/guides/software-development/) and
improved FAIR metadata.

### Metadata and documentation

- Added a persistent Zenodo archive (concept DOI
  [10.5281/zenodo.20842472](https://doi.org/10.5281/zenodo.20842472))
  and a DOI badge to the README.
- Added machine-readable metadata: `CITATION.cff` (CFF 1.2.0 with ORCID,
  DOI and qualified references), `codemeta.json` (CodeMeta/JSON-LD) and
  `.zenodo.json` (b3 community + Horizon Europe grant 101059592).
- Added a Darwin Core-aligned data dictionary
  (`inst/extdata/data_dictionary.csv`).
- DESCRIPTION now uses Title Case, includes the maintainer ORCID and a
  copyright holder (Stellenbosch University), and declares
  `scales`/`stringr` as imports.
- Rewrote the README as a concise, illustrated page with a minimal
  reproducible example using the bundled demo data.
- Restructured the documentation: one Getting-started vignette plus five
  tutorial articles under `vignettes/articles/`; renamed the bundled
  `inputs_vignettes.rds` to `invasimapr_vignettes.rds`.
- Added `CONTRIBUTING.md` and `LICENSE.md`; replaced the deprecated
  [`citEntry()`](https://rdrr.io/r/utils/citEntry.html) in
  `inst/CITATION` with
  [`bibentry()`](https://rdrr.io/r/utils/bibentry.html).

### Functions

- [`compute_invasion_fitness()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.md)
  gained an opt-in `standardise_inputs` argument to re-standardise raw
  predictors onto a common z-scale before computing fitness (off by
  default).

## invasimapr 0.1.0

First development release of invasimapr.

### Workflow

- Eight high-level wrapper functions covering the full pipeline:
  [`prepare_inputs()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_inputs.md),
  [`prepare_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.md),
  [`model_residents()`](https://b-cubed-eu.github.io/invasimapr/reference/model_residents.md),
  [`learn_sensitivities()`](https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.md),
  [`predict_invaders()`](https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.md),
  [`predict_establishment()`](https://b-cubed-eu.github.io/invasimapr/reference/predict_establishment.md)
  and
  [`summarise_results()`](https://b-cubed-eu.github.io/invasimapr/reference/summarise_results.md).
- 18 core functions, including
  [`assemble_matrices()`](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md),
  [`simulate_invaders()`](https://b-cubed-eu.github.io/invasimapr/reference/simulate_invaders.md),
  [`standardise_model_inputs()`](https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.md),
  [`compute_trait_space()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.md),
  [`compute_centrality_hull()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_centrality_hull.md),
  [`compute_resident_crowding()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.md),
  [`compute_site_saturation()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_site_saturation.md),
  [`derive_sensitivities()`](https://b-cubed-eu.github.io/invasimapr/reference/derive_sensitivities.md),
  [`site_varying_alpha_beta_gamma()`](https://b-cubed-eu.github.io/invasimapr/reference/site_varying_alpha_beta_gamma.md),
  [`build_invader_predictors()`](https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.md),
  [`compute_invasion_fitness()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.md),
  [`compute_establishment_probability()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md)
  and
  [`summarise_invasiveness_invasibility()`](https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.md).

### Features

- Decomposes invasion fitness into abiotic suitability, niche crowding
  and resident competition, following Hui et al. (2016).
- Multiple invasion-fitness model options (A-E) in
  [`compute_invasion_fitness()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.md).
- Establishment-probability mapping (probit, logistic and
  hard-threshold) via
  [`compute_establishment_probability()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md).
- Species invasiveness, site invasibility and trait-effect summaries via
  [`summarise_invasiveness_invasibility()`](https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.md).
- Trait data retrieval from external databases via
  [`get_trait_data()`](https://b-cubed-eu.github.io/invasimapr/reference/get_trait_data.md).
