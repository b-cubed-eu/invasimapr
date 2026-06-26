# invasimapr 0.2.0

Alignment with the
[B-Cubed software development guide](https://docs.b-cubed.eu/guides/software-development/)
and improved FAIR metadata.

## Metadata and documentation

* Added a persistent Zenodo archive (concept DOI
  [10.5281/zenodo.20842472](https://doi.org/10.5281/zenodo.20842472)) and a DOI
  badge to the README.
* Added machine-readable metadata: `CITATION.cff` (CFF 1.2.0 with ORCID, DOI and
  qualified references), `codemeta.json` (CodeMeta/JSON-LD) and `.zenodo.json`
  (b3 community + Horizon Europe grant 101059592).
* Added a Darwin Core-aligned data dictionary
  (`inst/extdata/data_dictionary.csv`).
* DESCRIPTION now uses Title Case, includes the maintainer ORCID and a copyright
  holder (Stellenbosch University), and declares `scales`/`stringr` as imports.
* Rewrote the README as a concise, illustrated page with a minimal reproducible
  example using the bundled demo data.
* Restructured the documentation: one Getting-started vignette plus five tutorial
  articles under `vignettes/articles/`; renamed the bundled
  `inputs_vignettes.rds` to `invasimapr_vignettes.rds`.
* Added `CONTRIBUTING.md` and `LICENSE.md`; replaced the deprecated `citEntry()`
  in `inst/CITATION` with `bibentry()`.

## Functions

* `compute_invasion_fitness()` gained an opt-in `standardise_inputs` argument to
  re-standardise raw predictors onto a common z-scale before computing fitness
  (off by default).

# invasimapr 0.1.0

First development release of invasimapr.

## Workflow

* Eight high-level wrapper functions covering the full pipeline:
  `prepare_inputs()`, `prepare_trait_space()`, `model_residents()`,
  `learn_sensitivities()`, `predict_invaders()`, `predict_establishment()`
  and `summarise_results()`.
* 18 core functions, including `assemble_matrices()`, `simulate_invaders()`,
  `standardise_model_inputs()`, `compute_trait_space()`,
  `compute_centrality_hull()`, `compute_resident_crowding()`,
  `compute_site_saturation()`, `derive_sensitivities()`,
  `site_varying_alpha_beta_gamma()`, `build_invader_predictors()`,
  `compute_invasion_fitness()`, `compute_establishment_probability()` and
  `summarise_invasiveness_invasibility()`.

## Features

* Decomposes invasion fitness into abiotic suitability, niche crowding and
  resident competition, following Hui et al. (2016).
* Multiple invasion-fitness model options (A-E) in `compute_invasion_fitness()`.
* Establishment-probability mapping (probit, logistic and hard-threshold) via
  `compute_establishment_probability()`.
* Species invasiveness, site invasibility and trait-effect summaries via
  `summarise_invasiveness_invasibility()`.
* Trait data retrieval from external databases via `get_trait_data()`.
