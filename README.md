<!-- README.md is generated from README.Rmd. Please edit that file -->

# **`invasimapr`**<a href="https://b-cubed-eu.github.io/invasimapr/"><img src="man/figures/logo.png" align="right" height="139" alt="invasimapr website" /></a>

## A Novel Framework to visualise trait dispersion and assess species invasiveness or site invasibility

<!-- badges: start -->

[![repo
status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Release](https://img.shields.io/github/release/b-cubed-eu/invasimapr.svg)](https://github.com/b-cubed-eu/invasimapr/releases)
[![invasimapr status
badge](https://b-cubed-eu.r-universe.dev/invasimapr/badges/version)](https://b-cubed-eu.r-universe.dev/invasimapr)
[![CRAN
status](https://www.r-pkg.org/badges/version/invasimapr)](https://CRAN.R-project.org/package=invasimapr)
[![R-CMD-check](https://github.com/b-cubed-eu/invasimapr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/b-cubed-eu/invasimapr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/b-cubed-eu/invasimapr/graph/badge.svg)](https://app.codecov.io/gh/b-cubed-eu/invasimapr)
[![DOI](https://img.shields.io/badge/DOI-awaiting_upload_to_zenodo-orange)](https://zenodo.org)
[![name status
badge](https://b-cubed-eu.r-universe.dev/badges/:name?color=6CDDB4)](https://b-cubed-eu.r-universe.dev/)
[![MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)

<!-- badges: end -->

------------------------------------------------------------------------

## Introduction

Biological invasions are a leading driver of biodiversity loss. Establishment success depends on a species’ functional traits, the suitability of local environments, and the competitive pressure from resident communities—so ad-hoc, single-component analyses are insufficient. **`invasimapr`** provides a transparent, trait- and site-specific framework that integrates these components into a single, reproducible workflow to estimate **invasion fitness** and derive decision-ready indicators of **species invasiveness** and **site invasibility**.

At its core, the package (i) models **intrinsic growth potential** from trait–environment responses, (ii) quantifies **competitive penalties** imposed by resident communities via trait overlap and environmental filtering, and (iii) combines these to compute a site- and species-resolved fitness surface that can be summarised and mapped. It relies on widely used statistical tools (e.g., GLMM/GAM) and standard distance measures, making it accessible and extensible for applied invasion ecology and conservation planning.

---

## Core concepts (what the framework estimates)

* **Invasion fitness (`λ`)** — Net potential for a species to increase when rare at a site: `λ = r − C`, where `r` is predicted intrinsic (abiotic) performance and `C` is the competitive penalty from residents.
* **Invasiveness (`Vᵢ`)** — Propensity of a species to establish across sites (spatial aggregation of `λ`).
* **Invasibility (`Vₛ`)** — Openness of a site to establishment by newcomers (aggregation of `λ` over candidate invaders).

These are built from three linked pillars:

1. **Trait space → competition:** species are embedded in a functional trait space; a kernel (e.g., Gaussian) converts pairwise trait distances to competition coefficients (higher similarity → stronger competition).
2. **Environmental filtering:** residents matter most where they are well matched to local conditions; an environmental kernel up- or down-weights their effect by site–resident mismatch.
3. **Resident context:** predicted/typical resident abundance further scales their suppressive effect.

Together these define an **interaction tensor** that aggregates to a site-level penalty `C` and, with `r`, yields `λ`.

---

## What the package does (high-level workflow)

* **Data preparation**

  * Harmonise traits, environments, and resident composition.
  * Optionally **simulate invaders** to test “what-if” scenarios.

* **Model trait–environment responses**

  * Fit a single trait–environment model to predict `r` (expected performance without competitors) for candidate invaders at each site.
  * Estimate resident optima and site–resident mismatch for environmental weighting.

* **Quantify competitive pressure**

  * Build trait space and compute pairwise similarity → competition kernel.
  * Combine trait overlap, environmental match, and resident context into interaction strengths; sum over residents to get `C`.

* **Compute and summarise outcomes**

  * **Invasion fitness `λ`** for every species × site.
  * Site-level **invasibility (`Vₛ`)** and species-level **invasiveness (`Vᵢ`)** for mapping, ranking, and prioritisation.

The pipeline is **modular** (each step inspectable/reusable) and designed for **reproducibility**.

---

## Typical outputs

* Matrices/data frames for `r`, `C`, and `λ` (species × sites).
* Site summaries (**`Vₛ`**) and species summaries (**`Vᵢ`**) for reporting and maps.
* Intermediate diagnostics (trait distances, kernels, resident optima/mismatch) to audit assumptions and perform sensitivity checks.

---

## When to use `invasimapr`

* Screening **candidate invaders** or pathways and ranking species by establishment potential.
* Identifying **vulnerable sites** and allocating surveillance/management effort.
* **Scenario analysis** under environmental change (e.g., altered climates, trait shifts).
* Creating consistent, repeatable **maps of invasion risk** across large landscapes.

---

## Data requirements (minimum viable inputs)

* **Traits** for residents (and invaders/simulated invaders).
* **Site environments** (e.g., climate, soils, habitat metrics).
* **Resident composition** (occurrence/abundance or a proxy).
* Consistent **species and site identifiers** for joins.
* Optional: curated trait tables and metadata for automated ingestion.

---

## Main functions (overview, not a tutorial)

**Data & simulation**

* `get_trait_data()` — Collect, clean, and standardise trait data; optionally augment with metadata.
* `simulate_invaders()` — Generate hypothetical invaders to probe scenarios.

**Trait–environment modelling**

* `compute_trait_space()` — Build trait space and competition coefficients from pairwise distances.
* `build_glmm_formula()` — Compose model formulae for trait–environment responses.
* `predict_invader_response()` — Estimate intrinsic growth potential `r` for species at sites.

**Competition & environment**

* `compute_environment_kernel()` — Weight resident effects by site–resident environmental mismatch.
* `compute_interaction_strength()` — Combine trait overlap, environmental match, and resident context into pairwise impacts; sum to get `C`.

**Outcomes & summaries**

* `compute_invasion_fitness()` — Compute `λ = r − C` (with optional scaled variants).
* Summaries: **`Vₛ`** (site invasibility) and **`Vᵢ`** (species invasiveness) for mapping and prioritisation.

> For full argument lists and return types see the package reference index.

---

## Design principles & assumptions (brief)

* **Single coherent model:** one trait–environment fit underpins both invader performance and resident context to keep assumptions aligned.
* **Distance/kernels are explicit:** choice of trait/environment distance and kernel bandwidths (e.g., `σ_t`, `σ_e`) is transparent and tunable.
* **Interpretation depends on response:** if `r` is predicted abundance/occurrence, `λ` is a **relative establishment proxy**, not a demographic rate—interpret accordingly.
* **Auditability:** intermediate objects are returned so you can inspect sensitivity to traits chosen, scaling, and kernel parameters.

---

## Interoperability

* Plays well with common R ecosystems for spatial data, modelling, and visualisation.
* Complements packages used for data access/prep and for downstream mapping and reporting (e.g., building site layers and trait tables prior to modelling).

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("b-cubed-eu/invasimapr")
```

---

## Citation

If you use **`invasimapr`**, please cite the package and associated methods. See `citation("invasimapr")` and the repository’s CITATION files.

------------------------------------------------------------------------
