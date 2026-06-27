# invasimapr_fit object

A lightweight S3 container used throughout **invasimapr** to store
assembled inputs, fitted models, sensitivities, invasion fitness, and
derived summaries.

## Details

Objects of class `invasimapr_fit` are created by
[prepare_inputs](https://b-cubed-eu.github.io/invasimapr/reference/prepare_inputs.md)
(via
[new_invasimapr_fit](https://b-cubed-eu.github.io/invasimapr/reference/utils_internal.md))
and are progressively enriched by downstream modelling and prediction
functions.

The object typically contains the following components:

- inputs:

  Core data structures returned by
  [assemble_matrices](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md).

- meta:

  Basic metadata such as number of sites, residents, and traits.

- residents, traits, sensitivities, fitness, prob, summary:

  Optional components added by downstream workflows.
