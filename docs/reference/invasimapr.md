# Run the full invasimapr pipeline

Convenience wrapper that chains the stages with sensible defaults.

## Usage

``` r
invasimapr(
  traits_inv,
  option = "D",
  prob_method = "probit",
  prob_args = list(sigma = 1),
  ...
)
```

## Arguments

- traits_inv:

  Invader traits (rows = invaders).

- option:

  invasion-fitness option; see
  [`predict_establishment()`](https://b-cubed-eu.github.io/invasimapr/reference/predict_establishment.md).

- prob_method:

  NULL or "probit"/"logit"/"hard".

- prob_args:

  list forwarded to
  [`compute_establishment_probability()`](https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.md).

- ...:

  arguments passed to
  [`assemble_matrices()`](https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.md)
  in the very first stage.

## Value

An object of class `invasimapr_fit` containing the results of the full
pipeline (inputs, trait space, resident model, sensitivities, invader
predictions, invasion fitness and summaries).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- invasimapr(traits_inv = my_invader_traits, long_df = my_long_table)
} # }
```
