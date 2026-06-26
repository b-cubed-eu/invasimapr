# Derive trait-varying sensitivities (αᵢ, βᵢ) and abiotic slope (θ, Γᵢ)

From an auxiliary GLMM on standardized predictors, extract fixed-effect
systems for `r_z`, `C_z`, and `S_z`, evaluate them at invader trait
coordinates, and (optionally) test whether the abiotic slope on `r_z`
should vary with traits (`tr1`, `tr2`). Biotic terms are enforced as
nonnegative penalties on invasion fitness.

## Usage

``` r
derive_sensitivities(fit_coeffs, Q_inv, inv_ids, lrt = TRUE)
```

## Arguments

- fit_coeffs:

  A fitted `glmmTMB` model on `(r_z + C_z + S_z) × (tr1 + tr2)`.

- Q_inv:

  Data frame with invader trait coordinates `tr1`, `tr2` (rownames =
  invader IDs).

- inv_ids:

  Character vector of invader IDs specifying the output order.

- lrt:

  Logical **or** character. If logical:

  - `TRUE` → Wald test; `FALSE` → no test. If character, one of
    `"wald"`, `"lrt"`, or `"none"`. Default is equivalent to `"wald"`
    for `TRUE` and `"none"` for `FALSE`.

## Value

A list with components:

- alpha_i, beta_i:

  Nonnegative penalties for `C_z` and `S_z`.

- alpha_signed_i, beta_signed_i:

  Signed (pre-clamp) slopes.

- theta0:

  Intercept for `r_z` slope system.

- theta_i:

  Trait-varying `r_z` slope evaluated at invader traits.

- gamma_i:

  Either `theta_i` (if test is significant) or a constant vector equal
  to `theta0` (if not).

- a0,a1,a2; b0,b1,b2:

  Coefficient components for `C_z` and `S_z` systems: main effect and
  interactions with `tr1`,`tr2`.

- df:

  Per-invader tidy data frame with traits, signed slopes, clamping
  flags, and `theta_i`/`gamma_i`.

- clamp_summary:

  Counts and proportions of clamped \\\alpha\\ and \\\beta\\.

- use_theta:

  Logical; `TRUE` if the chosen test was significant.

- lrt:

  **Back-compat alias** of `test_tab` (data frame with test result).

- test_tab:

  Data frame with test statistic, df, and p-value (or `NULL`).

- test_type:

  One of `"wald"`, `"lrt"`, or `"none"`.

- prior_note:

  Reminder of the nonnegative-penalty prior.

## Details

**Biotic penalties (prior/constraint).** Invasion fitness is modeled as
\$\$\lambda\_{is} = \Gamma\_{is} \\ r^{(z)}\_{is} \\-\\ \alpha_i \\
C^{(z)}\_{is} \\-\\ \beta_i \\ S^{(z)}\_{is}.\$\$ We therefore set
\\\alpha_i = \max(0, -\text{slope}\_{C,i})\\ and \\\beta_i = \max(0,
-\text{slope}\_{S,i})\\, where \\\text{slope}\_{C,i}\\ and
\\\text{slope}\_{S,i}\\ are the fitted (possibly trait-varying) slopes
for `C_z` and `S_z`. Signed versions (before clamping) are returned for
diagnostics.

**Testing trait-varying abiotic slope.** The argument `lrt` controls
whether and how to test the joint null \\H_0: r_z:tr1 = r_z:tr2 = 0\\.

- For backward compatibility, `lrt = TRUE` performs a **Wald** joint χ²
  test using the fixed-effects variance–covariance (no refit).

- `lrt = FALSE` performs no test.

- Alternatively, pass a character string: `"wald"`, `"lrt"`, or
  `"none"`.

  - `"lrt"` refits a reduced model (dropping both interactions) and runs
    a likelihood-ratio test (`anova(..., test = "Chisq")`).

If the chosen test is significant at \\\alpha = 0.05\\, the function
uses the trait-varying `theta_i`; otherwise a constant `theta0` is used
for `gamma_i`.

## See also

[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html),
[`stats::anova()`](https://rdrr.io/r/stats/anova.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Called inside learn_sensitivities(); see that wrapper for the typical path.
sens <- derive_sensitivities(fit_coeffs, Q_inv, inv_ids, lrt = TRUE)
} # }
```
