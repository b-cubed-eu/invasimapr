#' Learn sensitivities (αᵢ, βᵢ, θᵢ/γᵢ) and optional site-varying αᵢₛ, Γᵢₛ
#'
#' Fits an auxiliary GLMM on **resident** data to estimate invader-level
#' sensitivities to crowding and saturation (αᵢ, βᵢ), and abiotic conversion
#' slopes (θᵢ or γᵢ), with optional **site-varying** random slopes that yield
#' per-site adjustments (αᵢₛ, Γᵢₛ). Results are written into the
#' `fit$sensitivities` slot of an [`invasimapr_fit`] object for downstream
#' invasion-fitness and establishment steps.
#'
#' @param fit An object of class `invasimapr_fit` produced by
#'   [`prepare_inputs()`] / [`assemble_matrices()`], containing resident
#'   matrices (`r_js_z`, `C_js_z`, `S_js_z`), trait-space structures (`Q_res`,
#'   `Q_inv`), and resident community layout (`inputs$comm_res`).
#' @param use_site_random_slopes Logical; if `TRUE`, the auxiliary model is fit
#'   with site-level random slopes for the abiotic and crowding terms, enabling
#'   estimation of site-varying αᵢₛ and Γᵢₛ when supported by the data.
#'   Defaults to `TRUE`.
#' @param lrt Logical; if `TRUE`, compute Wald/LRT summaries for key contrasts
#'   (e.g., trait-varying vs global slopes) to guide model choice and reporting.
#'   Defaults to `TRUE`.
#'
#' @details
#' **Workflow**
#' 1. Build an **auxiliary GLMM** on resident responses via
#'    [`fit_auxiliary_residents_glmm()`], optionally including random slopes
#'    `(0 + r_z || site)` and `(0 + C_z || site)` when
#'    `use_site_random_slopes = TRUE`.
#' 2. Convert GLMM coefficients to **sensitivities** (αᵢ, βᵢ, θᵢ/γᵢ) with
#'    [`derive_sensitivities()`], returning signed/unsigned variants and
#'    summary tests (Wald/LRT).
#' 3. If supported, extract **site-varying parameters** (αᵢₛ, Γᵢₛ) using
#'    [`site_varying_alpha_beta_gamma()`], along with decompositions into
#'    invader slopes and site deltas.
#'
#' **Output structure written to `fit$sensitivities`**
#' \itemize{
#'   \item `fit_coeffs`, `data_used`, `formula`: artefacts from the auxiliary GLMM.
#'   \item `alpha_i`, `beta_i`, `theta0`, `theta_i`, `gamma_i`,
#'         `alpha_signed_i`, `beta_signed_i`: estimated sensitivities.
#'   \item `wald_lrt`, `sens_df`, `clamp_summary`, `prior_note`: inference
#'         and diagnostics from [`derive_sensitivities()`].
#'   \item `site_alpha_beta_gamma`: full list returned by
#'         [`site_varying_alpha_beta_gamma()`] (if available).
#'   \item `alpha_is`, `Gamma_is`: site-varying crowding and abiotic matrices
#'         (if available).
#'   \item `site_alpha`, `site_gamma`: compact bundles with decompositions
#'         and a narrow data frame for inspection/joins.
#'   \item `a0..a2`, `b0..b2`: trait-plane slope/intercept components used by
#'         calibration steps downstream.
#'   \item `abg_df`: tidy table of site × invader effects when available.
#' }
#'
#' **Robustness notes**
#' - Performs finite-value checks on resident predictors (`r_js_z`, `C_js_z`,
#'   `S_js_z`).
#' - Resolves site identifiers from multiple possible locations to ensure
#'   compatibility with older and newer `assemble_matrices()` outputs.
#'
#' @return The input `fit` (invisibly) with an updated `fit$sensitivities`
#'   list (see Details).
#'
#' @seealso
#' - Input assembly: \href{https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.html}{`assemble_matrices()`}
#' - GLMM fit: \href{https://b-cubed-eu.github.io/invasimapr/reference/fit_auxiliary_residents_glmm.html}{`fit_auxiliary_residents_glmm()`}
#' - Sensitivity derivation: \href{https://b-cubed-eu.github.io/invasimapr/reference/derive_sensitivities.html}{`derive_sensitivities()`}
#' - Site-varying effects: \href{https://b-cubed-eu.github.io/invasimapr/reference/site_varying_alpha_beta_gamma.html}{`site_varying_alpha_beta_gamma()`}
#'
#' @examples
#' \dontrun{
#' fit <- prepare_inputs(sites = site_df, residents = resident_df,
#'                       invaders = invader_df, traits = trait_df)
#'
#' fit <- learn_sensitivities(fit, use_site_random_slopes = TRUE, lrt = TRUE)
#' names(fit$sensitivities)
#' head(fit$sensitivities$sens_df)
#'
#' # Example: map of site-varying Γ (if present)
#' if (!is.null(fit$sensitivities$Gamma_is)) {
#'   Gamma_is <- fit$sensitivities$Gamma_is
#'   # user mapping code goes here...
#' }
#' }
#'
#' @export
learn_sensitivities = function(fit,
                               use_site_random_slopes = TRUE,
                               lrt = TRUE) {

  # --- Preconditions -----------------------------------------------------------
  # Ensure we have the right container and sane resident predictor matrices.
  stopifnot(inherits(fit, "invasimapr_fit"))
  stopifnot(all(is.finite(fit$residents$r_js_z)),
            all(is.finite(fit$residents$C_js_z)),
            all(is.finite(fit$residents$S_js_z)))

  `%||%` = function(a, b) if (!is.null(a)) a else b

  # Robust site lookup (supports older/newer assemble_matrices outputs)
  sites = (fit$meta$sites %||% fit$sites %||% rownames(fit$inputs$comm_res))

  # --- 1) Auxiliary GLMM on residents -----------------------------------------
  # Includes optional site-level random slopes so that downstream extraction of
  # α_is and Γ_is is possible when variance components are non-trivial.
  aux = fit_auxiliary_residents_glmm(
    comm_res = fit$inputs$comm_res,
    r_js_z   = fit$residents$r_js_z,
    C_js_z   = fit$residents$C_js_z,
    S_js_z   = fit$residents$S_js_z,
    Q_res    = fit$traits$Q_res,
    use_site_random_slopes = use_site_random_slopes,
    verbose  = TRUE
  )

  # --- 2) Convert GLMM coefficients to sensitivities --------------------------
  # Produces α_i / β_i (signed & unsigned), θ_0 / θ_i (or γ_i), plus tests.
  sens = derive_sensitivities(
    fit_coeffs = aux$fit,
    Q_inv      = fit$traits$Q_inv,
    inv_ids    = rownames(fit$traits$Q_inv),
    lrt        = lrt
  )

  # --- 3) Optional site-varying effects (α_is, Γ_is) --------------------------
  # Quietly attempt extraction; not all models/datasets support these.
  abg = try(
    site_varying_alpha_beta_gamma(
      fit_coeffs = aux$fit,
      Q_inv      = fit$traits$Q_inv,
      sites      = sites,
      inv_ids    = rownames(fit$traits$Q_inv),
      quiet      = TRUE
    ),
    silent = TRUE
  )
  have_abg = !inherits(abg, "try-error")

  # Compact, communication-friendly bundles (subset + decompositions) ----------
  site_alpha_bc = if (have_abg) {
    list(
      alpha_is   = abg$alpha_is,
      slope_C_i  = abg$slope_C_i,
      delta_C_s  = abg$delta_C_s,
      df         = if ("df" %in% names(abg))
        abg$df[, intersect(c("site","invader","alpha_is","slope_C_i","delta_C_s"),
                           names(abg$df)), drop = FALSE]
      else NULL
    )
  } else NULL

  site_gamma_bc = if (have_abg) {
    list(
      Gamma_is   = abg$Gamma_is,
      theta_i    = abg$theta_i,
      delta_r_s  = abg$delta_r_s,
      df         = if ("df" %in% names(abg))
        abg$df[, intersect(c("site","invader","Gamma_is","theta_i","delta_r_s"),
                           names(abg$df)), drop = FALSE]
      else NULL
    )
  } else NULL

  # --- 4) Pack results into `fit$sensitivities` -------------------------------
  fit$sensitivities = list(
    # Model artefacts
    fit_coeffs     = aux$fit,
    data_used      = aux$data,
    formula        = aux$formula,

    # Global/trait-varying sensitivities
    alpha_i        = sens$alpha_i,
    alpha_signed_i = sens$alpha_signed_i,
    beta_i         = sens$beta_i,
    beta_signed_i  = sens$beta_signed_i,
    theta0         = sens$theta0,
    theta_i        = sens$theta_i,
    gamma_i        = sens$gamma_i,

    # Inference and diagnostics
    wald_lrt       = sens$lrt,
    sens_df        = sens$df,
    clamp_summary  = sens$clamp_summary,
    prior_note     = sens$prior_note,

    # Site-varying effects (if available)
    site_alpha_beta_gamma = if (have_abg) abg else NULL,
    alpha_is       = if (have_abg) abg$alpha_is else NULL,
    Gamma_is       = if (have_abg) abg$Gamma_is else NULL,
    site_alpha     = site_alpha_bc,
    site_gamma     = site_gamma_bc,

    # Trait-plane components (used by calibration elsewhere)
    a0 = sens$a0, a1 = sens$a1, a2 = sens$a2,
    b0 = sens$b0, b1 = sens$b1, b2 = sens$b2,

    # Tidy table for joins/plots
    abg_df = if (have_abg) abg$df else NULL
  )

  # Return updated fit for chaining
  fit
}
