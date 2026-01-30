#' Learn sensitivities (alpha_i, beta_i, theta_i/gamma_i) and optional site-varying alpha_is, Gamma_is
#'
#' Fits an auxiliary GLMM on **resident** data to estimate invader-level
#' sensitivities to crowding and saturation (alpha_i, beta_i), and abiotic conversion
#' slopes (theta_i or gamma_i), with optional **site-varying** random slopes that yield
#' per-site adjustments (alpha_is, Gamma_is). Results are written into the
#' `fit$sensitivities` slot of an [`invasimapr_fit`] object for downstream
#' invasion-fitness and establishment steps.
#'
#' @param fit An object of class `invasimapr_fit` produced by
#'   [`prepare_inputs()`] / [`assemble_matrices()`], containing resident
#'   matrices (`r_js_z`, `C_js_z`, `S_js_z`), trait-space structures (`Q_res`,
#'   `Q_inv`), and resident community layout (`inputs$comm_res`).
#' @param use_site_random_slopes Logical; if `TRUE`, the auxiliary model is fit
#'   with site-level random slopes for the abiotic and crowding terms, enabling
#'   estimation of site-varying alpha_is and Gamma_is when supported by the data.
#'   Defaults to `TRUE`.
#' @param lrt Logical; if `TRUE`, compute Wald/LRT summaries for key contrasts
#'   (e.g., trait-varying vs global slopes) to guide model c#' Learn sensitivities (alpha_i, beta_i, theta_i / gamma_i) and optional site effects
#'
#' @description
#' Fits an auxiliary GLMM on resident data to estimate invader-level sensitivities
#' to crowding and saturation (\eqn{\alpha_i}, \eqn{\beta_i}), and abiotic conversion
#' slopes (\eqn{\theta_i} or \eqn{\gamma_i}). When supported by the data, optional
#' site-level random slopes yield site-varying adjustments
#' (\eqn{\alpha_{is}}, \eqn{\Gamma_{is}}).
#'
#' Results are written into the `fit$sensitivities` slot of an
#' \link{new_invasimapr_fit} object for downstream invasion-fitness and
#' establishment calculations.
#'
#' @param fit An object produced by \link{prepare_inputs} and
#'   \link{assemble_matrices}, containing resident predictor matrices
#'   (`r_js_z`, `C_js_z`, `S_js_z`), trait-space structures (`Q_res`, `Q_inv`),
#'   and the resident community layout (`inputs$comm_res`).
#' @param use_site_random_slopes Logical; if `TRUE`, the auxiliary model includes
#'   site-level random slopes for abiotic and crowding terms, enabling estimation
#'   of site-varying \eqn{\alpha_{is}} and \eqn{\Gamma_{is}} when supported by the data.
#'   Defaults to `TRUE`.
#' @param lrt Logical; if `TRUE`, compute Wald or likelihood-ratio tests for key
#'   contrasts (e.g., trait-varying versus global slopes). Defaults to `TRUE`.
#'
#' @details
#' \strong{Workflow}
#' \enumerate{
#'   \item Fit an auxiliary GLMM on resident responses using
#'         \link{fit_auxiliary_residents_glmm}, optionally including site-level
#'         random slopes for \eqn{r_z} and \eqn{C_z}.
#'   \item Convert GLMM coefficients to sensitivities
#'         (\eqn{\alpha_i}, \eqn{\beta_i}, \eqn{\theta_i} or \eqn{\gamma_i}) using
#'         \link{derive_sensitivities}, returning signed and unsigned variants
#'         plus inference summaries.
#'   \item When supported, extract site-varying effects
#'         (\eqn{\alpha_{is}}, \eqn{\Gamma_{is}}) via
#'         \link{site_varying_alpha_beta_gamma}.
#' }
#'
#' The resulting components are stored in `fit$sensitivities`, including:
#' \itemize{
#'   \item global and trait-varying sensitivities;
#'   \item inference diagnostics and clamping summaries;
#'   \item optional site-varying matrices and compact decomposition tables.
#' }
#'
#' @return The input `fit` object (invisibly), with an updated
#'   `fit$sensitivities` list.
#'
#' @seealso
#' \link{prepare_inputs},
#' \link{assemble_matrices},
#' \link{fit_auxiliary_residents_glmm},
#' \link{derive_sensitivities},
#' \link{site_varying_alpha_beta_gamma}
#'
#' @examples
#' \dontrun{
#' fit <- prepare_inputs(
#'   sites = site_df,
#'   residents = resident_df,
#'   invaders = invader_df,
#'   traits = trait_df
#' )
#'
#' fit <- learn_sensitivities(fit, use_site_random_slopes = TRUE)
#' names(fit$sensitivities)
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
  # alpha_is and Gamma_is are possible when variance components are non-trivial.
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
  # Produces alpha_i / beta_i (signed & unsigned), theta_0 / theta_i (or gamma_i), plus tests.
  sens = derive_sensitivities(
    fit_coeffs = aux$fit,
    Q_inv      = fit$traits$Q_inv,
    inv_ids    = rownames(fit$traits$Q_inv),
    lrt        = lrt
  )

  # --- 3) Optional site-varying effects (alpha_is, Gamma_is) --------------------------
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
