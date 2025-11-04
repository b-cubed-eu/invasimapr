# ======================================================================
# 5) LEARN SENSITIVITIES  (minor robustness update)
# ======================================================================

#' Learn sensitivities (αᵢ, βᵢ, θᵢ/γᵢ) and optional site-varying αᵢₛ, Γᵢₛ
#' @export
learn_sensitivities = function(fit,
                                use_site_random_slopes = TRUE,
                                lrt = TRUE) {

  # (loaders for aux functions omitted for brevity; keep your originals)
  stopifnot(inherits(fit, "invasimapr_fit"))

  stopifnot(all(is.finite(fit$residents$r_js_z)),
            all(is.finite(fit$residents$C_js_z)),
            all(is.finite(fit$residents$S_js_z)))


  `%||%` = function(a, b) if (!is.null(a)) a else b

  # Robust site lookup (new assemble_matrices returns $sites at top-level)
  sites = (fit$meta$sites %||% fit$sites %||% rownames(fit$inputs$comm_res))

  aux = fit_auxiliary_residents_glmm(
    comm_res = fit$inputs$comm_res,
    r_js_z   = fit$residents$r_js_z,
    C_js_z   = fit$residents$C_js_z,
    S_js_z   = fit$residents$S_js_z,
    Q_res    = fit$traits$Q_res,
    use_site_random_slopes = use_site_random_slopes,
    verbose  = TRUE
  )

  sens = derive_sensitivities(
    fit_coeffs = aux$fit,
    Q_inv      = fit$traits$Q_inv,
    inv_ids    = rownames(fit$traits$Q_inv),
    lrt        = lrt
  )

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

  site_alpha_bc = if (have_abg) {
    list(
      alpha_is   = abg$alpha_is,
      slope_C_i  = abg$slope_C_i,
      delta_C_s  = abg$delta_C_s,
      df         = if ("df" %in% names(abg))
        abg$df[, intersect(c("site","invader","alpha_is","slope_C_i","delta_C_s"), names(abg$df)), drop = FALSE]
      else NULL
    )
  } else NULL

  site_gamma_bc = if (have_abg) {
    list(
      Gamma_is   = abg$Gamma_is,
      theta_i    = abg$theta_i,
      delta_r_s  = abg$delta_r_s,
      df         = if ("df" %in% names(abg))
        abg$df[, intersect(c("site","invader","Gamma_is","theta_i","delta_r_s"), names(abg$df)), drop = FALSE]
      else NULL
    )
  } else NULL

  fit$sensitivities = list(
    fit_coeffs     = aux$fit,
    data_used      = aux$data,
    formula        = aux$formula,

    alpha_i        = sens$alpha_i,
    alpha_signed_i = sens$alpha_signed_i,
    beta_i         = sens$beta_i,
    beta_signed_i  = sens$beta_signed_i,
    theta0         = sens$theta0,
    theta_i        = sens$theta_i,
    gamma_i        = sens$gamma_i,
    wald_lrt       = sens$lrt,
    sens_df        = sens$df,
    clamp_summary  = sens$clamp_summary,
    prior_note     = sens$prior_note,

    site_alpha_beta_gamma = if (have_abg) abg else NULL,

    alpha_is       = if (have_abg) abg$alpha_is else NULL,
    Gamma_is       = if (have_abg) abg$Gamma_is else NULL,
    site_alpha     = site_alpha_bc,
    site_gamma     = site_gamma_bc,

    a0 = sens$a0, a1 = sens$a1, a2 = sens$a2,
    b0 = sens$b0, b1 = sens$b1, b2 = sens$b2,

    abg_df = if (have_abg) abg$df else NULL
  )

  fit
}
